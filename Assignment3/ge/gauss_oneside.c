#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

int main(int argc, char** argv) {
	//	A * x = b
	//	Matrix Names
	char matrix_name[200], vector_name[200], solution_name[200];
	//	Global rows and columns of A, # processes and rank
	int rows, columns, size, rank;
	//	**pointer to matrix A, *pointer to matrix A in vector form, respective *pointers to rhs b and solution x
	double **matrix_2d_mapped, *matrix_1D_mapped, *rhs, *solution;
	//	Times to measure
	double total_time, io_time = 0, setup_time, kernel_time, mpi_time = 0;
	double total_start, io_start, setup_start, kernel_start, mpi_start;
	//	Files to read
	FILE *matrix_file, *vector_file, *solution_file;
	MPI_Status status;     
	MPI_Request req;
	MPI_Status m_status[2];
	MPI_Request req_rec[2];
        MPI_Request req_send[2];

	if( argc != 2 ) { 
		perror("The base name of the input matrix and vector files must be given\n"); 
		exit(-1);
	}

	int print_a = 0;
	int print_b = 0;
	int print_x = 0;

	sprintf(matrix_name,   "%s.mat", argv[1]);
	sprintf(vector_name,   "%s.vec", argv[1]);
	sprintf(solution_name, "%s.sol", argv[1]);

	//	Basic MPI initialization
	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if(rank == 0){
		printf("Solving the Ax=b system with Gaussian Elimination:\n");
		printf("READ:  Matrix   (A) file name: \"%s\"\n", matrix_name);
		printf("READ:  RHS      (b) file name: \"%s\"\n", vector_name);
		printf("WRITE: Solution (x) file name: \"%s\"\n", solution_name);
	}

	total_start = MPI_Wtime();

	//	Index for loops in next if-statement
	int row, column, index;
	//	Some memory allocations, reads and basic assertions
	if(rank == 0) {
		io_start = MPI_Wtime();
		if ((matrix_file = fopen (matrix_name, "r")) == NULL) {
			perror("Could not open the specified matrix file");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		fscanf(matrix_file, "%d %d", &rows, &columns);     
		if(rows != columns) {
			perror("Only square matrices are allowed\n");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}  	
		if(rows % size != 0) {
			perror("The matrix should be divisible by the number of processes\n");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}  	

		matrix_2d_mapped = (double **) malloc(rows * sizeof(double *));
		for(row = 0; row < rows; row++){
			matrix_2d_mapped[row] = (double *) malloc(rows * sizeof(double));
			for(column = 0; column < columns; column++){
				fscanf(matrix_file, "%lf", &matrix_2d_mapped[row][column]);
			}
		}
		fclose(matrix_file);

		if ((vector_file = fopen (vector_name, "r")) == NULL){
			perror("Could not open the specified vector file");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		int rhs_rows;
		fscanf(vector_file, "%d", &rhs_rows);     
		if(rhs_rows != rows){
			perror("RHS rows must match the sizes of A");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		rhs  = (double *)malloc(rows * sizeof(double));
		for (row = 0; row < rows; row++){
			fscanf(vector_file, "%lf",&rhs[row]);
		}
		fclose(vector_file); 
		io_time += MPI_Wtime() - io_start;

		matrix_1D_mapped = (double *) malloc(rows * rows * sizeof(double));
		index = 0;
		for(row=0; row < rows; row++){
			for(column=0; column < columns; column++){
				matrix_1D_mapped[index++] = matrix_2d_mapped[row][column];
			}
		}
		solution = (double *) malloc (rows * sizeof(double));
	}

	setup_start = MPI_Wtime();

        // --- PREPARE WINDOS FOR COMMUNICATION ---
        MPI_Win win1;
        MPI_Win win2;
        //MPI_Win_create(&rows, 1*sizeof(int), 0, MPI_INFO_NULL, MPI_COMM_WORLD, &win1 );
        //MPI_Win_create(&columns, 1*sizeof(int), 0, MPI_INFO_NULL, MPI_COMM_WORLD, &win2 );


	int i;
	//	Rank 0 sends number of rows and columns and the other processes receive them
	if(rank == 0) {
		for(i = 1; i < size; i++){
		MPI_Isend(&rows, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &req_send[0]);
                MPI_Isend(&columns, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &req_send[1]);
		MPI_Request_free(&req_send[0]);
MPI_Request_free(&req_send[1]);
		}
	} else {
		MPI_Irecv(&rows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &req_rec[0]);
		MPI_Irecv(&columns, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &req_rec[1]);
MPI_Waitall(2,req_rec,m_status);
	}	// overlap communication with different buffers, do row transmission while doing columns one

	//	Matrix A will be divided in groups of rows(chunks) and each group sent to each process
	//	The same is done for the vectors x and b

	//	local_block_size is the number of rows each process is going to work with
	int local_block_size = rows / size;

	int process, column_pivot;
	double tmp, pivot;
	//	contains all entries of its chunk of the matrix A
	double *matrix_local_block = (double *) malloc(local_block_size * rows * sizeof(double));
	//	contains the chunk of the rhs b
	double *rhs_local_block = (double *) malloc(local_block_size * sizeof(double));
	//	(to be sent/received) contains in order: rank, chunk of A, chunk of rhs
	double *pivots = (double *) malloc((local_block_size + (rows * local_block_size) + 1) * sizeof(double));
	//	contains normalized rhs b after one step of GE (divided by pivot)
	double *local_work_buffer = (double *) malloc(local_block_size * sizeof(double));
	//	(to be sent/received) contains in order: position in rows, rhs b in such positions (for example [3,4,b_3,b_4])
	double *accumulation_buffer = (double *) malloc(local_block_size * 2 * sizeof(double));
	//	contains chunk of already computed solution x
	double *solution_local_block = (double *) malloc(local_block_size * sizeof(double));
	
	
	// --- PUT DATA INTO WINDOWS ---
	MPI_Win_allocate(rows * rows * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &matrix_1D_mapped, &win1);
	MPI_Win_fence(0,win1);
	MPI_Win_allocate(rows * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &rhs, &win2);
	MPI_Win_fence(0,win2);

	//	send/receive respective chunk of data of A and rhs b to each process
	if(rank == 0) {
		for(i = 1; i < size; i++){
			//MPI_Isend((matrix_1D_mapped + (i * (local_block_size * rows))), (local_block_size * rows), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &req);
			//MPI_Isend((rhs + (i * local_block_size)), local_block_size, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &req);
			//MPI_Wait(&req,&status);
		}
		for(i = 0; i < local_block_size * rows; i++){
			matrix_local_block[i] = matrix_1D_mapped[i];
		}
		for(i = 0; i < local_block_size; i++){
			rhs_local_block[i] = rhs[i];
		}
	} else {
		//MPI_Irecv(matrix_local_block, local_block_size * rows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &req);
		//MPI_Irecv(rhs_local_block, local_block_size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &req);
		//MPI_Wait(&req,&status);
		
		// --- RECEIVE FROM SHARED ---
		MPI_Get(matrix_local_block, local_block_size * rows, MPI_DOUBLE, 0, 0, local_block_size * rows, MPI_DOUBLE, win1 );
		MPI_Win_fence(0,win1);
		MPI_Get(rhs_local_block, local_block_size, MPI_DOUBLE, 0, 0, local_block_size, MPI_DOUBLE, win2 );
		MPI_Win_fence(0,win2);
		
		
	}// Here Irecv for the two recv's since allocating data takes time but they are unrelated!

	setup_time = MPI_Wtime() - setup_start;
	kernel_start = MPI_Wtime();

	//	receive *pivots from previous ranks, make its chunk of A upper triangular and recompute rhs b
	for(process = 0; process < rank; process++) {
		mpi_start = MPI_Wtime();
		
		// --- GET PIVOTS ---
		//MPI_Recv(pivots, (local_block_size * rows + local_block_size + 1), MPI_DOUBLE, process, process, MPI_COMM_WORLD, &status);
		MPI_Get(pivots, (local_block_size * rows + local_block_size + 1), MPI_DOUBLE, process,  process, (local_block_size * rows + local_block_size + 1), MPI_DOUBLE, win1 );
		MPI_Win_fence(0,win1);
		
		mpi_time += MPI_Wtime() - mpi_start;
                
		for(row = 0; row < local_block_size; row++){
			column_pivot = ((int) pivots[0]) * local_block_size + row;
			for (i = 0; i < local_block_size; i++){
				index = i * rows;
				tmp = matrix_local_block[index + column_pivot];
				for (column = column_pivot; column < columns; column++){
					matrix_local_block[index + column] -=  tmp * pivots[(row * rows) + (column + local_block_size + 1)];
				}
				rhs_local_block[i] -= tmp * pivots[row + 1];
				matrix_local_block[index + column_pivot] = 0.0;
			}
		}
	}

	//	performs GE for its chunk of A and rhs b
	for(row = 0; row < local_block_size; row++){
		column_pivot = (rank * local_block_size) + row;
		index = row * rows;
		pivot = matrix_local_block[index + column_pivot];
		assert(pivot!= 0);

		for (column = column_pivot; column < columns; column++){
			matrix_local_block[index + column] = matrix_local_block[index + column]/pivot; 
			pivots[index + column + local_block_size + 1] = matrix_local_block[index + column];
		}

		local_work_buffer[row] = (rhs_local_block[row])/pivot;
		pivots[row+1] =  local_work_buffer[row];

		for (i = (row + 1); i < local_block_size; i++) {
			tmp = matrix_local_block[i*rows + column_pivot];
			for (column = column_pivot+1; column < columns; column++){
				matrix_local_block[i*rows+column] -=  tmp * pivots[index + column + local_block_size + 1];
			}
			rhs_local_block[i] -= tmp * local_work_buffer[row];
			matrix_local_block[i * rows + row] = 0;
		}
	}

	//	send *pivots
	for (process = (rank + 1); process < size; process++) {
		pivots[0] = (double) rank;
		mpi_start = MPI_Wtime();
		
		// --- ALLOCATE PIVOTS ---
		//MPI_Send( pivots, (local_block_size * rows + local_block_size + 1), MPI_DOUBLE, process, rank, MPI_COMM_WORLD);
		MPI_Win_allocate((local_block_size * rows + local_block_size + 1) * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &pivots, &win1);
		MPI_Win_fence(0,win1);
		
		
		mpi_time += MPI_Wtime() - mpi_start;
	} 

	//	receive chunks of rhs b after GE
	for (process = (rank + 1); process<size; ++process) {
		mpi_start = MPI_Wtime();
		
		// --- GET RHS CHUNKS ---
		//MPI_Recv( accumulation_buffer, (2 * local_block_size), MPI_DOUBLE, process, process, MPI_COMM_WORLD, &status); 
		MPI_Get( accumulation_buffer, (2 * local_block_size), MPI_DOUBLE, process, process, (2 * local_block_size), MPI_DOUBLE, win1 );
		MPI_Win_fence(0,win1);
		
		mpi_time += MPI_Wtime() - mpi_start;

		for (row  = (local_block_size - 1); row >= 0; row--) {
			for (column = (local_block_size - 1);column >= 0; column--) {
				index = (int) accumulation_buffer[column];
				local_work_buffer[row] -= accumulation_buffer[local_block_size + column] * matrix_local_block[row * rows + index];
			}
		}
	}

	//	compute local solutions(chunk of x)
	for (row = (local_block_size - 1); row >= 0; row--) {
		index = rank * local_block_size + row;
		accumulation_buffer[row] = (double) index;
		accumulation_buffer[local_block_size+row] = solution_local_block[row] = local_work_buffer[row];
		for (i = (row - 1); i >= 0; i--){
			local_work_buffer[i] -= solution_local_block[row] * matrix_local_block[ (i * rows) + index];
		}
	}

	//	send chunks of rhs b after GE to other ranks
	for (process = 0; process < rank; process++){
		mpi_start = MPI_Wtime();
		
		// --- ALLOC RHS CHUNKS ---
		//MPI_Send( accumulation_buffer, (2 * local_block_size), MPI_DOUBLE, process, rank, MPI_COMM_WORLD);
		MPI_Win_allocate((2 * local_block_size) * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &accumulation_buffer, &win1);
		MPI_Win_fence(0,win1); 
		
		mpi_time += MPI_Wtime() - mpi_start;
	}

	//	send/receive solutions
	// --- ALLOC AND GET ---
	MPI_Win_fence(0,win1); 
	if(rank == 0) {
		for(i = 0; i < local_block_size; i++){
			solution[i] = solution_local_block[i];
		}
		mpi_start = MPI_Wtime();
		for(i = 1; i < size; i++){
			//MPI_Recv(solution + (i * local_block_size), local_block_size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
			MPI_Get( solution + (i * local_block_size), local_block_size, MPI_DOUBLE, i, 0,  local_block_size, MPI_DOUBLE, win1 );
			
		}
		mpi_time += MPI_Wtime() - mpi_start;
	} else {
		mpi_start = MPI_Wtime();
		//MPI_Send(solution_local_block, local_block_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MPI_Win_allocate(local_block_size * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &solution_local_block, &win1);
		mpi_time += MPI_Wtime() - mpi_start;
	}
        MPI_Win_fence(0,win1); 
        
	kernel_time = MPI_Wtime() - kernel_start;

	//	basic assertions and results
	if (rank == 0) {
		io_start = MPI_Wtime();
		if ((solution_file = fopen(solution_name, "w+")) == NULL) {
			perror("Could not open the solution file");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		fprintf(solution_file, "%d\n", rows);
		for(i = 0; i < rows; i++) {
			fprintf(solution_file, "%f ", solution[i]);
		}
		fprintf(solution_file, "\n");
		fclose(solution_file);
		io_time += MPI_Wtime() - io_start;

		if(print_a){
			printf("\nSystem Matrix (A):\n");
			for (row = 0; row < rows; row++) {
				for (column = 0; column < columns; column++){
					printf("%4.1f ", matrix_2d_mapped[row][column]);
				}
				printf("\n");
			}
		}

		if(print_b){
			printf("\nRHS Vector (b):\n");
			for (row = 0; row < rows; row++) {
				printf("%4.1f\n", rhs[row]);
			}
		}

		if(print_x){
			printf("\n\nSolution Vector (x):\n");
			for(row = 0; row < rows; row++){
				printf("%4.4f\n",solution[row]);
			}
		}
	}

	total_time = MPI_Wtime() - total_start;

	printf("[R%02d] Times: IO: %f; Setup: %f; Compute: %f; MPI: %f; Total: %f;\n", 
			rank, io_time, setup_time, kernel_time, mpi_time, total_time);

	if(rank == 0){
		for(i = 0; i < rows; i++){
			free(matrix_2d_mapped[i]);
		}
		free(matrix_2d_mapped);
		free(rhs);
		free(solution);
	}
	free(matrix_local_block);
	free(rhs_local_block);
	free(pivots);
	free(local_work_buffer);
	free(accumulation_buffer);
	free(solution_local_block);

	MPI_Finalize(); 
	return 0;
}

