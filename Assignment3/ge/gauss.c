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

	MPI_Request req;
	MPI_Status status;
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

	MPI_Group comm_group;
	MPI_Comm_group(MPI_COMM_WORLD, &comm_group);

	if(rank == 0){
		printf("Solving the Ax=b system with Gaussian Elimination:\n");
		printf("READ:  Matrix   (A) file name: \"%s\"\n", matrix_name);
		printf("READ:  RHS      (b) file name: \"%s\"\n", vector_name);
		printf("WRITE: Solution (x) file name: \"%s\"\n", solution_name);
	}


	total_start = MPI_Wtime();

	//	Index for loops in next if-statement
	int row, column, index, i, ii;
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
			printf("rows = %d and size = %d", rows, size);
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

	MPI_Win test_win; //Window for sending #rows #columns
	int dim[2];
	//	Rank 0 sends number of rows and columns and the other processes receive them
	if(rank == 0) {
		dim[0] = rows; dim[1]=columns;
		MPI_Win_create(dim, 2*sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &test_win);
		MPI_Win_fence(0, test_win);
		//Wait for message to be read with GET from the other ranks!
		MPI_Win_fence(0, test_win);
	} else {
		MPI_Win_create(NULL, 0, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &test_win);
		MPI_Win_fence(0, test_win);
		MPI_Get(dim,2, MPI_INT,0,0,2,MPI_INT,test_win);
		MPI_Win_fence(0, test_win);
		rows = dim[0];
		columns = dim[1];
		//MPI_Irecv(&rows, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &req_rec[0]);
		//MPI_Irecv(&columns, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &req_rec[1]);
		//MPI_Waitall(2,req_rec,m_status);
	}	// overlap communication with different buffers, do row transmission while doing columns one
	//MPI_Win_free(&test_win);

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


	MPI_Win test_win_chunks; //Window for sending chunk of data of A and rhs b to each process
	if(rank == 0) {
		int big_dim = rows * rows + rows;
		double *matrix_1D_mapped_rhs = (double *) malloc(big_dim * sizeof(double));
		int counter = 0;
		//Making new big array that ill be window
		for(i=0; i < (rows * rows); i++){
			matrix_1D_mapped_rhs[counter++] = matrix_1D_mapped[i];
		}
		for(i=0; i < rows; i++){
			matrix_1D_mapped_rhs[counter++] = rhs[i];
		}
		//Copying local part of process 0
		for(i = 0; i < local_block_size * rows; i++){
			matrix_local_block[i] = matrix_1D_mapped[i];
		}
		for(i = 0; i < local_block_size; i++){
			rhs_local_block[i] = rhs[i];
		}
		MPI_Win_create(matrix_1D_mapped_rhs, big_dim * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &test_win_chunks);
		MPI_Win_fence(0, test_win_chunks);
		//Wait for message to be read with GET from the other ranks!
		MPI_Win_fence(0, test_win_chunks);
		free(matrix_1D_mapped_rhs);
	} else {
		MPI_Win_create(NULL, 0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &test_win_chunks);
		MPI_Win_fence(0, test_win_chunks);
		int target_disp = rank * (local_block_size * rows);
		int target_disp_rhs = rows * rows + rank * local_block_size ;
		int count = local_block_size * rows;
		MPI_Get(matrix_local_block, count, MPI_DOUBLE, 0, target_disp, count, MPI_DOUBLE, test_win_chunks);
		MPI_Get(rhs_local_block, local_block_size, MPI_DOUBLE, 0, target_disp_rhs, local_block_size, MPI_DOUBLE, test_win_chunks);
		MPI_Win_fence(0, test_win_chunks);
	}
	//MPI_Win_free(&test_win_chunks);

	//TIME CALCULATION
	setup_time = MPI_Wtime() - setup_start;
	kernel_start = MPI_Wtime();

//SOLVING - GE PERFORMING
	//Receive *pivots from previous ranks, and update its chunk of A and rhs b with respect to the other chunks
	MPI_Win test_win_ge; //Window for receiving pivots from processes before
	int size_pivots = local_block_size + rows * local_block_size + 1;
	if(rank==0){
		MPI_Win_create(NULL, 0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &test_win_ge);
		//printf("\nInside rank %d in part SOLVING - GE PERFORMING\n", rank);
	}
	else
    {
		MPI_Win_create(pivots, size_pivots * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &test_win_ge);
		//printf("\nInside rank %d in part Win Creation - SOLVING - GE PERFORMING\n", rank);
    }
	//MPI_Group *all_groups = (MPI_Group *)malloc(sizeof(MPI_Group) * (size));
	//MPI_Group all_groups2[size];

	//1 - First step in this subsection - calculation based on pivots from processes above
	for(process = 0; process < rank; process++) {
		/*int ranks_in_group_to_get_from[size-process-1];
		//int *ranks_in_group_to_get_from = (int *) malloc(sizeof(int) * (size-process-1));
		int ii = 0;
		int process_index;
		for (process_index = process; process_index < rank; process_index++){
			ranks_in_group_to_get_from[ii++]=process_index;
		}
		for (process_index = rank+1; process_index < size; process_index++){
			ranks_in_group_to_get_from[ii++]=process_index;
		}*/
		mpi_start = MPI_Wtime();
		//MPI_Group_incl(comm_group, size-process-1, ranks_in_group_to_get_from, all_groups2+process);
		MPI_Group any_group;
		MPI_Group_incl(comm_group, 1, &process, &any_group);
		MPI_Win_post(any_group, 0, test_win_ge);
		MPI_Win_wait(test_win_ge);
		//printf("\nInside rank %d in part THAY ALL TAKE PIVOTS FROM ME\n", rank);
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
		MPI_Group_free(&any_group);
		//free(ranks_in_group_to_get_from);
	}

	//2 - Second step in this subsection - calculation based own pivots
	for(row = 0; row < local_block_size; row++){
		column_pivot = (rank * local_block_size) + row;
		index = row * rows;
		pivot = matrix_local_block[index + column_pivot];
		assert(pivot!= 0);

		for (column = column_pivot; column < columns; column++){
			matrix_local_block[index + column] = matrix_local_block[index + column]/pivot; 
			// I change/overwrite last pivots I got and this will be the one I'm sending....
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
	pivots[0] = (double) rank;


	//3 - Third step in this subsection - sending pivots to the processe below
	//int ranks_in_group_to_put[size-rank-1];
	int *ranks_in_group_to_put = (int *) malloc(sizeof(int) * (size-rank-1));
	ii = 0;
	for (process = (rank + 1); process < size; process++){
		ranks_in_group_to_put[ii++]=process;
	}
	mpi_start = MPI_Wtime();
	//MPI_Group_incl(comm_group, size-rank-1, ranks_in_group_to_put, all_groups2+rank);
	//MPI_Win_start(all_groups2[rank], 0, test_win_ge);
	MPI_Group anny_group2;
	MPI_Group_incl(comm_group, size-rank-1, ranks_in_group_to_put, &anny_group2);
	MPI_Win_start(anny_group2, 0, test_win_ge);
	for (process = (rank + 1); process < size; process++) {
		MPI_Put(pivots, size_pivots, MPI_DOUBLE, process, 0, size_pivots, MPI_DOUBLE, test_win_ge);
	} 
	MPI_Win_complete(test_win_ge);
	//MPI_Win_free(&test_win_ge);
	MPI_Group_free(&anny_group2);
	mpi_time += MPI_Wtime() - mpi_start;
	free(ranks_in_group_to_put);

//SOLVING - BACK SUBSTITUTION
	MPI_Win test_win_bs; //Window for receiving pivots from processes before
	int size_accumulation_buffer = 2 * local_block_size;
	if (rank != size-1){
		MPI_Win_create(accumulation_buffer, size_accumulation_buffer * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &test_win_bs);
	}else{
		MPI_Win_create(NULL, 0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &test_win_bs);
	}
	//Receive chunks of rhs b after GE. Now vice verse I'm receiving from all the processes beneath me
	//for (process = (rank + 1); process<size; ++process)
	for (process = size-1; process>rank; process--) {
		mpi_start = MPI_Wtime();
		MPI_Group any_group3;
		MPI_Group_incl(comm_group, 1, &process, &any_group3);
		MPI_Win_post(any_group3, 0, test_win_bs);
		MPI_Win_wait(test_win_bs);
		mpi_time += MPI_Wtime() - mpi_start;

		for (row  = (local_block_size - 1); row >= 0; row--) {
			for (column = (local_block_size - 1);column >= 0; column--) {
				index = (int) accumulation_buffer[column];
				local_work_buffer[row] -= accumulation_buffer[local_block_size + column] * matrix_local_block[row * rows + index];
			}
		}
		MPI_Group_free(&any_group3);
	}

	//	compute local solutions(chunk of x) First for zero process!!!
	for (row = (local_block_size - 1); row >= 0; row--) {
		index = rank * local_block_size + row;
		accumulation_buffer[row] = (double) index;
		accumulation_buffer[local_block_size+row] = solution_local_block[row] = local_work_buffer[row];
		for (i = (row - 1); i >= 0; i--){
			local_work_buffer[i] -= solution_local_block[row] * matrix_local_block[ (i * rows) + index];
		}
	}


	//	send chunks of rhs b after GE to other ranks, Now vice verse I'm sending to all the processes above me
	int *ranks_in_group_to_put_assbuf = (int *) malloc(sizeof(int) * (rank));
	ii = 0;
	for (process = 0; process < rank; process++){
		ranks_in_group_to_put_assbuf[ii++]=process;
	}
	/*for (process = rank-1; process > 0; process--){
			ranks_in_group_to_put_assbuf[ii++]=process;
	}*/
	mpi_start = MPI_Wtime();
	MPI_Group anny_group4;
	MPI_Group_incl(comm_group, rank, ranks_in_group_to_put_assbuf, &anny_group4);
    MPI_Win_start(anny_group4, 0, test_win_bs);
	for (process = 0; process < rank; process++){
		MPI_Put(accumulation_buffer, size_accumulation_buffer, MPI_DOUBLE, process, 0, size_accumulation_buffer, MPI_DOUBLE, test_win_bs);
	}
	MPI_Win_complete(test_win_bs);
	//MPI_Win_free(&test_win_bs);
	mpi_time += MPI_Wtime() - mpi_start;
	MPI_Group_free(&anny_group4);
	free(ranks_in_group_to_put_assbuf);


//Gathering the solution
	MPI_Win test_win_soultion;
	if(rank == 0) {
		mpi_start = MPI_Wtime();
		MPI_Win_create(solution, rows * sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &test_win_soultion);
		MPI_Win_fence(0, test_win_soultion);
		//Wait for solution to be put with PUT from the other ranks!
	    MPI_Win_fence(0, test_win_soultion);
		mpi_time += MPI_Wtime() - mpi_start;
		for(i = 0; i < local_block_size; i++){
			solution[i] = solution_local_block[i];
		}
	} else {
		mpi_start = MPI_Wtime();
		MPI_Win_create(NULL, 0, sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &test_win_soultion);
		MPI_Win_fence(0, test_win_soultion);
		int target_disp_solutio  = rank * local_block_size ;
		MPI_Put(solution_local_block, local_block_size, MPI_DOUBLE, 0, target_disp_solutio, local_block_size, MPI_DOUBLE, test_win_soultion);
		MPI_Win_fence(0, test_win_soultion);
		mpi_time += MPI_Wtime() - mpi_start;
	}
	//MPI_Win_free(&test_win_soultion);

	//TIME CALCULATION
	kernel_time = MPI_Wtime() - kernel_start;

	//Basic assertions and results
	if (rank == 0) {
		io_start = MPI_Wtime();
		if ((solution_file = fopen(solution_name, "w+")) == NULL) {
			perror("Could not open the solution file");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
		//MPI_Wait(&req, &status);
		//MPI_Waitall(size-1,many_req,many_status);
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

	//free allocated requests and statuses
	//free(many_req);
	//free(many_status);

	MPI_Win_free(&test_win);
	MPI_Win_free(&test_win_chunks);
	MPI_Win_free(&test_win_ge); //Window for receiving pivots from processes before
	MPI_Win_free(&test_win_bs);
	MPI_Win_free(&test_win_soultion);
   // for (process = 0; process < size; process++){
    	//MPI_Group_free(all_groups+process);
    	//MPI_Group_free(all_groups2+process);
    //}
	MPI_Group_free(&comm_group);

	MPI_Finalize(); 
	return 0;
}

