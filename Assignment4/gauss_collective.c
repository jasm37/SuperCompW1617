#include <assert.h>
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
	
	MPI_Group cur_group;
	MPI_Comm cur_comm;
	MPI_Comm_group(MPI_COMM_WORLD, &cur_group);
	MPI_Comm_create(MPI_COMM_WORLD, cur_group, &cur_comm);
	int zerorank = 0;

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

	int i;
	//	Rank 0 sends number of rows and columns and the other processes receive them
	MPI_Bcast(&rows, 1, MPI_INT, 0 , MPI_COMM_WORLD);
	MPI_Bcast(&columns, 1, MPI_INT, 0 , MPI_COMM_WORLD);

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

	//	send/receive respective chunk of data of A and rhs b to each process

	MPI_Scatter( matrix_1D_mapped, local_block_size * rows, MPI_DOUBLE, matrix_local_block, local_block_size * rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter( rhs, local_block_size, MPI_DOUBLE, rhs_local_block, local_block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	setup_time = MPI_Wtime() - setup_start;
	kernel_start = MPI_Wtime();

 //--- PIVOT PART ---

        for(process = 0; process < size; process++){   
                if(process == rank){
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
		        pivots[0] = (double) rank;
		}
		
		//communicate pivots
		if (MPI_COMM_NULL != cur_comm){
		        //	send *pivots
		        mpi_start = MPI_Wtime();
		        MPI_Bcast(pivots, local_block_size * rows + local_block_size + 1, MPI_DOUBLE,0, cur_comm);
		        MPI_Group_excl(cur_group, 1, &zerorank, &cur_group);
		        MPI_Comm_create(cur_comm, cur_group, &cur_comm);
		        mpi_time += MPI_Wtime() - mpi_start; 
                }
                
                if (process < rank) {
                // make its chunk of A upper triangular and recompute rhs b
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
        }

        // --- SOLUTION PART ---
               	
	MPI_Comm_group(MPI_COMM_WORLD, &cur_group);
	MPI_Comm_create(MPI_COMM_WORLD, cur_group, &cur_comm);

        for(process = size-1; process>-1; process--){
	        if(rank == process){
	        //	compute local solutions(chunk of x)
	                for (row = (local_block_size - 1); row >= 0; row--) {
		                index = rank * local_block_size + row;
		                accumulation_buffer[row] = (double) index;
		                accumulation_buffer[local_block_size+row] = solution_local_block[row] = local_work_buffer[row];
		                for (i = (row - 1); i >= 0; i--){
			                local_work_buffer[i] -= solution_local_block[row] * matrix_local_block[ (i * rows) + index];
		                }
	                }
	        }
	        //communicate local solutions
	        if (MPI_COMM_NULL != cur_comm){
		                mpi_start = MPI_Wtime();
		                MPI_Bcast( accumulation_buffer, 2 * local_block_size, MPI_DOUBLE, process, cur_comm);
		                MPI_Group_excl(cur_group, 1, &process, &cur_group);
		                MPI_Comm_create(cur_comm, cur_group, &cur_comm);
		                mpi_time += MPI_Wtime() - mpi_start; 
                        }
                if(process > rank){
		        for (row  = (local_block_size - 1); row >= 0; row--) {
			        for (column = (local_block_size - 1);column >= 0; column--) {
				        index = (int) accumulation_buffer[column];
				        local_work_buffer[row] -= accumulation_buffer[local_block_size + column] * matrix_local_block[row * rows + index];
			        }
		        }
                }
        }

        // --- FINAL PART ---
		
	//	send/receive solutions
	mpi_start = MPI_Wtime();
	MPI_Gather(solution_local_block, local_block_size, MPI_DOUBLE, solution, local_block_size , MPI_DOUBLE, 0, MPI_COMM_WORLD);
	mpi_time += MPI_Wtime() - mpi_start;
	
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

