#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

int main(int argc, char** argv) {
	//	A * x = b
	//	Matrix Names
	char matrix_name[200], vector_name[200], solution_name[200];
	char solution_name2[200] ="sample_sol.sol";
	char bin_matrix_name[200], bin_vector_name[200], bin_solution_name[200];
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

	if(rank == 0){
		printf("Solving the Ax=b system with Gaussian Elimination:\n");
		printf("READ:  Matrix   (A) file name: \"%s\"\n", matrix_name);
		printf("READ:  RHS      (b) file name: \"%s\"\n", vector_name);
		printf("WRITE: Solution (x) file name: \"%s\"\n", solution_name);
	}
	total_start = MPI_Wtime();

	int row, column, index;

	//Names for binary files
	sprintf(bin_matrix_name,   "%s_bin.mat", argv[1]);
	sprintf(bin_vector_name,   "%s_bin.vec", argv[1]);
	sprintf(bin_solution_name,	"%s_bin.sol", argv[1]);

	io_start = MPI_Wtime();
	//MPI Input: matrix A dimensions
	MPI_File fd;
	double *dim = (double*) malloc(2*sizeof(double));// contains dimensions of A : columns and rows
	//double *dim;
	MPI_File_open(MPI_COMM_WORLD, bin_matrix_name,MPI_MODE_RDONLY, MPI_INFO_NULL, &fd);
	MPI_File_read_all(fd, dim, 2,MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_close(&fd);
	columns = (int) *(dim);
	rows = (int) *(dim+1);

	//MPI Input: vector dim
	double double_rhs_rows;
	MPI_File fs;
	MPI_File_open(MPI_COMM_WORLD, bin_vector_name,MPI_MODE_RDONLY, MPI_INFO_NULL, &fs);
	MPI_File_read_all(fs, &double_rhs_rows, 1,MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_close(&fs);
	int rhs_rows = (int) double_rhs_rows;

	//	Some memory allocations and basic assertions
	if(rank == 0) {
		if(rows != columns) {
			perror("Only square matrices are allowed\n");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
		if(rows % size != 0) {
			perror("The matrix should be divisible by the number of processes\n");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
		if(rhs_rows != rows){
			perror("RHS rows must match the sizes of A");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
		index = 0;
	}
	/** Setup time will be added to IO time to avoid confusions**/
	//setup_start = MPI_Wtime();

	int i;
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


	//Initializing MPI-Input: matrix Agl
	MPI_File fh;
	MPI_File_open(MPI_COMM_WORLD, bin_matrix_name,MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_File_read_at_all(fh,  sizeof(double)*(rank*local_block_size*rows+2) , matrix_local_block, local_block_size * rows,MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_close(&fh);

	//Initializing MPI-Input: rhs b
	MPI_File fg;
	MPI_File_open(MPI_COMM_WORLD, bin_vector_name,MPI_MODE_RDONLY, MPI_INFO_NULL, &fg);
	MPI_File_read_at_all(fg,  sizeof(double)*(rank*local_block_size+1) , rhs_local_block, local_block_size,MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_close(&fg);


	io_time += MPI_Wtime() - io_start;
	//setup_time = MPI_Wtime() - setup_start;
	/** Setup time will be added to IO time to avoid confusions**/
	kernel_start = MPI_Wtime();

	//	receive *pivots from previous ranks, make its chunk of A upper triangular and recompute rhs b
	for(process = 0; process < rank; process++) {
		mpi_start = MPI_Wtime();
		MPI_Recv(pivots, (local_block_size * rows + local_block_size + 1), MPI_DOUBLE, process, process, MPI_COMM_WORLD, &status);
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
		MPI_Send( pivots, (local_block_size * rows + local_block_size + 1), MPI_DOUBLE, process, rank, MPI_COMM_WORLD);
		mpi_time += MPI_Wtime() - mpi_start;
	}

	//	receive chunks of rhs b after GE
	for (process = (rank + 1); process<size; ++process) {
		mpi_start = MPI_Wtime();
		MPI_Recv( accumulation_buffer, (2 * local_block_size), MPI_DOUBLE, process, process, MPI_COMM_WORLD, &status);
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
		MPI_Send( accumulation_buffer, (2 * local_block_size), MPI_DOUBLE, process, rank, MPI_COMM_WORLD);
		mpi_time += MPI_Wtime() - mpi_start;
	}

	kernel_time = MPI_Wtime() - kernel_start;

	//	basic assertions
	//	This part is not needed when solution and rhs are not stored
/**
	if (rank == 0) {
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
**/
	//MPI OUTPUT
	io_start = MPI_Wtime();
	MPI_File outfile;
	MPI_File_open(MPI_COMM_WORLD, bin_solution_name, MPI_MODE_CREATE | MPI_MODE_WRONLY,
			MPI_INFO_NULL, &outfile);
	/**
	MPI_File_set_view(outfile, sizeof(double)*rank*local_block_size, MPI_DOUBLE, MPI_DOUBLE,
			"native", MPI_INFO_NULL);
	MPI_File_write(outfile, solution_local_block, local_block_size,MPI_DOUBLE, MPI_STATUS_IGNORE);
	**/
	MPI_File_write_at_all(outfile, sizeof(double)*rank*local_block_size, solution_local_block, local_block_size,MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_close(&outfile);
	io_time += MPI_Wtime() - io_start;

	total_time = MPI_Wtime() - total_start;

	printf("[R%02d] Times: IO: %f; Setup: %f; Compute: %f; MPI: %f; Total: %f;\n",
			rank, io_time, setup_time, kernel_time, mpi_time, total_time);
/**
	if(rank == 0){
		for(i = 0; i < rows; i++){
			free(matrix_2d_mapped[i]);
		}
		free(matrix_2d_mapped);
		free(rhs);
		free(solution);
	}
**/
	free(matrix_local_block);
	free(rhs_local_block);
	free(pivots);
	free(local_work_buffer);
	free(accumulation_buffer);
	free(solution_local_block);

	MPI_Finalize();
	return 0;
}

