#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

int main (int argc, char **argv) {
	FILE *fp;
	double **A = NULL, **B = NULL, **C = NULL, *A_array = NULL, *B_array = NULL, *C_array = NULL;
	double *A_local_block = NULL, *B_local_block = NULL, *C_local_block = NULL;
	int A_rows, A_columns, A_local_block_rows, A_local_block_columns, A_local_block_size;
	int B_rows, B_columns, B_local_block_rows, B_local_block_columns, B_local_block_size;
	int rank, size, sqrt_size, matrices_a_b_dimensions[4];
	MPI_Comm cartesian_grid_communicator, row_communicator, column_communicator;
	MPI_Status status;
	MPI_File fhA, fhB, fhC;
	MPI_Status read_statusA, read_statusB, write_status;

	// used to manage the cartesian grid
	int dimensions[2], periods[2], coordinates[2], remain_dims[2];

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* For square mesh */
	sqrt_size = (int)sqrt((double) size);
	if(sqrt_size * sqrt_size != size){
		if( rank == 0 ) perror("need to run mpiexec with a perfect square number of processes\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	// create a 2D cartesian grid
	dimensions[0] = dimensions[1] = sqrt_size;
	periods[0] = periods[1] = 1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, periods, 1, &cartesian_grid_communicator);
	MPI_Cart_coords(cartesian_grid_communicator, rank, 2, coordinates);

	// create a row communicator
	remain_dims[0] = 0;
	remain_dims[1] = 1;
	MPI_Cart_sub(cartesian_grid_communicator, remain_dims, &row_communicator);

	// create a column communicator
	remain_dims[0] = 1;
	remain_dims[1] = 0;
	MPI_Cart_sub(cartesian_grid_communicator, remain_dims, &column_communicator);

	// getting matrices from files at rank 0 only
	// example: mpiexec -n 64 ./cannon matrix1 matrix2 [test]

	//bufsize = FILESIZE/nprocs;
	//nints = bufsize/sizeof(int);
	successfulA = MPI_File_open(cartesian_grid_communicator, argv[1],
              							  MPI_MODE_RDONLY, MPI_INFO_NULL, &fhA);
	successfulB = MPI_File_open(cartesian_grid_communicator, argv[2],
              							  MPI_MODE_RDONLY, MPI_INFO_NULL, &fhB);
	if (rank==0){
		if (successfulA!=MPI_SUCCESS){
			fprintf(stderr, "error opening file for matrix A (%s)\n", argv[1]);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
		else if (successfulB!=MPI_SUCCESS){
			fprintf(stderr, "error opening file for matrix B (%s)\n", argv[2]);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
		else {
			MPI_File_read(fhA, &matrices_a_b_dimensions[0], 2/sizeof(int), MPI_INT, &read_statusA);
			MPI_File_read(fhB, &matrices_a_b_dimensions[2], 2/sizeof(int), MPI_INT, &read_statusB);

			//checking, if matrices are suitable
			if(matrices_a_b_dimensions[1] != matrices_a_b_dimensions[2]){
				fprintf(stderr, "A's column size (%d) must match B's row size (%d)\n",
								matrices_a_b_dimensions[1], matrices_a_b_dimensions[2]);
				MPI_Abort(MPI_COMM_WORLD, -1);
			}
			// this implementation is limited to cases where thematrices can be partitioned perfectly
			if( matrices_a_b_dimensions[0] % sqrt_size != 0
					|| matrices_a_b_dimensions[1] % sqrt_size != 0
					|| matrices_a_b_dimensions[2] % sqrt_size != 0
					|| matrices_a_b_dimensions[3] % sqrt_size != 0 ){
				fprintf(stderr, "cannot distribute work evenly among processes\n"
						"all dimensions (A: r:%d c:%d; B: r:%d c:%d) need to be divisible by %d\n",
						matrices_a_b_dimensions[0],matrices_a_b_dimensions[1],
						matrices_a_b_dimensions[2],matrices_a_b_dimensions[3], sqrt_size );
				MPI_Abort(MPI_COMM_WORLD, -1);
			}

			//send info about length that needs to be read to others
			MPI_Bcast(matrices_a_b_dimensions, 4, MPI_INT, 0, cartesian_grid_communicator);
		}
	}
	else {
		//get info about size of matrices from root
		MPI_Bcast(matrices_a_b_dimensions, 4, MPI_INT, 0, cartesian_grid_communicator);
	}
  // set variables and allocate memory
	A_rows = matrices_a_b_dimensions[0];
	A_columns = matrices_a_b_dimensions[1];
	B_rows = matrices_a_b_dimensions[2];
	B_columns = matrices_a_b_dimensions[3];

	// local metadata for A
	A_local_block_rows = A_rows / sqrt_size;
	A_local_block_columns = A_columns / sqrt_size;
	A_local_block_size = A_local_block_rows * A_local_block_columns;
	A_local_block = (double *) malloc (A_local_block_size * sizeof(double));

	// local metadata for B
	B_local_block_rows = B_rows / sqrt_size;
	B_local_block_columns = B_columns / sqrt_size;
	B_local_block_size = B_local_block_rows * B_local_block_columns;
	B_local_block = (double *) malloc (B_local_block_size * sizeof(double));

	// local metadata for C
	C_local_block = (double *) malloc (A_local_block_rows * B_local_block_columns * sizeof(double));
	// C needs to be initialized at 0 (accumulates partial dot-products)
	int i;
	for(i=0; i < A_local_block_rows * B_local_block_columns; i++){
		C_local_block[i] = 0;
	}
	//now read the rest of file (the actual matrices)
	MPI_File_seek(fhA, rank * A_local_block_size + 2, MPI_SEEK_SET); //+2 for the dimension we already read
	MPI_File_seek(fhB, rank * B_local_block_size + 2, MPI_SEEK_SET); //+2 for the dimension we already read
	MPI_File_read(fhA, A_local_block, A_local_block_size/sizeof(int), MPI_INT, &read_statusA);
	MPI_File_read(fhB, B_local_block, B_local_block_size/sizeof(int), MPI_INT, &read_statusB);

	//now close the files
	MPI_File_close(&fhA);
	MPI_File_close(&fhB);


	// fix initial arrangements before the core algorithm starts <- tega se lahko znebis, ce pravilno nastavis branje fajla!
	if(coordinates[0] != 0){
		MPI_Sendrecv_replace(A_local_block, A_local_block_size, MPI_DOUBLE,
				(coordinates[1] + sqrt_size - coordinates[0]) % sqrt_size, 0,
				(coordinates[1] + coordinates[0]) % sqrt_size, 0, row_communicator, &status);
	}
	if(coordinates[1] != 0){
		MPI_Sendrecv_replace(B_local_block, B_local_block_size, MPI_DOUBLE,
				(coordinates[0] + sqrt_size - coordinates[1]) % sqrt_size, 0,
				(coordinates[0] + coordinates[1]) % sqrt_size, 0, column_communicator, &status);
	}

	// cannon's algorithm
	int cannon_block_cycle;
	double compute_time = 0, mpi_time = 0, start;
	int C_index, A_row, A_column, B_column;
	for(cannon_block_cycle = 0; cannon_block_cycle < sqrt_size; cannon_block_cycle++){
		// compute partial result for this block cycle
		start = MPI_Wtime();
		for(C_index = 0, A_row = 0; A_row < A_local_block_rows; A_row++){
			for(B_column = 0; B_column < B_local_block_columns; B_column++, C_index++){
				for(A_column = 0; A_column < A_local_block_columns; A_column++){
					C_local_block[C_index] += A_local_block[A_row * A_local_block_columns + A_column] *
						B_local_block[A_column * B_local_block_columns + B_column];
				}
			}
		}
		compute_time += MPI_Wtime() - start;
		start = MPI_Wtime();
		// rotate blocks horizontally
		MPI_Sendrecv_replace(A_local_block, A_local_block_size, MPI_DOUBLE,
				(coordinates[1] + sqrt_size - 1) % sqrt_size, 0,
				(coordinates[1] + 1) % sqrt_size, 0, row_communicator, &status);
		// rotate blocks vertically
		MPI_Sendrecv_replace(B_local_block, B_local_block_size, MPI_DOUBLE,
				(coordinates[0] + sqrt_size - 1) % sqrt_size, 0,
				(coordinates[0] + 1) % sqrt_size, 0, column_communicator, &status);
		mpi_time += MPI_Wtime() - start;
	}

	// get C parts from other processes at rank 0
	//TODO WRITING TO FILE!
	if(rank == 0) {
		for(i = 0; i < A_local_block_rows * B_local_block_columns; i++){
			C_array[i] = C_local_block[i];
		}
		int i;
		for(i = 1; i < size; i++){
			MPI_Recv(C_array + (i * A_local_block_rows * B_local_block_columns), A_local_block_rows * B_local_block_columns,
				MPI_DOUBLE, i, 0, cartesian_grid_communicator, &status);
		}
	} else {
		MPI_Send(C_local_block, A_local_block_rows * B_local_block_columns, MPI_DOUBLE, 0, 0, cartesian_grid_communicator);
	}

	// generating output at rank 0
	if (rank == 0) {
		printf("(%d,%d)x(%d,%d)=(%d,%d)\n", A_rows, A_columns, B_rows, B_columns, A_rows, B_columns);
		printf("Computation time: %lf\n", compute_time);
		printf("MPI time:         %lf\n", mpi_time);

		if (argc == 4){
			// present results on the screen
			printf("\nA( %d x %d ):\n", A_rows, A_columns);
			for(row = 0; row < A_rows; row++) {
				for(column = 0; column < A_columns; column++)
					printf ("%7.3f ", A[row][column]);
				printf ("\n");
			}
			printf("\nB( %d x %d ):\n", B_rows, B_columns);
			for(row = 0; row < B_rows; row++){
				for(column = 0; column < B_columns; column++)
					printf("%7.3f ", B[row][column]);
				printf("\n");
			}
			printf("\nC( %d x %d ) = AxB:\n", A_rows, B_columns);
			for(row = 0; row < A_rows; row++){
				for(column = 0; column < B_columns; column++)
					printf("%7.3f ",C[row][column]);
				printf("\n");
			}


			printf("\nPerforming serial consistency check. Be patient...\n");
			fflush(stdout);
			int pass = 1;
			double temp;
			for(i=0; i<A_rows; i++){
				for(j=0; j<B_columns; j++){
					temp = 0;
					for(k=0; k<B_rows; k++){
						temp += A[i][k] * B[k][j];
					}
					printf("%7.3f ", temp);
					if(temp != C[i][j]){
						pass = 0;
					}
				}
				printf("\n");
			}
			if (pass) printf("Consistency check: PASS\n");
			else printf("Consistency check: FAIL\n");
		}
	}

	// free all memory
	free(A_local_block);
	free(B_local_block);
	free(C_local_block);

	// finalize MPI
	MPI_Finalize();
}
