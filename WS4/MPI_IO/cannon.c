#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
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

	// used to manage the cartesian grid
	int dimensions[2], periods[2], coordinates[2], remain_dims[2];

	double time_init;
	clock_t begin, end;

	begin = clock();
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	end = clock();
	if (rank == 0) time_init = (double)(end - begin);

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
	MPI_Cart_coords(cartesian_grid_communicator, rank, 2, coordinates); //v COORDINATES imas shranjene koordinate procesa RANK

	// create a row communicator
	remain_dims[0] = 0;
	remain_dims[1] = 1;
	MPI_Cart_sub(cartesian_grid_communicator, remain_dims, &row_communicator);

	// create a column communicator
	remain_dims[0] = 1;
	remain_dims[1] = 0;
	MPI_Cart_sub(cartesian_grid_communicator, remain_dims, &column_communicator);


	// set time variables for different MPI parts
	double read_time, send_dim_time, send_blocks_time, gather_time, write_time;


	read_time = MPI_Wtime();
	// getting matrices from files at rank 0 only
	// example: mpiexec -n 64 ./cannon matrix1 matrix2 [test]

	//the following reading of the dimensions is still needed, to check if everything is okay + to know how much each process will need to read.
	if (rank == 0){

		int row, column;
		if ((fp = fopen (argv[1], "r")) != NULL){
			fscanf(fp, "%d %d\n", &matrices_a_b_dimensions[0], &matrices_a_b_dimensions[1]);
			fclose(fp);
		} else {
			fprintf(stderr, "error opening file for matrix A (%s)\n", argv[1]);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
		if((fp = fopen (argv[2], "r")) != NULL){
			fscanf(fp, "%d %d\n", &matrices_a_b_dimensions[2], &matrices_a_b_dimensions[3]);
			fclose(fp);
		} else {
			fprintf(stderr, "error opening file for matrix B (%s)\n", argv[2]);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		// need to check that the multiplication is possible given dimensions
		// matrices_a_b_dimensions[0] = row size of A
		// matrices_a_b_dimensions[1] = column size of A
		// matrices_a_b_dimensions[2] = row size of B
		// matrices_a_b_dimensions[3] = column size of B
		if(matrices_a_b_dimensions[1] != matrices_a_b_dimensions[2]){
			fprintf(stderr, "A's column size (%d) must match B's row size (%d)\n", matrices_a_b_dimensions[1], matrices_a_b_dimensions[2]);
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
	}
	//to continue with reading, we need info on matrix dimensions. So we do first broadcasting...
	read_time -= MPI_Wtime();


	send_dim_time = MPI_Wtime();
	// send dimensions to all peers
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//has to be blocking, bcs data is used right afterwards...
	MPI_Bcast(matrices_a_b_dimensions, 4, MPI_INT, 0, cartesian_grid_communicator);
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	send_dim_time -= MPI_Wtime();


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
	int i,j;
	for(i=0; i < A_local_block_rows * B_local_block_columns; i++){
		C_local_block[i] = 0;
	}


	
	//here the actual IO reading takes place;
	// MODE = readonly + unique(to get rid of locks and make it faster). Definition: The MPI_MODE_UNIQUE_OPEN mode allows an 
	//                 implementation to optimize access by eliminating the overhead of file locking. It is erroneous to open 
	//                 a file in this mode unless the file will not be concurrently opened elsewhere.
	// INFO = no particularly useful info for read (what about max read buff?), so set to null
 	// FILE POINTER = output, returns pointer to the read file, which is given as string (arg[.]).
	// DISP = 2 = displacement, or the length of header that is to be ignored in the matrix files
	// SUBARRAY_TYPE = type of the underlying data structure, a custom datatype that we create based on how our matrix file looks like
	// STARTS = the upper left corner of each rank's subarray 

	MPI_File fhA, fhB;
	MPI_Datatype subarray_typeA, subarray_typeB;
	int lsizesA[] = {A_local_block_rows, 4*A_local_block_columns};
	int lsizesB[] = {B_local_block_rows, 4*B_local_block_columns};
	int startsA[] = {coordinates[0]*A_local_block_rows, ((coordinates[0] + coordinates[1])%sqrt_size)*A_local_block_columns * 4};
	int startsB[] = {coordinates[0]*B_local_block_rows, ((coordinates[0] + coordinates[1])%sqrt_size)*B_local_block_columns * 4};

	// Correct the file style matrix dimension
	int file_matrices_a_b_dimensions[] = { matrices_a_b_dimensions[0], matrices_a_b_dimensions[1] * 4 }; 

	// Reservation for array where file is going to be read. We need 4 (number.number[space]) times bigger array than the number of read values.
	char *larrayA = (char *) malloc (A_local_block_size * 4 * sizeof(char));
	char *larrayB = (char *) malloc (B_local_block_size * 4 * sizeof(char));

	MPI_Type_create_subarray(2, file_matrices_a_b_dimensions, lsizesA, startsA, MPI_ORDER_C, MPI_CHAR, &subarray_typeA);
	MPI_Type_create_subarray(2, file_matrices_a_b_dimensions, lsizesB, startsB, MPI_ORDER_C, MPI_CHAR, &subarray_typeB);

	//now we need to commit created datatypes
	MPI_Type_commit(&subarray_typeA);
	MPI_Type_commit(&subarray_typeB);

        /* We need to know how long the first line of a file is (number x number <- pattern), so that we know where to start reading actual data. */
	int displacement;
	switch( (int) (log(matrices_a_b_dimensions[0]) / log(10)) + 1) { // How many ciphres does the number have?
		case 1: displacement = 4; break;
		case 2: displacement = 6; break;
		case 3: displacement = 8; break;
		case 4: displacement = 10; break;
		case 5: displacement = 12; break;
	}

	read_time += MPI_Wtime(); //we now continue with the reading part

	//now we read each file and return info about reading success when failed
	if(MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &fhA) != MPI_SUCCESS){
		fprintf(stderr, "error opening file for matrix A (%s) at rank %d\n", argv[1], rank );
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	//after oppening, we need to set proper file views for the ranks:
	MPI_File_set_view (fhA, displacement*sizeof(char), MPI_CHAR, subarray_typeA, "native", MPI_INFO_NULL);

	//now we can actually read our part
	MPI_File_read_all(fhA, larrayA, A_local_block_size * 4, MPI_CHAR, MPI_STATUS_IGNORE);

	//now we can close the file
	MPI_File_close(&fhA);

	// postprocess the ASCII data for A_local_block
	for (i=0; i < A_local_block_size; i++){
		char tmpString[4];
		strncpy( tmpString, &larrayA[i*4], 4 );

		double number;
		sscanf(tmpString, "%lf ", &number);

		A_local_block[i] = number;
	}

	// now do the same (reading) things for matrix B
	if(MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &fhB) != MPI_SUCCESS){
		fprintf(stderr, "error opening file for matrix B (%s) at rank %d\n", argv[2], rank );
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	MPI_File_set_view (fhB, displacement*sizeof(char), MPI_CHAR, subarray_typeB, "native", MPI_INFO_NULL);
	MPI_File_read_all(fhB, larrayB, B_local_block_size * 4, MPI_CHAR, MPI_STATUS_IGNORE);
	MPI_File_close(&fhB);

        // postprocess the ASCII data for B_local_block
        for (i=0; i < B_local_block_size; i++){
                char tmpString[4];
                strncpy( tmpString, &larrayB[i*4], 4 );

                double number;
                sscanf(tmpString, "%lf ", &number);

                B_local_block[i] = number;
        }

	read_time -= MPI_Wtime();

	// C arrays init
	if (rank == 0){
		// !!! This is needed to gather everything at the end !!!
		C_array = (double *) malloc(sizeof(double) * A_rows * B_columns);
                
		// allocate output matrix C
                C = (double **) malloc(A_rows * sizeof(double *));
                for(i=0; i<A_rows ;i++){
                        C[i] = (double *) malloc(B_columns * sizeof(double));
                }
	}


	// !!!!!! cannon's algorithm !!!!!!!
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
		MPI_Sendrecv_replace(A_local_block, A_local_block_size, MPI_DOUBLE, //to bi slo z MPI_alltoallv, in tisto variablo za replacing. ampak bi blo inefficient - glej komentarje!
				(coordinates[1] + sqrt_size - 1) % sqrt_size, 0,
				(coordinates[1] + 1) % sqrt_size, 0, row_communicator, &status);
		// rotate blocks vertically
		MPI_Sendrecv_replace(B_local_block, B_local_block_size, MPI_DOUBLE,
				(coordinates[0] + sqrt_size - 1) % sqrt_size, 0,
				(coordinates[0] + 1) % sqrt_size, 0, column_communicator, &status);
		mpi_time += MPI_Wtime() - start;
	}


	// get C parts from other processes at rank 0
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	gather_time =  MPI_Wtime();
	MPI_Gather(C_local_block, A_local_block_rows * B_local_block_columns, MPI_DOUBLE,
               		 C_array, A_local_block_rows * B_local_block_columns, MPI_DOUBLE,
                   0, cartesian_grid_communicator);  //blocking, ker gres takoj nekaj delat s tem pol... right?
	gather_time -= MPI_Wtime();
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	// -------- Writing to the File ----------
	
	MPI_File fhOutput;
	char outputFile[20];
	sprintf(outputFile, "resultFiles/%dx%d.output", matrices_a_b_dimensions[2], matrices_a_b_dimensions[3]);

	// Make a file type for C matrix
	MPI_Datatype C_style;
	int C_sizes[] = {matrices_a_b_dimensions[0], matrices_a_b_dimensions[3]};
	int C_whereToStart[] = {coordinates[0] * A_local_block_rows, coordinates[1] * B_local_block_columns};
	int C_local_sizes[] = {A_local_block_rows, B_local_block_columns};

	MPI_Type_create_subarray(2, C_sizes, C_local_sizes, C_whereToStart, MPI_ORDER_C, MPI_DOUBLE, &C_style);
	MPI_Type_commit(&C_style);

	
	// Open, set view, write and close. Done. ;)
	write_time = MPI_Wtime();

	MPI_File_open(MPI_COMM_WORLD, outputFile, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fhOutput);
	MPI_File_set_view(fhOutput, 0, MPI_DOUBLE, C_style, "native", MPI_INFO_NULL);

	MPI_File_write(fhOutput, C_local_block, A_local_block_rows * B_local_block_columns, MPI_DOUBLE, MPI_STATUS_IGNORE);

	MPI_File_close(&fhOutput);

	write_time -= MPI_Wtime();
	// ----------------------- END writing to the file -----------------

	// generating output at rank 0
	if (rank == 0) {
		printf("\n");
		printf("(%d,%d)x(%d,%d)=(%d,%d)\n", A_rows, A_columns, B_rows, B_columns, A_rows, B_columns);
		printf("Computation time: %lf\n", compute_time);
		printf("MPI init time: %lf\n", -time_init);
		printf("MPI time:         %lf\n", mpi_time);
		printf("Read time:        %lf\n", -read_time);
		printf("Send dims time:   %lf\n", -send_dim_time);
		printf("Send blocks time: %lf\n", -send_blocks_time);
		printf("Gather time:      %lf\n", -gather_time);
		printf("Write time:       %lf\n", -write_time);

	}

	// free all memory
	if(rank == 0){
		int i;
	
		for(i = 0; i < A_rows; i++){
			free(C[i]);
	}
		free(A);
		free(B);
		free(C);
		free(A_array);
		free(B_array);
		free(C_array);
	}
	free(A_local_block);
	free(B_local_block);
	free(C_local_block);

	// finalize MPI
	MPI_Finalize();
}

