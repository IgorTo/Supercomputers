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
	MPI_Cart_coords(cartesian_grid_communicator, rank, 2, coordinates); //v COORDINATES imas shranjene koordinate procesa RANK

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
	if (rank == 0){
		int row, column;
		if ((fp = fopen (argv[1], "r")) != NULL){
			fscanf(fp, "%d %d\n", &matrices_a_b_dimensions[0], &matrices_a_b_dimensions[1]);
			A = (double **) malloc (matrices_a_b_dimensions[0] * sizeof(double *));
			for (row = 0; row < matrices_a_b_dimensions[0]; row++){
				A[row] = (double *) malloc(matrices_a_b_dimensions[1] * sizeof(double));
				for (column = 0; column < matrices_a_b_dimensions[1]; column++)
					fscanf(fp, "%lf", &A[row][column]);
			}
			fclose(fp);
		} else {
			if(rank == 0) fprintf(stderr, "error opening file for matrix A (%s)\n", argv[1]);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
		if((fp = fopen (argv[2], "r")) != NULL){
			fscanf(fp, "%d %d\n", &matrices_a_b_dimensions[2], &matrices_a_b_dimensions[3]);
			B = (double **) malloc (matrices_a_b_dimensions[2] * sizeof(double *));
			for(row = 0; row < matrices_a_b_dimensions[2]; row++){
				B[row] = (double *) malloc(matrices_a_b_dimensions[3] * sizeof(double *));
				for(column = 0; column < matrices_a_b_dimensions[3]; column++)
					fscanf(fp, "%lf", &B[row][column]);
			}
			fclose(fp);
		} else {
			if(rank == 0) fprintf(stderr, "error opening file for matrix B (%s)\n", argv[2]);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		// need to check that the multiplication is possible given dimensions
		// matrices_a_b_dimensions[0] = row size of A
		// matrices_a_b_dimensions[1] = column size of A
		// matrices_a_b_dimensions[2] = row size of B
		// matrices_a_b_dimensions[3] = column size of B
		if(matrices_a_b_dimensions[1] != matrices_a_b_dimensions[2]){
			if(rank == 0) fprintf(stderr, "A's column size (%d) must match B's row size (%d)\n",
					matrices_a_b_dimensions[1], matrices_a_b_dimensions[2]);
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		// this implementation is limited to cases where thematrices can be partitioned perfectly
		if( matrices_a_b_dimensions[0] % sqrt_size != 0
				|| matrices_a_b_dimensions[1] % sqrt_size != 0
				|| matrices_a_b_dimensions[2] % sqrt_size != 0
				|| matrices_a_b_dimensions[3] % sqrt_size != 0 ){
			if(rank == 0) fprintf(stderr, "cannot distribute work evenly among processe\n"
					"all dimensions (A: r:%d c:%d; B: r:%d c:%d) need to be divisible by %d\n",
					matrices_a_b_dimensions[0],matrices_a_b_dimensions[1],
					matrices_a_b_dimensions[2],matrices_a_b_dimensions[3], sqrt_size );
			MPI_Abort(MPI_COMM_WORLD, -1);
		}
	}

	// send dimensions to all peers
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	/*if(rank == 0) {
		int i;
		for(i = 1; i < size; i++){
			MPI_Send(matrices_a_b_dimensions, 4, MPI_INT, i, 0, cartesian_grid_communicator);
		}
	} else {
		MPI_Recv(matrices_a_b_dimensions, 4, MPI_INT, 0, 0, cartesian_grid_communicator, &status);
	}*/
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//has to be blocking, bcs it is used right afterwards...
	MPI_Bcast(matrices_a_b_dimensions, 4, MPI_INT, 0, cartesian_grid_communicator);
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

	// full arrays only needed at root
	if(rank == 0){
		A_array = (double *) malloc(sizeof(double) * A_rows * A_columns);
		B_array = (double *) malloc(sizeof(double) * B_rows * B_columns);
		C_array = (double *) malloc(sizeof(double) * A_rows * B_columns);
		// generate the 1D arrays of the matrices at root
		int row, column, i, j;
		for (i = 0; i < sqrt_size; i++){
			for (j = 0; j < sqrt_size; j++){
				for (row = 0; row < A_local_block_rows; row++){
					for (column = 0; column < A_local_block_columns; column++){
						A_array[((i * sqrt_size + j) * A_local_block_size) + (row * A_local_block_columns) + column]
							= A[i * A_local_block_rows + row][j * A_local_block_columns + column];
					}
				}
				for (row = 0; row < B_local_block_rows; row++){
					for (column = 0; column < B_local_block_columns; column++){
						B_array[((i * sqrt_size + j) * B_local_block_size) + (row * B_local_block_columns) + column]
							= B[i * B_local_block_rows + row][j * B_local_block_columns + column];
					}
				}
			}
		}
		// allocate output matrix C
		C = (double **) malloc(A_rows * sizeof(double *));
		for(i=0; i<A_rows ;i++){
			C[i] = (double *) malloc(B_columns * sizeof(double));
		}
	}

	// send a block to each process
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	/*if(rank == 0) {
		int i;
		for(i = 1; i < size; i++){
			MPI_Send((A_array + (i * A_local_block_size)), A_local_block_size, MPI_DOUBLE, i, 0, cartesian_grid_communicator);
			MPI_Send((B_array + (i * B_local_block_size)), B_local_block_size, MPI_DOUBLE, i, 0, cartesian_grid_communicator);

		}
		for(i = 0; i < A_local_block_size; i++){
			A_local_block[i] = A_array[i];
		}
		for(i = 0; i < B_local_block_size; i++){
			B_local_block[i] = B_array[i];
		}
	} else {
		MPI_Recv(A_local_block, A_local_block_size, MPI_DOUBLE, 0, 0, cartesian_grid_communicator, &status);
		MPI_Recv(B_local_block, B_local_block_size, MPI_DOUBLE, 0, 0, cartesian_grid_communicator, &status);
	}*/
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	/*MPI_Scatter(A_array,
    A_local_block_size, //int send_count,
    MPI_DOUBLE,
    A_local_block,
    A_local_block_size,
    MPI_DOUBLE,
    0,
    cartesian_grid_communicator);*/
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


	// fix initial arrangements before the core algorithm starts - fora je, da se preden se prvic zacne computational part of algo, moras ze bloke zamenjat...
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	/*if(coordinates[0] != 0){
		MPI_Sendrecv_replace(A_local_block, A_local_block_size, MPI_DOUBLE,
				(coordinates[1] + sqrt_size - coordinates[0]) % sqrt_size, 0,
				(coordinates[1] + coordinates[0]) % sqrt_size, 0, row_communicator, &status);
	}
	if(coordinates[1] != 0){
		MPI_Sendrecv_replace(B_local_block, B_local_block_size, MPI_DOUBLE,
				(coordinates[0] + sqrt_size - coordinates[1]) % sqrt_size, 0,
				(coordinates[0] + coordinates[1]) % sqrt_size, 0, column_communicator, &status);
	}*/
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// two independent scattervs one after another, so can be non-blocking, but barrier needed right afterwards, since data is than used...
	int *displsA[size];
	int *displsB[size];
	int *localblsizA[size];
	int *localblsizB[size];
	MPI_Request requests[2];
	MPI_Status statuses[2];
	for (int i=0; i<sqrt_size; i++){
		for (int j=0; j<sqrt_size; j++){
			displsA[i*sqrt_size + j] = (i*sqrt_size + (j+i)%sqrt_size)*A_local_block_size;
			displsB[i*sqrt_size + j] = (j + ((j+i)%size)*sqrt_size)*B_local_block_size;
			localblsizA[i*sqrt_size+j] = A_local_block_size;
			localblsizB[i*sqrt_size+j] = B_local_block_size;
		}
	}
	MPI_IScatterv(A_array, localblsizA, displsA,
                 MPI_DOUBLE, A_local_block, A_local_block_size,
                 MPI_DOUBLE, 0, cartesian_grid_communicator, &requests[0]);
  MPI_IScatterv(B_array, localblsizB, displsB,
                 MPI_DOUBLE, B_local_block, B_local_block_size,
                 MPI_DOUBLE, 0, cartesian_grid_communicator, &requests[1]);
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	//tega initial realignmenta ne bo treba, če boš že v A_array na zacetku alignano napisal data... ker te bo dovolj zgolj scatter.
	//Je pa Isaias tudi reku, da lahko to pustima in samo povema, da to pac ne gre prepisat, ker je if stavek in ne sodelujejo vsi ranki;
	//    	pogoj za collective je pa ravno to, da sodelujejo vsi!
	//
	//Je pa še ena moznost: scatterv! pa das primerne displacemente!
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	// cannon's algorithm
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555ta del je za SCATTER + GATHER
 	  int *dispA[sqrt_size];
	  int *dispB[sqrt_size];
		int *localsizesA[sqrt_size];
		int *localsizesB[sqrt_size];

		if (coordinates[0]==0) {
			double *B_rowarray[sqrt_size*B_local_block];
		}
		if (coordinates[1]==0){
			double *A_rowarray[sqrt_size*A_local_block];
		}

		for (int i=0; i<sqrt_size; i++){
			dispA[i] = ((i+1)%sqrt_size)*A_local_block_size;
			dispB[i] = ((i+1)%sqrt_size)*B_local_block_size;
			localsizesA[i] = A_local_block_size;
			localsizesB[i] = B_local_block_size;
    }
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
	int cannon_block_cycle;
	double compute_time = 0, mpi_time = 0, start;
	int C_index, A_row, A_column, B_column;

	MPI_Waitall(2, requests, statuses);
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
		//start = MPI_Wtime();
		// rotate blocks horizontally
		/*MPI_Sendrecv_replace(A_local_block, A_local_block_size, MPI_DOUBLE, //to bi slo z MPI_alltoallv, in tisto variablo za replacing. ampak bi blo inefficient - glej komentarje!
				(coordinates[1] + sqrt_size - 1) % sqrt_size, 0,
				(coordinates[1] + 1) % sqrt_size, 0, row_communicator, &status);
		// rotate blocks vertically
		MPI_Sendrecv_replace(B_local_block, B_local_block_size, MPI_DOUBLE,
				(coordinates[0] + sqrt_size - 1) % sqrt_size, 0,
				(coordinates[0] + 1) % sqrt_size, 0, column_communicator, &status);
		mpi_time += MPI_Wtime() - start; */



		//ce uporabis sendrecv, imas skupno v vsaki vrsti/stolpu sqrt_size komunikacij  (oz SQRT_SIZE posiljanj + SQRT_SIZE prejemanj),,
		//ce bi dala alltoall pa bi jih (kljub temu, da bi ponekod posiljal bloke velikosti 0) imel SIZE. Kar je pa tut ful inefficient(glej komentarje):
		/*This is allowed by the standard, but be warned that it is likely to perform
		poorly compared to what could be done with point-to-point or one-sided
		operations if most links are empty. ! ! ! ! ! ! ! ! */
		//lahko pa nardis to z gather+scatter, pa po en rank v vsaki vrsti/stolpu gathera vse, in jih zashiftano nazaj poslje. But that would still mean
		// 2*SQRT_SIZE communications, and it would have to be blocking, since data is used right afterwards. It might sound good to do this because even though
		//you need same amount of communication, the collectives are optimized; so the sole communcation should in this case take less time. However it's probably
		//not that big of a difference... lahko pa vseeno probata?
		//An even better idea seems to be if you figure out the pattern in which the blocks are shifted, and only use the A_array to scatter it from rnk 0
		//in the right order to all other ranks... This way we would all in all need SIZE communications (rnk 0 with everyone else) while with the previous
		//way we would all together need num_rows/colums*2*SQRT_SIZE, which is twice more. However in this last way, we would also need to compute the right
		//indeces for scatter everytime?
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCATTER + GATHER
		//rabis se rezervacijo prostora, za tiste ranke k si shranjo vse bloke neke vrste/stolpa
		start = MPI_Wtime();

		MPI_Gather(A_local_block, A_local_block_size, MPI_DOUBLE,
    					 A_rowarray, A_local_block_size, MPI_DOUBLE, coordinates[0]*sqrt_size, row_communicator);
		MPI_Gather(B_local_block, B_local_block_size, MPI_DOUBLE,
    					 B_rowarray, B_local_block_size, MPI_DOUBLE, coordinates[1], column_communicator);

		MPI_Scatterv(A_rowarray, localsizesA, dispA,
	                 MPI_DOUBLE, A_local_block, A_local_block_size,
	                 MPI_DOUBLE, coordinates[0]*sqrt_size, row_communicator);
		MPI_Scatterv(B_rowarray, localsizesB, dispB,
	                 MPI_DOUBLE, B_local_block, B_local_block_size,
	                 MPI_DOUBLE, coordinates[1], column_communicator);

 	  mpi_time += MPI_Wtime() - start;
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCATTER from origin
		start = MPI_Wtime(); //al bi se to moglo merit sele spodej, za zankami?

		for (int i=0; i<sqrt_size-1; i++){
			for (int j=0; j<sqrt_size-1; j++){
				displsA[i*sqrt_size + j] += A_local_block_size;
				displsB[i*sqrt_size + j] += B_local_block_size*size_sqrt;
			}
		}
		for (int i=0; i<sqrt_size; i++){
			displsA[size - sqrt_size + i] -= A_local_block_size*(sqrt_size-1);
			displsB[size - sqrt_size + i] -= B_local_block_size*(sqrt_size-1)*size_sqrt;
		}
		MPI_Scatterv(A_array, localblsizA, displsA,
	                 MPI_DOUBLE, A_local_block, A_local_block_size,
	                 MPI_DOUBLE, 0, cartesian_grid_communicator);
  	MPI_Scatterv(B_array, localblsizB, displsB,
 	                 MPI_DOUBLE, B_local_block, B_local_block_size,
 	                 MPI_DOUBLE, 0, cartesian_grid_communicator);

		mpi_time += MPI_Wtime() - start;
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	}


	// get C parts from other processes at rank 0
	/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
	}*/
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	MPI_Gather(C_local_block, A_local_block_rows * B_local_block_columns, MPI_DOUBLE,
               		 C_array, A_local_block_rows * B_local_block_columns, MPI_DOUBLE,
                   0, cartesian_grid_communicator);  //blocking, ker gres takoj nekaj delat s tem pol... right?
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	// generating output at rank 0
	if (rank == 0) {
		// convert the ID array into the actual C matrix
		int i, j, k, row, column;
		for (i = 0; i < sqrt_size; i++){  // block row index
			for (j = 0; j < sqrt_size; j++){ // block column index
				for (row = 0; row < A_local_block_rows; row++){
					for (column = 0; column < B_local_block_columns; column++){
						C[i * A_local_block_rows + row] [j * B_local_block_columns + column] =
							C_array[((i * sqrt_size + j) * A_local_block_rows * B_local_block_columns)
							+ (row * B_local_block_columns) + column];
					}
				}
			}
		}

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
	if(rank == 0){
		int i;
		for(i = 0; i < A_rows; i++){
			free(A[i]);
		}
		for(i = 0; i < B_rows; i++){
			free(B[i]);
		}
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