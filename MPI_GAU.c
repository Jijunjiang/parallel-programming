#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <sys/time.h>
#include <math.h>
#include <limits.h>


#define MAXN 5000
#define N 5000

volatile float A[MAXN][MAXN], B[MAXN], X[MAXN];

unsigned int time_seed() {
  struct timeval t;
  struct timezone tzdummy;

  gettimeofday(&t, &tzdummy);
  return (unsigned int)(t.tv_usec);
}


void initialize_inputs() {
  int row, col;
  srand(5555);
  printf("\nInitializing...\n");
  for (col = 0; col < N; col++) {
    for (row = 0; row < N; row++) {
      A[row][col] = (float)rand() / 32768.0;
    }
    B[col] = (float)rand() / 32768.0;
    X[col] = 0.0;
  }

}



void print_inputs() {
  int row, col;

  if (N < 10) {
    printf("\nA =\n\t");
    for (row = 0; row < N; row++) {
      for (col = 0; col < N; col++) {
   printf("%5.2f%s", A[row][col], (col < N-1) ? ", " : ";\n\t");
      }
    }
    printf("\nB = [");
    for (col = 0; col < N; col++) {
      printf("%5.2f%s", B[col], (col < N-1) ? "; " : "]\n");
    }
  }
}

void print_X() {
  int row;

  if (N < 10) {
    printf("\nX = [");
    for (row = 0; row < N; row++) {
      printf("%5.2f%s", X[row], (row < N-1) ? "; " : "]\n");
    }
  }
}




int main( int argc, char *argv[])
{  
  
   double startwtime = 0.0, endwtime;
   int numprocs, myid, namelen;
   float multiplier;
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);
   MPI_Get_processor_name(processor_name,&namelen);
   fprintf(stdout,"Process %d on %s.\n",myid, processor_name);


   if (myid == 0) {
    initialize_inputs();
    print_inputs();
    printf("\nStarting clock.\n");
    startwtime = MPI_Wtime();
   }


   MPI_Bcast((void*)&A[0][0], MAXN*MAXN, MPI_FLOAT, 0, MPI_COMM_WORLD);
   MPI_Bcast((void*)&B[0], MAXN, MPI_FLOAT, 0, MPI_COMM_WORLD);
   
   
   int norms, row, col;
   for (norms = 0; norms < MAXN - 1; norms++) {
     // each loop, one column of the matrix woudl be the final result
      

      int i;
      //Since MPI_Bcast for other processer would wait for message, MPI_Bcast acts as a barriter
      MPI_Bcast((void*)&A[norms][0], MAXN, MPI_FLOAT, norms % numprocs, MPI_COMM_WORLD);
      MPI_Bcast((void*)&B[norms], 1, MPI_FLOAT, norms % numprocs, MPI_COMM_WORLD);

      for (i = myid; i < MAXN ; i += numprocs){
          if (i > norms){
            multiplier = A[i][norms] / A[norms][norms];
            for (col = norms; col < MAXN; col++) {
                       A[i][col] -= A[norms][col] * multiplier;
            }
            B[i] -= B[norms] * multiplier;
          } else {
            continue;
          }
      }   
  }

  MPI_Bcast((void*)&A[MAXN - 1][0], MAXN, MPI_FLOAT, (MAXN - 1) % numprocs, MPI_COMM_WORLD);
  MPI_Bcast((void*)&B[MAXN - 1], 1, MPI_FLOAT, (MAXN - 1) % numprocs, MPI_COMM_WORLD);
  
  
   if ( myid == 0) {
     endwtime = MPI_Wtime();
     for (row = MAXN - 1; row >= 0; row--) {
        X[row] = B[row];
        for (col = MAXN-1; col > row; col--) {
          X[row] -= A[row][col] * X[col];
        }
        X[row] /= A[row][row];
      }
      print_X();
      printf("wall clock time = %f\n", endwtime-startwtime); 
    }
  MPI_Finalize();
  return 0;

}


