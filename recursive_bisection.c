#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NUM_POINTS 524288 

unsigned int X_axis[NUM_POINTS];
unsigned int Y_axis[NUM_POINTS];
double glo_cost = 0;
double starttime = 0;
double endtime;

void find_quadrants (int num_quadrants, int myid, int numprocs);
void cost_cal(int myid, int numprocs, int num_quadrants, int num_points);
void split_x(int split_points[], int quadrants, int points, unsigned int X_axis[], unsigned int Y_axis[]);
void split_y(int split_points[], int quadrants, int points, unsigned int X_axis[], unsigned int Y_axis[]);
unsigned int find_median(unsigned int *x, unsigned int *y, int num, int mid);
void selectSort(int m, unsigned int *medians_array, unsigned int *x, unsigned int *y, int num);
void cost_cal(int myid, int numprocs, int num_quadrants, int num_points);
void swap(unsigned int array[], int i, int j);
void print_result(int split_points[], int num_quadrants);




void find_quadrants (int num_quadrants, int myid, int numprocs) 
{
  /* YOU NEED TO FILL IN HERE */
  int split_points[num_quadrants];
  if (myid == 0) {
    int x_cut_flag = 0;  //flag for switch of x and y coordinary
    int quadrants = 1;

    

    // record the bisection coordinate

    while (num_quadrants > quadrants) {
      int num_points = NUM_POINTS / quadrants;

      if (!x_cut_flag) {
        split_x(split_points, quadrants, num_points, X_axis, Y_axis);
        x_cut_flag = 1;
      } else {
        split_y(split_points, quadrants, num_points, X_axis, Y_axis);
        x_cut_flag = 0;
        }
        quadrants *= 2;
    }
  }
  
  int chunkpoints = NUM_POINTS / num_quadrants;
  cost_cal(myid, numprocs, num_quadrants, chunkpoints);
  if (myid == 0) {
    endtime = MPI_Wtime();
    print_result(split_points, num_quadrants);
    printf("time cost: %f\n", endtime - starttime);
    printf("total cost: %lf\n", glo_cost);
  }
}


void cost_cal(int myid, int numprocs, int num_quadrants, int chunk_points) {
  
  MPI_Bcast(&X_axis, NUM_POINTS, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Y_axis, NUM_POINTS, MPI_INT, 0, MPI_COMM_WORLD);
  int i, j, k;
  double cost = 0;
  for (i = myid; i < num_quadrants; i += numprocs) {
    for (j = 0; j < chunk_points - 1; j++) {
      for (k = j + 1; k < chunk_points; k++) {
        double disx = (double)(abs(X_axis[chunk_points * i + j] - X_axis[chunk_points * i + k]));
        double disy = (double)(abs(Y_axis[chunk_points * i + j] - Y_axis[chunk_points * i + k]));
        double distance = sqrt(disx*disx + disy*disy);
        cost += distance;
        }
      }
    }
  MPI_Reduce(&cost, &glo_cost, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}



void split_x(int split_points[], int quadrants, int points,  unsigned int X_axis[],  unsigned int Y_axis[]) {
  int i, j, k;
  for (i = 0; i < quadrants; i++) {
        unsigned int x_pivot = find_median(X_axis + i * points, Y_axis + i * points, points, points / 2 - 1);;
        k = i * points;
        j = i * points + points - 1;
        while (k < j && k < k + points / 2) {
          if (X_axis[j] > x_pivot) {
            while (x_pivot < X_axis[j]) {
              j--;
            }
            if (k < j) {
              swap(X_axis, k, j);
              swap(Y_axis, k, j);
            }
          }
          k++;
        }
        split_points[i + quadrants - 1] = x_pivot;
  }
}

void split_y(int split_points[], int quadrants, int points,  unsigned int X_axis[],  unsigned int Y_axis[]) {
  split_x(split_points, quadrants, points, Y_axis, X_axis);
}

unsigned int find_median(unsigned int *x, unsigned int *y, int num, int mid) {
  if (num == 1) {
    return x[0];
  }

  if (mid == 0) {
    return x[0];
  }

  int chunk = 5;
  //divide the array in to 
  int m = (num + 4) / chunk;
  unsigned int *medians = (unsigned int *)malloc(m * sizeof(int));
  selectSort(m, medians, x, y, num);
  int pivot = find_median(medians, medians, m, m/2);
  free(medians);
  int i;
  int memory = 0;
  for (i = 0; i < num; i++) {
    if (x[i] == pivot) {
      swap(x, i, num - 1);
      swap(y, i, num - 1);
      break;
    }
  }
  for (i = 0; i < num - 1; i++) {
    if (x[i] < pivot) {
      swap(x, i, memory);
      swap(y, i, memory);
      memory++;
    }
  }
  swap(x, memory, num - 1);
  swap(y, memory, num - 1);

  if (memory == mid) {
    return pivot;
  } else if (memory > mid) {
    return find_median(x, y, memory, mid);
  } else {
    return find_median(x + memory + 1, y + memory + 1, num - memory - 1, mid - memory - 1);
  }
}


void selectSort(int m, unsigned int *medians_array, unsigned int *x, unsigned  int *y, int num) {
  int i, j, k;
  for (i = 0; i < m ; i++) {
    if (5*i + 4 < num) {
      unsigned int *w1 = x + 5 * i;
      unsigned int *w2 = y + 5 * i;
      for (j = 0; j < 3; j ++) {
        int min = j;
        for (k = j + 1; k < 5; k++) {
          if (w1[k] < w1[min]) {
            min = k;
          }
        }
        swap(w1, j, k);
        swap(w2, j, k);
      }
      medians_array[i] = w1[2];
    } else {
      medians_array[i] = x[5*i];
    }
  }
}

void swap(unsigned int array[], int i, int j) {
  int temp = array[i];
  array[i] = array[j];
  array[j] = temp;
}

void print_result(int split_points[], int num_quadrants) {
    int i, j, k;
    int top[num_quadrants];
    int left[num_quadrants];
    int right[num_quadrants];
    int bot[num_quadrants];


    int min_x, min_y, max_x, max_y;
    
    
    for (i = 1; i < NUM_POINTS; i++) {
      min_x = X_axis[i] < X_axis[0] ? X_axis[i] : X_axis[0];
      max_x = X_axis[i] > X_axis[0] ? X_axis[i] : X_axis[0];
      min_y = Y_axis[i] < Y_axis[0] ? Y_axis[i] : Y_axis[0];
      max_y = Y_axis[i] > Y_axis[0] ? Y_axis[i] : Y_axis[0];
    }

    int X_turn_flag = 0;
    int quadrants = 1;
    i = 0;

    printf("\nthe location of quadrants(bot-left, top-right):\n");
    while (i < num_quadrants - 1) {
      k = i;
      for (j = 0; j < quadrants; j ++) {
        top[j + quadrants] = top[j];
        bot[j + quadrants] = bot[j];
        left[j + quadrants] = left[j];
        right[j + quadrants] = right[j];
      }

      if (!X_turn_flag) {
        for (j = 0; j < 2 * quadrants; j += 2) {
          right[j] = split_points[k];
          left[j + 1] = split_points[k];
          
          k++;
        }
        X_turn_flag = 1;
      } else {
        for (j = 0; j < 2 * quadrants; j += 2) {
          bot[j] = split_points[k];
          top[j + 1] = split_points[k];
          k++;
        }
        X_turn_flag = 0;
      }

      i += quadrants;
      quadrants *= 2;
    }
    for (i = 0; i < num_quadrants; i++) {
      printf("quadrants %d :", i);
      printf("(%d, %d) ",left[i], bot[i]);
      printf("(%d, %d) ",left[i], top[i]);
      printf("(%d, %d) \n",right[i], top[i]);
      printf("(%d, %d) \n",right[i], bot[i]);
    }
  }



int main(argc,argv)
  int argc;
 char *argv[];
{
  int num_quadrants;
  int myid, numprocs;
  int  namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  
    
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Get_processor_name(processor_name,&namelen);

  if (argc != 2)
    {
      fprintf (stderr, "Usage: recursive_bisection <#of quadrants>\n");
      MPI_Finalize();
      exit (0);
    }

  fprintf (stderr,"Process %d on %s\n", myid, processor_name);

  num_quadrants = atoi (argv[1]);

  if (myid == 0)
    fprintf (stdout, "Extracting %d quadrants with %d processors \n", num_quadrants, numprocs);

  if (myid == 0)
    {
      int i;

      srand (10000);
      
      for (i = 0; i < NUM_POINTS; i++) {X_axis[i] = (unsigned int)rand();}
	         

      for (i = 0; i < NUM_POINTS; i++) {Y_axis[i] = (unsigned int)rand();}
	         

      starttime = MPI_Wtime();
    }

  MPI_Bcast(&X_axis, NUM_POINTS, MPI_INT, 0, MPI_COMM_WORLD);  
  MPI_Bcast(&Y_axis, NUM_POINTS, MPI_INT, 0, MPI_COMM_WORLD);  

  find_quadrants(num_quadrants, myid, numprocs);

 
  MPI_Finalize();
  return 0;
}
  

