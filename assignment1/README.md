In this problem, you are asked to write a parallel algorithm for solving a set of dense linear equations of the form A*x = b, where A is an N x N matrix and x and b are column vectors. You will use Gaussian elimination without pivoting. The algorithm has two parts:
(a) Gaussian Elimination: The original system of equations is reduced to an upper triangular form
U x =y
where U is a matrix of size N x N in which all elements below the diagonal are zero, which are modified values of the A matrix. In addition, the diagonal elements have the value 1. The column vector y is the modified version of the b vector when you do the updates on the matrix and the vector in the Gaussian elimination stage.
(b) Back substitution: The new system of equations is solved to obtain the values of x.
The Gaussian elimination stage of the algorithm comprises N-1 steps. In the algorithm, the ith step eliminates nonzero subdiagonal elements in column i by subtracting the ith row from row j
in the range [i+1,n], in each case scaling the ith row by the factor Aji / Aii sa as to make the element Aji zero. See the figure below for illustration:
       0
Normalization
 Update Values
   Hence, the algorithm sweeps down the matrix from the top corner to the bottom right corner, leaving subdiagonal elements behind it.
The serial code for the algorithm is provided in the same hw1 directory for the previous questions. The file name is gauss.c
In the parallel algorithm, you must use dynamic scheduling. The whole point of parallel programming is performance, so you will be graded partially on the efficiency of your algorithm.
Suggestions:
Forking processes is very expensive. Your program does not need to fork more than once.
Each process should grab tasks until all tasks have been completed. Tasks of O(n) operations would be considered sufficiently large to hide the cost of obtaining a task, so a logical breakdown of tasks would be to let each task be responsible for the creation of one zero element in the lower diagonal of A.
You may observe a slowdown from 1 processor to 2 because the locking overhead is minimal with 1 processor.
Consider carefully the data dependencies in Gaussian elimination and the order in which tasks may be distributed.
Gaussian elimination involves O(n3) operations. The back substitution requires O(n2) operations, so you are not expected to parallelize back substitution.
     
The algorithm should scale, assuming n is much larger than the number of processors.
The machine will likely be overloaded the night(s) before the assignment is due. Do the
assignment early.
Send specific questions to the TA and the instructor.
a) [20 pts.] Write a shared memory parallel algorithm using explicit parallel programming constructs (pthread). Begin with the provided serial code “gauss.c”. Compile with
% g++ -o gauss gauss.c -lpthread
To test your program, use a command like (1024x1024 matrix, 4 processors, seed=5555):
./gauss 1024 4 5555
b) [20 pts.] Beginning with the same serial code, write a shared memory parallel algorithm using the implicit parallel programming directives ($DOACROSS or #pragma parallel, etc.). Compile with
% g++ -o gauss gauss.c -fopenmp
c) [10 pts.] Evaluate the performance of both algorithms on a given matrix size (5000x5000). Hand in the timing results of your code for 1, 4, and 8 processors.
In addition to a report describing the timing results, please submit your codes for Gaussian elimination (modified gauss.c for section a and b) using Canvas. There should be clear comments at the beginning of the code explaining your algorithm.
