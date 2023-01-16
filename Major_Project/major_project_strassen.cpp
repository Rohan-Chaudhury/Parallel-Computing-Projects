#include <omp.h>
#include <stdio.h>   
#include <stdlib.h>  
#include <time.h>
#include <fstream>
#include <iostream>
#include <ctime>
#include <math.h>
#include <cmath>



int **Allocate2DArray( int rowSize, int colSize) {

    // Allocating the memory for array of column elements 

    int **a = new int*[rowSize];

    // Allocating memory for array of elements of each row

    int *currentPtr = new int [rowSize * colSize];

    for( int i = 0; i < rowSize; ++i)
    {
        *(a + i) = currentPtr;
         currentPtr += colSize;
    }
    return a;
}


void Free2DArray(int** Array) {
    delete [] *Array;
    delete [] Array;
}


void initializeMatrices(int **A, int **B, int size){
	int i, j;

    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j) {
	        B[i][j] = 0;
			/*Randomly initializing the matrix values*/
	        A[i][j] = rand()%50;
	    }

    for (int i = 0; i < size; ++i)
        B[i][i] = 1;
}

int Threshold;
double Size;

void matrixMultiplication(int afirst, int alast, int bfirst, int blast, int cfirst, int clast, int **A, int **B, int **C) {
/*  
Function uses a triple loop to multiply a submatrix of A with 
a submatrix of B and stores it in submatrix of C   
*/  
// afirst is the first i index */
// alast is the last+1 i index */   
// bfirst is the first j index */ 
// blast is the last+1 j index */ 
// cfirst is the first k index */  
// clast is the last+1 k index */  

	for (int i = afirst; i < alast; i++) {

		for (int j = bfirst; j < blast; j++) {

			C[i][j] = 0.0;  

			for (int k = cfirst; k < clast; k++)   {

				C[i][j] += A[i][k]*B[k][j];  

			}
		} 
	}
}   
  

void matrixCopy(int **X, int size, int **Y, int xfirst, int yfirst) {
	/*
	Function to copy matrix
	*/
	for (int i = 0; i < size; i++) {
		X[i] = &Y[xfirst+i][yfirst];
	}
}

void SubPartialMat(int **T, int m, int n, int **X, int **Y) {

	for (int i = 0; i < m; i++) {

		for (int j = 0; j < n; j++) {

			T[i][j] = X[i][j] - Y[i][j];

		}
	}
}


void AddPartialMat(int **T, int m, int n, int **X, int **Y) {

	for (int i = 0; i < m; i++) {

		for (int j = 0; j < n; j++) {

			T[i][j] = X[i][j] + Y[i][j];

		}
	}
}



void strassenMatrixMulAlgo(int mat1Size, int mat2Size, int mat3Size, int **A, int **B, int **C) {
	if (((float)mat1Size) <= (Size / pow(2,Threshold)) ) {
		matrixMultiplication(0, mat1Size, 0, mat2Size, 0, mat3Size, A, B, C);
	} else {

		int mat1Mid = mat1Size/2;
		int mat2Mid = mat2Size/2;
		int mat3Mid = mat3Size/2;


		int **M1 = Allocate2DArray(mat1Mid, mat2Mid);
		int **M2 = Allocate2DArray(mat1Mid, mat2Mid);
		int **M3 = Allocate2DArray(mat1Mid, mat2Mid);
		int **M4 = Allocate2DArray(mat1Mid, mat2Mid);
		int **M5 = Allocate2DArray(mat1Mid, mat2Mid);
		int **M6 = Allocate2DArray(mat1Mid, mat2Mid);
		int **M7 = Allocate2DArray(mat1Mid, mat2Mid);


		int **AMat_1 = Allocate2DArray(mat1Mid, mat3Mid);
		int **BMat_1 = Allocate2DArray(mat3Mid, mat2Mid);
		int **AMat_2 = Allocate2DArray(mat1Mid, mat3Mid);
		int **BMat_3 = Allocate2DArray(mat3Mid, mat2Mid);
		int **BMat_4 = Allocate2DArray(mat3Mid, mat2Mid);
		int **AMat_5 = Allocate2DArray(mat1Mid, mat3Mid);
		int **AMat_6 = Allocate2DArray(mat1Mid, mat3Mid);
		int **BMat_6 = Allocate2DArray(mat3Mid, mat2Mid);
		int **AMat_7 = Allocate2DArray(mat1Mid, mat3Mid);
		int **BMat_7 = Allocate2DArray(mat3Mid, mat2Mid);


		int **A11 = new int*[mat1Mid];
		int **A12 = new int*[mat1Mid];
		int **A21 = new int*[mat1Mid];
		int **A22 = new int*[mat1Mid];


		int **B11 = new int*[mat3Mid];
		int **B12 = new int*[mat3Mid];
		int **B21 = new int*[mat3Mid];
		int **B22 = new int*[mat3Mid];


		int **C11 = new int*[mat1Mid];
		int **C12 = new int*[mat1Mid];
		int **C21 = new int*[mat1Mid];
		int **C22 = new int*[mat1Mid];


		matrixCopy(A11, mat1Mid, A,  0,  0);
		matrixCopy(A12, mat1Mid, A,  0, mat3Mid);
		matrixCopy(A21, mat1Mid, A, mat1Mid,  0);
		matrixCopy(A22, mat1Mid, A, mat1Mid, mat3Mid);


		matrixCopy(B11, mat3Mid, B,  0,  0);
		matrixCopy(B12, mat3Mid, B,  0, mat2Mid);
		matrixCopy(B21, mat3Mid, B, mat3Mid,  0);
		matrixCopy(B22, mat3Mid, B, mat3Mid, mat2Mid);


		matrixCopy(C11, mat1Mid, C,  0,  0);
		matrixCopy(C12, mat1Mid, C,  0, mat2Mid);
		matrixCopy(C21, mat1Mid, C, mat1Mid,  0);
		matrixCopy(C22, mat1Mid, C, mat1Mid, mat2Mid);


#pragma omp task
		{
	/* M1 = (A11 + A22)*(B11 + B22) */

		AddPartialMat(AMat_1, mat1Mid, mat3Mid, A11, A22);
		AddPartialMat(BMat_1, mat3Mid, mat2Mid, B11, B22);
		strassenMatrixMulAlgo(mat1Mid, mat2Mid, mat3Mid, AMat_1, BMat_1, M1);
		}

#pragma omp task
		{
	/* M2 = (A21 + A22)*B11 */

		AddPartialMat(AMat_2, mat1Mid, mat3Mid, A21, A22);
		strassenMatrixMulAlgo(mat1Mid, mat2Mid, mat3Mid, AMat_2, B11, M2);
		}

#pragma omp task
		{
	/* M3 = A11*(B12 - B22) */

		SubPartialMat(BMat_3, mat3Mid, mat2Mid, B12, B22);
		strassenMatrixMulAlgo(mat1Mid, mat2Mid, mat3Mid, A11, BMat_3, M3);
		}

#pragma omp task
		{
	/* M4 = A22*(B21 - B11) */

		SubPartialMat(BMat_4, mat3Mid, mat2Mid, B21, B11);
		strassenMatrixMulAlgo(mat1Mid, mat2Mid, mat3Mid, A22, BMat_4, M4);
		}

#pragma omp task
		{
	/* M5 = (A11 + A12)*B22 */

		AddPartialMat(AMat_5, mat1Mid, mat3Mid, A11, A12);
		strassenMatrixMulAlgo(mat1Mid, mat2Mid, mat3Mid, AMat_5, B22, M5);
		}

#pragma omp task
		{
	/* M6 = (A21 - A11)*(B11 + B12) */

		SubPartialMat(AMat_6, mat1Mid, mat3Mid, A21, A11);
		AddPartialMat(BMat_6, mat3Mid, mat2Mid, B11, B12);
		strassenMatrixMulAlgo(mat1Mid, mat2Mid, mat3Mid, AMat_6, BMat_6, M6);
		}

#pragma omp task
		{
	/* M7 = (A12 - A22)*(B21 + B22) */
		SubPartialMat(AMat_7, mat1Mid, mat3Mid, A12, A22);
		AddPartialMat(BMat_7, mat3Mid, mat2Mid, B21, B22);
		strassenMatrixMulAlgo(mat1Mid, mat2Mid, mat3Mid, AMat_7, BMat_7, M7);
		}
#pragma omp taskwait

		for (int i = 0; i < mat1Mid; i++)
			for (int j = 0; j < mat2Mid; j++) {
				C11[i][j] = M1[i][j] + M4[i][j] - M5[i][j] + M7[i][j];
				C12[i][j] = M3[i][j] + M5[i][j];
				C21[i][j] = M2[i][j] + M4[i][j];
				C22[i][j] = M1[i][j] - M2[i][j] + M3[i][j] + M6[i][j];
			}

		Free2DArray(M1);
		Free2DArray(M2);
		Free2DArray(M3);
		Free2DArray(M4);
		Free2DArray(M5);
		Free2DArray(M6);
		Free2DArray(M7);

		Free2DArray(AMat_1);
		Free2DArray(BMat_1);
		Free2DArray(AMat_2);
		Free2DArray(BMat_3);
		Free2DArray(BMat_4);
		Free2DArray(AMat_5);
		Free2DArray(AMat_6);
		Free2DArray(BMat_6);
		Free2DArray(AMat_7);
		Free2DArray(BMat_7);

		delete[] A11; delete[] A12; delete[] A21; delete[] A22;
		delete[] B11; delete[] B12; delete[] B21; delete[] B22;
		delete[] C11; delete[] C12; delete[] C21; delete[] C22;
	}
}
              
void matMultiplicationStrassen(int aSize, int bSize, int cSize, int **A, int **B, int **C) {   
#pragma omp parallel 
  {
#pragma omp single
	  {
    strassenMatrixMulAlgo(aSize, bSize, cSize, A, B, C);
	  }
  }
}  

bool checkAlgoCorrectness(int size, int **D, int **C){
			bool algoCorrectness;
	        for (int i = 0; i < size; ++i){
            for (int j = 0; j < size; ++j){
                if(C[i][j] != D[i][j]){ return false;
				}
			}
			}
			return true;
}

  
int main(int argc, char* argv[]) {   

  int k = atoi(argv[1]);

  int size = pow(2,k);

  Size = size;

  Threshold = atoi(argv[2]);  
  int proceses_count = atoi(argv[3]);
  

  int **A = Allocate2DArray(size, size);
  int **B = Allocate2DArray(size, size);
  int **C = Allocate2DArray(size, size);
  int **D = Allocate2DArray(size, size);

  int i, j;

 double startTime, endTime, startTimeSerial, endTimeSerial;
  initializeMatrices(A, B, size);

// startTimeSerial = omp_get_wtime();
// matrixMultiplication(0, size, 0, size, 0, size, A, B, D);
//   endTimeSerial = omp_get_wtime();

//     double secsForSerialMulti = endTimeSerial - startTimeSerial;

omp_set_dynamic(0);     // Explicitly disable dynamic teams
omp_set_num_threads(proceses_count);

startTime = omp_get_wtime();

  matMultiplicationStrassen(size, size, size, A, B, C);
          
  endTime = omp_get_wtime();

  double secsForParallelMulti = endTime - startTime;
  


        bool algoCorrectness;

		algoCorrectness=checkAlgoCorrectness(size, C, A);
		

        // if(algoCorrectness){
        //     printf("Correct Output: size (n=2^k) = %d, k' = %d, Threads = %d, Strassen Exec time = %lf sec, Serial Exec time = %lf sec \n",  size, Threshold, omp_get_max_threads(), secsForParallelMulti, secsForSerialMulti );
        // }
        // else{
        //     printf("Incorrect Output: size (n=2^k) = %d, k' = %d, Threads count = %d, Strassen Exec time = %lf sec, Serial Exec time = %lf sec \n",  size, Threshold, omp_get_max_threads(), secsForParallelMulti, secsForSerialMulti );
        // }


        if(algoCorrectness){
            printf("Correct Output: size (n=2^k) = %d, k' = %d, Threads = %d, Strassen Algo execution time = %lf sec\n",  size, Threshold, omp_get_max_threads(), secsForParallelMulti);
        }
        else{
            printf("Incorrect Output: size (n=2^k) = %d, k' = %d, Threads count = %d, Strassen Algo Exec time = %lf sec\n",  size, Threshold, omp_get_max_threads(), secsForParallelMulti);
        }


  Free2DArray(A);
  Free2DArray(B);
  Free2DArray(C);

  return 0;   
}  