#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <omp.h>

void printMatrix(int* matrix, int rows, int cols) {
  for (int i=0; i < rows; i++) {
    for (int j=0; j < cols; j++) {
      std::cout << matrix[i*cols+j] << "\t";
    }
    std::cout << std::endl;
  }
}
void initializeMatrix(int* &matrix, int rows, int cols) {
  for (int i=0; i < rows; i++) {
    for (int j=0; j < cols; j++) {
      matrix[i*cols+j] = 100;
    }
  }
}

int main(int argc, char* argv[]) {
  std::string X = argv[1];
  std::string Y = argv[2];
  int N = X.length();
  int M = Y.length();
  int DP_rows = M+1; int DP_cols = N+1;
  int* DP = new int [DP_rows*DP_cols];
  int i,j;
  initializeMatrix(DP, DP_rows, DP_cols);
  double start, end, diff_parallel;
  start = omp_get_wtime();

  for (i = 0; i < DP_rows; i++) {
    DP[i*DP_cols] = i;
  }
  for (j = 0; j < DP_cols; j++) {
    DP[j] = j;
  }

  for (i=1; i < DP_rows; i++) {
    for (j=1; j < DP_cols; j++) {
      DP[i*DP_cols+j] = std::min(std::min((DP[(i-1)*DP_cols+j]+1), (DP[i*DP_cols+(j-1)]+1)),(DP[(i-1)*DP_cols+(j-1)]+ ((X[i-1] != Y[j-1]) ? 1 : 0)));
    }
  }
  end = omp_get_wtime();
	diff_parallel = end - start;
  std::cout << X << "\t" << X.length() << std::endl << Y << "\t" << Y.length() << std::endl;
  printMatrix(DP, DP_rows, DP_cols);
  printf("Parallel Execution Time=%f\n",diff_parallel);

  return 0;
}
