#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>

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
  struct timespec start_serial, end_serial;
  double diff_serial;
  clock_gettime(CLOCK_MONOTONIC_RAW, &start_serial);

  for (i = 0; i < DP_rows; i++) {
    DP[i*DP_cols] = i;
  }
  for (j = 0; j < DP_cols; j++) {
    DP[j] = j;
  }
  for (i=1; i < DP_rows; i++) {
    for (j=1; j < DP_cols; j++) {
      //int val = (X[j] != Y[i]) ? 2 : 0;
      //std::cout << DP[(i-1)*DP_cols+j]+1 << "\t" << DP[i*DP_cols+(j-1)]+1 << "\t" << DP[(i-1)*DP_cols+(j-1)]+((X[i-1] != Y[j-1]) ? 1 : 0) << std::endl;
      DP[i*DP_cols+j] = std::min(std::min((DP[(i-1)*DP_cols+j]+1), (DP[i*DP_cols+(j-1)]+1)),(DP[(i-1)*DP_cols+(j-1)]+ ((X[j-1] != Y[i-1]) ? 1 : 0)));
    }
  }
  clock_gettime(CLOCK_MONOTONIC_RAW, &end_serial);
  diff_serial = double((end_serial.tv_sec - start_serial.tv_sec) + (end_serial.tv_nsec - start_serial.tv_nsec)/1000000000.0);

  std::cout << X << "\t" << X.length() << std::endl << Y << "\t" << Y.length() << std::endl;
  printMatrix(DP, DP_rows, DP_cols);

  printf("Serial Execution Time=%f\n",diff_serial);
  return 0;
}
