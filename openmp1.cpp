#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <omp.h>
#include <fstream>

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

void testResult(std::string X, std::string Y, int* DP, int DP_rows, int DP_cols) {
  int i, j;
  int* serialMatrix= new int [DP_rows*DP_cols];
  for (i = 0; i < DP_rows; i++) {
    serialMatrix[i*DP_cols] = i;
  }
  for (j = 0; j < DP_cols; j++) {
    serialMatrix[j] = j;
  }
  for (i=1; i < DP_rows; i++) {
    for (j=1; j < DP_cols; j++) {
      serialMatrix[i*DP_cols+j] = std::min(std::min((serialMatrix[(i-1)*DP_cols+j]+1),
                                                    (serialMatrix[i*DP_cols+(j-1)]+1)),
                                          (serialMatrix[(i-1)*DP_cols+(j-1)]+ ((X[j-1] != Y[i-1]) ? 1 : 0)));
    }
  }
  int flag=0;
  for (i=0; i < DP_rows; i++) {
    for (j=0; j < DP_cols; j++) {
      if (serialMatrix[i*DP_cols+j] != DP[i*DP_cols+j]) flag = -1;
    }
  }
  if (flag == -1) printf("ERROR!!!! Please check the code. The serial and parallel DP matrix are not the same.\n");
}

int main(int argc, char* argv[]) {
  std::string X = "";
  std::string Y = "";
  std::ifstream ifile1;
  ifile1.open(argv[1]);
	if(!ifile1)
	{
		std::cout<<"Error in opening file..!!";
		exit(0);
	}
	while(ifile1.eof()==0)
	{
        ifile1>>X;
	}
  std::ifstream ifile2;
  ifile2.open(argv[2]);
	if(!ifile2)
	{
		std::cout<<"Error in opening file..!!";
		exit(0);
	}
	while(ifile2.eof()==0)
	{
        ifile2>>Y;
	}
  //std::cout << X.size() << std::endl;
  int N = X.length();
  int M = Y.length();
  int DP_rows = M+1; int DP_cols = N+1;
  int* DP = new int [DP_rows*DP_cols];
  int i,j, it, k, w1, w2;
  initializeMatrix(DP, DP_rows, DP_cols);
  double start, end, diff_parallel;
  start = omp_get_wtime();

  for (i = 0; i < DP_rows; i++) {
    DP[i*DP_cols] = i;
  }
  for (j = 0; j < DP_cols; j++) {
    DP[j] = j;
  }
  int P=0, tid=0;
  #pragma omp parallel firstprivate(DP_rows, DP_cols, X, Y, P) shared(DP) private(it, tid, i, j)
  {
  P = omp_get_num_threads();
  for (it=1; it < (DP_cols+DP_rows); it++){
    tid = omp_get_thread_num();
    for (j=tid+1; j<=it; j+=P) {
      i = it-j+1;
      if (i >= DP_rows || j >= DP_cols) continue;
      DP[i*DP_cols+j] = std::min(std::min((DP[(i-1)*DP_cols+j]+1), (DP[i*DP_cols+(j-1)]+1)),(DP[(i-1)*DP_cols+(j-1)]+ ((X[j-1] != Y[i-1]) ? 1 : 0)));
    }
    #pragma omp barrier
  }
}
  end = omp_get_wtime();
	diff_parallel = end - start;
  testResult(X, Y, DP, DP_rows, DP_cols);
  //std::cout << X << "\t" << X.length() << std::endl << Y << "\t" << Y.length() << std::endl;
  //printMatrix(DP, DP_rows, DP_cols);
  printf("Parallel Execution Time=%f\n",diff_parallel);

  return 0;
}
