#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <mpi.h>
#include <fstream>

#define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW(id+1, p, n) - 1 )
#define BLOCK_SIZE(id, p, n) (BLOCK_LOW(id+1, p, n) - BLOCK_LOW(id, p, n))
#define BLOCK_OWNER(index, p, n) (((p) * ((index)+1)-1)/(n))

using namespace std;
int pNum=0, pRank=0;

void printMatrix(int* matrix, int rows, int cols) {
  for (int i=0; i < rows; i++) {
    for (int j=0; j < cols; j++) {
      std::cout << matrix[i*cols+j] << "\t";
    }
    std::cout << std::endl;
  }
}
void processTerminate(int* DP, int* DP_proc, int* pBlockInd, int* pBlockSize){
	if (pRank == 0){delete [] DP;}
	delete [] DP_proc;
	delete [] pBlockInd;
	delete [] pBlockSize;
}

void distributeData(std::string X, std::string &Y, std::string &pX, int* &DP_proc, int &M, int &N, int &pNumCols, int* &pBlockInd, int* &pBlockSize) {
  MPI_Status status;
  MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Y[0], M, MPI_CHAR, 0, MPI_COMM_WORLD);
  pBlockSize = new int [pNum];
  pBlockInd = new int [pNum];
  pBlockSize[0] = BLOCK_SIZE(0, pNum, N);
  pBlockInd[0] = 0;
  for (int i =1; i < pNum; i++) {
    pBlockInd[i] = pBlockInd[i-1] + pBlockSize[i-1];
    pBlockSize[i] = BLOCK_SIZE(i, pNum, N);
  }
  pNumCols = BLOCK_SIZE(pRank, pNum, N);
  char* pX_char = new char [pBlockSize[pRank]];
  if (pRank == 0) {
    for (int j = 0; j < pNum; j++) {
      MPI_Send(&X[pBlockInd[j]], pBlockSize[j], MPI_CHAR, j, j, MPI_COMM_WORLD);
    }
  }
  MPI_Recv(&pX_char[0], pBlockSize[pRank], MPI_CHAR, 0, pRank, MPI_COMM_WORLD, &status);
  pX = std::string(pX_char);
  DP_proc = new int [(M+1)*(pBlockSize[pRank]+1)];
  delete [] pX_char;
}

void createMatrix(std::string X, std::string &Y, std::string pX, int* DP_proc, int M, int N, int pNumCols, int* pBlockInd, int* pBlockSize) {
  int i, j;
  MPI_Status status;
  MPI_Request *request;
  request = (MPI_Request*)calloc(M, sizeof(MPI_Request));
  pNumCols+=1;
  if (pRank == 0) {
    for (i = 0; i < M+1; i++) {
      DP_proc[i*pNumCols] = i;
    }
  }
  for (j = 0; j < pNumCols ; j++) {
    DP_proc[j] = j+pBlockInd[pRank];
  }
  for (i=1; i < M+1; i++) {
     if (pRank != 0) {
       MPI_Recv(&DP_proc[i*pNumCols+0], 1, MPI_INT, pRank-1, pRank, MPI_COMM_WORLD, &status);
     }
     for (j = 1; j < pNumCols; j++) {
       DP_proc[i*pNumCols+j] = std::min(std::min((DP_proc[(i-1)*pNumCols+j]+1), (DP_proc[i*pNumCols+(j-1)]+1)),(DP_proc[(i-1)*pNumCols+(j-1)]+ ((pX[j-1] != Y[i-1])? 1 : 0)));
     }
     if (pRank < pNum-1) {
       MPI_Isend(&DP_proc[i*pNumCols+(j-1)], 1, MPI_INT, pRank+1, pRank+1, MPI_COMM_WORLD, &request[i-1]);
     }
   }
}

void gatherResult(int* DP_proc, int* DP, int M, int N, int* pBlockInd, int* pBlockSize, int pNumCols) {
  int i, j;
  pBlockSize[0]+=1;
  int* row;
  for (i=1; i < pNum; i++) {
    pBlockInd[i] = pBlockInd[i-1] + pBlockSize[i-1];
  }
  for (i=0; i < M+1; i++) {
    int* newRow = new int [N+1];
    if(pRank == 0) row = new int [pNumCols+1];
    if(pRank != 0) row = new int [pNumCols];
    if (pRank == 0) {
      for (j=0; j < pNumCols+1; j++) {
        row[j] = DP_proc[i*(pNumCols+1)+j];
      }
    }
    if (pRank) {
      for (j=0; j < pNumCols; j++) {
        row[j] = DP_proc[i*(pNumCols+1)+j+1];
      }
    }

    MPI_Gatherv(row, pBlockSize[pRank], MPI_INT, newRow, pBlockSize, pBlockInd, MPI_INT, 0, MPI_COMM_WORLD);
    for (j=0; j < N+1; j++) {
      if (pRank == 0) DP[i*(N+1)+j] = newRow[j];
    }
    delete [] row;
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
  else printf("%s\n", "Good Job!");
}

int main(int argc, char* argv[]) {
  std::string X = "";
  std::string Y = "";
  std::string pX= "";
  int N, M, pNumCols;
  int* DP;
  int* DP_proc;
  int *pBlockSize, *pBlockInd;
  MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &pNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &pRank);
  if (pRank == 0) {
    X = "Hellowerssssajhk";
    Y= "Yellowersssssjib";
    N = X.length();
    M = Y.length();
    DP = new int [(M+1)*(N+1)];
  }
  double start, end, diff_parallel;
  start = MPI_Wtime();
  distributeData(X, Y, pX, DP_proc, M, N, pNumCols, pBlockInd, pBlockSize);
  createMatrix(X, Y, pX, DP_proc, M, N, pNumCols, pBlockInd, pBlockSize);
  gatherResult(DP_proc, DP, M, N, pBlockInd, pBlockSize, pNumCols);

  end = MPI_Wtime();
  if (pRank == 0) testResult(X, Y, DP, M+1, N+1);
  diff_parallel = end - start;
  if (pRank == 0) printf("MPI Time = %fs\n", diff_parallel);
  processTerminate(DP, DP_proc, pBlockInd, pBlockSize);
  MPI_Finalize();
  return 0;
}
