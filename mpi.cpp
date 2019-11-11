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

void distributeData(std::string X, std::string Y, std::string &pX, int* &DP_proc, int DP_rows, int DP_cols, int &pNumCols) {
  MPI_Status status;
  MPI_Bcast(&DP_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&DP_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Y[0], DP_rows, MPI_CHAR, 0, MPI_COMM_WORLD);
  int* pBlockSize = new int [pNum];
  int* pBlockInd = new int [pNum];
  pBlockSize[0] = BLOCK_SIZE(0, pNum, DP_cols);
  pBlockInd[0] = 0;
  for (int i =1; i < pNum; i++) {
    pBlockInd[i] = pBlockInd[i-1] + pBlockSize[i-1];
    pBlockSize[i] = BLOCK_SIZE(i, pNum, DP_cols);
  }
  pNumCols = BLOCK_SIZE(pRank, pNum, DP_cols);
  char* pX_char = new char [pBlockSize[pRank]];

  if (pRank == 1) {
    cout << Y[498] << "\t" << DP_rows << "\t" << DP_cols<< endl;
  }
  if (pRank == 0) {
    for (int j = 0; j < pNum; j++) {
      MPI_Send(&X[pBlockInd[j]], pBlockSize[j], MPI_CHAR, j, j, MPI_COMM_WORLD);
    }
  }
  MPI_Recv(&pX_char[0], pBlockSize[pRank], MPI_CHAR, 0, pRank, MPI_COMM_WORLD, &status);
  pX = std::string(pX_char);
  //if (pRank == 0) cout << pX_char[0] << endl;
  //if (pRank == 1) cout << pX_char[0] << endl;
  DP_proc = new int [DP_rows*pNumCols];
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
  std::string pX;
  int N, M, DP_cols, DP_rows, pNumCols;
  int* DP;
  int* DP_proc;
  MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &pNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &pRank);
  if (pRank == 0) {
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
    N = X.length();
    M = Y.length();
    DP_rows = M+1;
    DP_cols = N+1;
    DP = new int [DP_rows*DP_cols];
  }
  distributeData(X, Y, pX, DP_proc, DP_rows, DP_cols, pNumCols);
  if (pRank == 0) cout << pNumCols << endl;
  //if (pRank == 1) cout << pX[0] << endl;
  //if (pRank==0) cout << pRank << "\t" << N << "\t" << M << endl;
  MPI_Finalize();
  return 0;
}
