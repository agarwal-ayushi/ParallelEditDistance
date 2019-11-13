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
  cout.flush();
}
void processTerminate(int* DP, int* DP_proc, int* pBlockInd, int* pBlockSize, char* pX_char){
	if (pRank == 0){delete [] DP;}
	delete [] DP_proc;
	delete [] pBlockInd;
	delete [] pBlockSize;
  delete [] pX_char;
}

void distributeData(std::string &X, std::string &Y, char* &pX_char, int* &DP_proc, int &M, int &N, int &pNumCols, int* &pBlockInd, int* &pBlockSize) {
  MPI_Status status;
  int recv_amount;
  MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Y[0], M, MPI_CHAR, 0, MPI_COMM_WORLD);

  //if (pRank == 7) cout << Y[300] << "\n\n\n";
  pBlockSize = new int [pNum];
  pBlockInd = new int [pNum];
  pBlockSize[0] = BLOCK_SIZE(0, pNum, N);
  pBlockInd[0] = 0;
  for (int i =1; i < pNum; i++) {
    pBlockInd[i] = pBlockInd[i-1] + pBlockSize[i-1];
    pBlockSize[i] = BLOCK_SIZE(i, pNum, N);
  }
  pNumCols = BLOCK_SIZE(pRank, pNum, N);
  pX_char= new char [pBlockSize[pRank]];

  if (pRank == 0) {
    for (int j = 0; j < pNum; j++) {
      MPI_Send(&X[pBlockInd[j]], pBlockSize[j], MPI_CHAR, j, j, MPI_COMM_WORLD);
    }
  }
  MPI_Recv(&pX_char[0], pBlockSize[pRank], MPI_CHAR, 0, pRank, MPI_COMM_WORLD, &status);
  MPI_Get_count(&status, MPI_CHAR, &recv_amount);
  //cout << status.MPI_SOURCE << status.MPI_TAG << "\t" << number_amount<< "\n\n\n\n";
  DP_proc = new int [(M+1)*(pBlockSize[pRank]+1)];
}

void createMatrix(std::string X, std::string &Y, char* &pX_char, int* &DP_proc, int M, int N, int pNumCols, int* pBlockInd, int* pBlockSize) {
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
       DP_proc[i*pNumCols+j] = std::min(std::min((DP_proc[(i-1)*pNumCols+j]+1), (DP_proc[i*pNumCols+(j-1)]+1)),(DP_proc[(i-1)*pNumCols+(j-1)]+ ((pX_char[j-1] != Y[i-1])? 1 : 0)));
     }
     if (pRank < pNum-1) {
       MPI_Isend(&DP_proc[i*pNumCols+(j-1)], 1, MPI_INT, pRank+1, pRank+1, MPI_COMM_WORLD, &request[i-1]);
     }
   }

  // if (pRank == 7) cout << i << "\t" << DP_proc[M*pNumCols+(pNumCols-1)] << "\n\n\n\n";
}

void gatherResult(int* &DP_proc, int* &DP, int M, int N, int* pBlockInd, int* pBlockSize, int pNumCols) {
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
    delete [] newRow;
  }
}
void testResult(std::string X, std::string &Y, int* &DP, int DP_rows, int DP_cols) {
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
  cout << serialMatrix[(DP_rows-1)*DP_cols+(DP_cols-1)] << "\t " << DP[(DP_rows-1)*DP_cols+(DP_cols-1)]<< "\n\n\n\n";
  int flag=0;
  for (i=0; i < DP_rows; i++) {
    for (j=0; j < DP_cols; j++) {
      if (serialMatrix[i*DP_cols+j] != DP[i*DP_cols+j]) flag++ ;
    }
  }
  if (flag != 0) printf("ERROR!!!! Please check the code. The serial and parallel DP matrix are not the same - mismatches= %d\n", flag);
  else printf("%s\n", "Good Job!");
  delete [] serialMatrix;
}

int main(int argc, char* argv[]) {
  std::string X = "";
  std::string Y = "";
  int N, M, pNumCols;
  int* DP;
  int* DP_proc;
  int *pBlockSize, *pBlockInd;
  char* pX_char;

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
    //X = "Hellower";
    //Y= "Yellow";
    N = X.length();
    M = Y.length();
    DP = new int [(M+1)*(N+1)];
  }


  double start, end, diff_parallel, average=0.0;
  MPI_Barrier(MPI_COMM_WORLD);
  for (int i=0; i < 100; i++) {
    diff_parallel=0.0;
    start = MPI_Wtime();
    distributeData(X, Y, pX_char, DP_proc, M, N, pNumCols, pBlockInd, pBlockSize);
    createMatrix(X, Y, pX_char, DP_proc, M, N, pNumCols, pBlockInd, pBlockSize);
    gatherResult(DP_proc, DP, M, N, pBlockInd, pBlockSize, pNumCols);
    end = MPI_Wtime();
    diff_parallel = end - start;
    average += diff_parallel;
  }
  average /= 100;
  if (pRank == 0) testResult(X, Y, DP, M+1, N+1);
  if (pRank == 0) printf("MPI Time = %fs\n", average);
  processTerminate(DP, DP_proc, pBlockInd, pBlockSize, pX_char);
  MPI_Finalize();
  return 0;
}
