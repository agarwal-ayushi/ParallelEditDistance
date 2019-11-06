#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>

int main(int argc, char* argv[]) {
  std::string X = "Hello";
  std::string Y = "Yellow";
  int N = X.length();
  int M = Y.length();
  int* DP = new int(N*M);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j <M; j++) {


    }


  }

  std::cout << X << "\t" << X.length() << std::endl << Y << "\t" << Y.length() << std::endl;
  return 0;
}
