#!/bin/bash
make clean
make all

# Please make sure you've set the ENV variables on the terminal running the scripts
rm -rf results_openmp.txt
rm -rf results_mpi.txt
{

  echo "Edit Distance Calculation for 2 strings with OpenMP\n"
  export OMP_NUM_THREADS=1
  echo "\nOMP_NUM_THREADS =1"
  ./ed_openmp openmp/f1.txt openmp/f2.txt
  echo "\nOMP_NUM_THREADS =2"
  export OMP_NUM_THREADS=2
  ./ed_openmp openmp/f1.txt openmp/f2.txt
  echo "\nOMP_NUM_THREADS =3"
  export OMP_NUM_THREADS=3
  ./ed_openmp openmp/f1.txt openmp/f2.txt
  echo "\nOMP_NUM_THREADS =4"
  export OMP_NUM_THREADS=4
  ./ed_openmp openmp/f1.txt openmp/f2.txt
  echo "\nOMP_NUM_THREADS =8"
  export OMP_NUM_THREADS=8
  ./ed_openmp openmp/f1.txt openmp/f2.txt
  echo "\nOMP_NUM_THREADS =16"
  export OMP_NUM_THREADS=16
  ./ed_openmp openmp/f1.txt openmp/f2.txt
} > results_openmp.txt

{

  echo "Edit Distance Calculation for 2 strings with MPI\n"

  echo "\nNP =4"
  mpirun -np 4 ./ed_mpi mpi/f1.txt mpi/f2.txt

  echo "\nNP =8"
  mpirun -np 8 ./ed_mpi mpi/f1.txt mpi/f2.txt

  echo "\nNP =16"
  mpirun -np 16 ./ed_mpi mpi/f1.txt mpi/f2.txt
} > results_mpi.txt
