CC = gcc
CXX = g++
MPICXX = mpicxx
CFLAGS = -g -Wall
CFLAGS_opemp = -g -Wall -fopenmp

all: ed_serial ed_openmp ed_mpi

serial: ed_serial.cpp
	$(CXX) $(CFLAGS) ed_serial.cpp -o ed_serial

ed_openmp: openmp/ed_openmp.cpp
	$(CXX) $(CFLAGS_opemp) openmp/ed_openmp.cpp -o ed_openmp

ed_mpi: mpi/mpi_edit_distance.cpp
	$(MPICXX) mpi/mpi_edit_distance.cpp -o ed_mpi

clean:
	rm ed_serial ed_openmp ed_mpi
