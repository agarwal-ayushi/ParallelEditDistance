CC = gcc
CXX = g++
MPICXX = mpicxx
CFLAGS = -g -Wall
CFLAGS_opemp = -g -Wall -fopenmp

all: ed_serial ed_openmp ed_mpi plots

serial: ed_serial.cpp
	$(CXX) $(CFLAGS) ed_serial.cpp -o ed_serial

ed_openmp: openmp/ed_openmp.cpp
	$(CXX) $(CFLAGS_opemp) openmp/ed_openmp.cpp -o ed_openmp

ed_mpi: mpi/mpi_edit_distance.cpp
	$(MPICXX) mpi/mpi_edit_distance.cpp -o ed_mpi

clean:
	rm ed_serial ed_openmp ed_mpi

plots:
	gnuplot -e "set terminal jpeg; set xlabel 'No of threads'; set ylabel 'Time'; set title 'Performance comparison'; set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5; set pointintervalbox 3; plot 'results_openmp.txt' with linespoints ls 1" > openmp.jpeg
	eog openmp.jpeg
	gnuplot -e "set terminal jpeg; set xlabel 'No of threads'; set ylabel 'Time'; set title 'Performance comparison'; set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5; set pointintervalbox 3; plot 'results_mpi.txt' with linespoints ls 1" > mpi.jpeg
	eog mpi.jpeg
	gnuplot -e "set terminal jpeg; set xlabel 'No of threads'; set ylabel 'Time'; set title 'Performance comparison'; set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5; set pointintervalbox 3; plot 'results_cuda.txt' with linespoints ls 1" > cuda.jpeg
	eog cuda.jpeg
