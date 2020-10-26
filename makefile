fftw = /Users/yang/Downloads/fftw-3.3.6-pl2/

default:
	mpicc -Wall -O2 -c pfcpure3d.c -I$(fftw)/include 
	mpicc -Wall -lm -o pfcpure3d pfcpure3d.o -L/usr/local/lib -lfftw3 -lfftw3_mpi

