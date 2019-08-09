

gfortran -O0 -g -frecord-marker=8 -c tesmodule.f90 
gcc -O0 -g -c main.c
gcc -O0 -g -frecord-marker=8 -o tescfortran main.o tesmodule.o -lgfortran
./tescfortran

gfortran -c tesmodule.f90
gcc -c main.c
gcc -o tescfortran main.o tesmodule.o
./tescfortran


gcc -c test_sparse_matrix.c
gfortran -c sparse_matrix.f90
gcc -o main.o test_sparse_matrix.o sparse_matrix.o -lgfortran
./main.o
