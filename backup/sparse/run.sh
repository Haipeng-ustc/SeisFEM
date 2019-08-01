#!/bin/bash

gfortran -c sparse_matrix.f90
gcc -c test_sparse_matrix.c
gcc -o main.o test_sparse_matrix.o sparse_matrix.o -lgfortran
./main.o
