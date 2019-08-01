#!/bin/bash


gfortran -c coo2csr_lib.f90
gcc -c test_all.c -lm
gcc -o main.o test_all.o coo2csr_lib.o -lgfortran -lm  #  -lgfortran and -lm link packages
./main.o
