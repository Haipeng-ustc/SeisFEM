#!/bin/bash

gfortran -c coo2csr_lib.f90
gcc -c test_coo2csr_lib.c
gcc -o main.o test_coo2csr_lib.o coo2csr_lib.o -lgfortran
./main.o
