#!/bin/bash

rm *.o
rm *.dat
gfortran -c coo2csr_lib.f90
gcc -fopenmp -c test_all.c -lm
gcc -fopenmp -o main.o test_all.o coo2csr_lib.o -lgfortran -lm  #  -lgfortran and -lm link packages

if [ $? -ne 0 ]; then
echo "Compile error."
exit
fi
./main.o

echo "Normal end of execution."
