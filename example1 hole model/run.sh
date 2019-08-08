#!/bin/bash

rm *.o
rm *.dat
gfortran -c coo2csr_lib.f90
gcc -fopenmp -c hole_model.c -lm
gcc -fopenmp -o main.o hole_model.o coo2csr_lib.o -lgfortran -lm  #  -lgfortran and -lm link packages

if [ $? -ne 0 ]; then
echo "Compile error."
exit
fi
./main.o

echo "Normal end of execution."
