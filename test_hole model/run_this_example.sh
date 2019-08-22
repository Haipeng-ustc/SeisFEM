#!/bin/bash

rm *.o
rm *.dat
gfortran -c coo2csr_lib.f90
gcc -fopenmp -c hole_model.c -lm
gcc  -o main.o hole_model.o coo2csr_lib.o -L /Users/haipeng/Desktop/seisfem/solver/pardiso/lib/ -lpardiso600-MACOS-X86-64 -I/Users/haipeng/Desktop/seisfem/solver/SuperLU_5.2.1/SRC/ /Users/haipeng/Desktop/seisfem/solver/SuperLU_5.2.1/lib/libsuperlu_5.2.1.a /Users/haipeng/Desktop/seisfem/solver/SuperLU_5.2.1/lib/libblas.a -fopenmp -lgfortran -lm
#  -lgfortran and -lm link math packages

if [ $? -ne 0 ]; then
echo "Compile error."
exit
fi
./main.o

echo "Normal end of execution."
