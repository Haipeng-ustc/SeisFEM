#!/bin/bash

cp ../../sparse_matrix/coo2csr_lib.f90 .
cp ../../solver/pardiso/libpardiso600-MACOS-X86-64.dylib .
cp ../../solver/pardiso/pardiso.lic .
cp ../../seisfem/seisfem.c .
gfortran -c coo2csr_lib.f90
gcc -fopenmp -c seisfem.c -lm
gcc  -o seisfem seisfem.o coo2csr_lib.o -L /Users/haipeng/Desktop/seisfem/solver/pardiso/lib/ -lpardiso600-MACOS-X86-64 -I/Users/haipeng/Desktop/seisfem/solver/SuperLU_5.2.1/SRC/ /Users/haipeng/Desktop/seisfem/solver/SuperLU_5.2.1/lib/libsuperlu_5.2.1.a /Users/haipeng/Desktop/seisfem/solver/SuperLU_5.2.1/lib/libblas.a -fopenmp -lgfortran -lm
#  -lgfortran and -lm link math packages


if [ $? -ne 0 ]; then
echo "Compile error."
exit
fi

./seisfem

rm seisfem.c
rm seisfem.o
rm coo2csr_lib.o
rm coo2csr_lib.f90
rm coo2csr_lib.mod
rm libpardiso600-MACOS-X86-64.dylib
rm pardiso.lic
