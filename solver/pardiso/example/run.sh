#!/bin/bash
export PARDISOLICMESSAGE=1
export OMP_NUM_THREADS=1
gcc 
gcc pardiso_main.c -o test -L /Users/haipeng/Desktop/seisfem/solver/pardiso/lib/ -lpardiso600-MACOS-X86-64 -I/Users/haipeng/Desktop/seisfem/solver/SuperLU_5.2.1/SRC/ /Users/haipeng/Desktop/seisfem/solver/SuperLU_5.2.1/lib/libsuperlu_5.2.1.a /Users/haipeng/Desktop/seisfem/solver/SuperLU_5.2.1/lib/libblas.a -fopenmp  -lm -lgfortran
