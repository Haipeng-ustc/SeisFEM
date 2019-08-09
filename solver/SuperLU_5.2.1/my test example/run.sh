#!/bin/bash

gcc superlu.c -o superlu -I/Users/haipeng/Desktop/superlu_test/SuperLU_5.2.1/SRC/ /Users/haipeng/Desktop/superlu_test/SuperLU_5.2.1/lib/libblas.a /Users/haipeng/Desktop/superlu_test/SuperLU_5.2.1/lib/libsuperlu_5.2.1.a  -lm -lgfortran


