#! /bin/bash
#
cp mgmres.h /$HOME/include
#
gcc -c -Wall mgmres.c
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv mgmres.o ~/libc/mgmres.o
#
echo "Normal end of execution."
