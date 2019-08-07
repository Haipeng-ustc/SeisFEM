#! /bin/bash
#
cp fem2d_pack.h /$HOME/include
#
gcc -c -Wall -I/$HOME/include fem2d_pack.c
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv fem2d_pack.o ~/libc/fem2d_pack.o
#
echo "Normal end of execution."
