#!/bin/bash

ulimit -s  102400

cd ..
make
cd ./src


#./hydrocode.out Sod_test Sod_test/Sod_test_ROE 1 1_Roe free_1D
./hydrocode.out Sod_10_test Sod_10_test/Sod_10_test_ROE 2 1_Roe Sod
#./hydrocode.out Sod_10_test2 Sod_10_test2/Sod_10_test_ROE 2 1_Roe t1

./hydrocode.out odd_even odd_even/odd_even_ROE 2 1_Roe odd_even


exit 0
