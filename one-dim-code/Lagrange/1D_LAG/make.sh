#!/bin/bash

ulimit -s  102400


INITIAL_PATH=$PWD

#Compile

cd ../../lib/file_io/
gcc -c ./comman.c
gcc -c ./_1D_f_io.c -I ../
ar crv file_io.a comman.o _1D_f_io.o

cd ../Riemann_solver/
gcc -c ./Riemann_solver_exact.c -g


cd $INITIAL_PATH/finite_volume_solver/
gcc -c ./first_order_solver.c -I ../ -g
gcc -c ./linear_GRP_solver_LAG.c -I ../-g
gcc -c ./second_order_solver.c -I ../ -g
ar crv finite_volume_solver.a first_order_solver.o linear_GRP_solver_LAG.o second_order_solver.o
ranlib finite_volume_solver.a

cd ../
gcc -c ./LAG_source.c -g
gcc -o LAG_source.out ./LAG_source.o ../../lib/file_io/file_io.a ./finite_volume_solver/finite_volume_solver.a ../../lib/Riemann_solver/Riemann_solver_exact.o -lm

#Initilize

cd ../../../data_in/one-dim/
./data_initilize.sh
cd $INITIAL_PATH

#Run

./LAG_source.out Sod Sod GRP -1 0.43 second_order
./LAG_source.out Shock_Contact Shock_Contact GRP -1 0.49 second_order

./LAG_source.out Sod Sod Riemann_exact -1 0.43
./LAG_source.out Shock_Contact Shock_Contact Riemann_exact -1 0.49


exit 0
