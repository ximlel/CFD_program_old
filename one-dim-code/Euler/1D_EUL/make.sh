#!/bin/bash

ulimit -s  102400


INITIAL_PATH=$PWD

#Compile

cd ../../lib/file_io/
gcc -c ./comman.c
gcc -c ./_1D_f_io.c -I ../
ar crv file_io.a comman.o _1D_f_io.o

cd ../Riemann_solver/
gcc -c ./Riemann_solver_exact.c
gcc -c ./ROE_solver.c
ar crv Riemann_solver.a Riemann_solver_exact.o ROE_solver.o


cd $INITIAL_PATH/finite_volume_solver/
gcc -c ./linear_GRP_solver_Edir.c -I ../
gcc -c ./first_order_solver.c -I ../
ar crv finite_volume_solver.a linear_GRP_solver_Edir.o first_order_solver.o
ranlib finite_volume_solver.a

cd ../
gcc -c ./EUL_source.c
gcc -o EUL_source.out ./EUL_source.o ../../lib/file_io/file_io.a ./finite_volume_solver/finite_volume_solver.a ../../lib/Riemann_solver/Riemann_solver.a -lm

#Initilize

cd ../../../data_in/one-dim/
./data_initilize.sh
cd $INITIAL_PATH

#Run

## first order
./EUL_source.out Sod Sod Riemann_exact -1 0.43
./EUL_source.out Shock_Contact Shock_Contact Riemann_exact -1 0.49

./EUL_source.out Sod Sod_ROE ROE -1 0.43
./EUL_source.out Shock_Contact Shock_Contact_ROE ROE -1 0.49

./EUL_source.out Test1 Test1 Riemann_exact -1 0.49
./EUL_source.out Test1 Test1_ROE ROE -1 0.49

./EUL_source.out Contact_only Contact_only Riemann_exact -1 0.49
./EUL_source.out Contact_only Contact_only_ROE ROE -1 0.49


exit 0
