#!/bin/bash

ulimit -s  10240


INITIAL_PATH=$PWD

#Compile

cd ../../lib/file_io/
gcc -c ./comman.c
gcc -c ./file_io.c -I ../
ar crv file_io.a comman.o file_io.o

cd ../custom/
gcc -c ./initialize_memory.c
gcc -c ./miu.c
ar crv custom.a initialize_memory.o miu.o

cd ../Riemann_solver/
gcc -c ./GLACE.c -I ../
gcc -c ./EUCCLHYD.c -I ../
gcc -c ./EUCCLHYD_dt.c -I ../
gcc -c ./GLACE_dt.c -I ../
gcc -c ./EUCCLHYD_iteration.c -I ../
gcc -c ./EUCCLHYD_dt_iteration.c -I ../
ar crv Riemann_solver.a GLACE.o EUCCLHYD.o EUCCLHYD_dt.o GLACE_dt.o EUCCLHYD_iteration.o EUCCLHYD_dt_iteration.o

cd $INITIAL_PATH/meshing/
gcc -c ./Sedov_mesh.c -I ../
gcc -c ./Sod_mesh.c -I ../
gcc -c ./Sod_2material_mesh.c -I ../
gcc -c ./Saltzman_mesh.c -I ../
gcc -c ./Riemann_2D3_Tria_mesh.c -I ../
gcc -c ./Riemann_2D3_Quad_mesh.c -I ../
gcc -c ./Blunt_mesh.c -I ../
ar crv meshing.a Sedov_mesh.o Sod_mesh.o Sod_2material_mesh.o Saltzman_mesh.o Riemann_2D3_Tria_mesh.o Riemann_2D3_Quad_mesh.o Blunt_mesh.o

cd ../cell_centered_scheme/
gcc -c ./second_order_solver.c -I ../
gcc -c ./second_order_solver_GLACE.c -I ../
gcc -c ./second_order_iteration_solver.c -I ../
ar crv cell_centered_scheme.a second_order_solver.o second_order_solver_GLACE.o second_order_iteration_solver.o

cd ../
gcc -c ./LAG_source.c
gcc -o LAG_source.out ./LAG_source.o ../../lib/file_io/file_io.a ./meshing/meshing.a ./cell_centered_scheme/cell_centered_scheme.a ../../lib/Riemann_solver/Riemann_solver.a ../../lib/custom/custom.a -lm

#Initilize

cd ../../../data_in/two-dim/
./data_initilize.sh
cd $INITIAL_PATH

#RUN

./LAG_source.out Sedov Sedov_EUC Sedov_mesh nonlinear Ven -1 0.25

./LAG_source.out Saltzman Saltzman_0.65_GLA Saltzman_mesh GLACE Ven 6500 -0.0001
./LAG_source.out Saltzman Saltzman_0.95 Saltzman_mesh linear Ven 9500 -0.0001
./LAG_source.out Saltzman Saltzman_0.99 Saltzman_mesh linear Ven 99000 -0.00001

./LAG_source.out Sod Sod Sod_mesh linear BJ -1 0.25
./LAG_source.out Sod_10 Sod_10 Sod_mesh linear BJ -1 0.25
./LAG_source.out Sod Sod_2material Sod_2material_mesh linear BJ -1 0.25
./LAG_source.out Sod_10 Sod_10_2material Sod_2material_mesh linear BJ -1 0.25

./LAG_source.out Riemann_2D3_Tria Riemann_2D3_Tria Riemann_2D3_Tria_mesh linear BJ -1 0.25
./LAG_source.out Riemann_2D3_Quad Riemann_2D3_Quad Riemann_2D3_Quad_mesh linear BJ -1 0.25

./LAG_source.out Blunt Blunt Blunt_mesh linear Ven -1 0.1



exit 0
