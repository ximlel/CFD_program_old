#!/bin/bash

ulimit -s  102400


INITIAL_PATH=$PWD

#Compile

cd ../../lib/file_io/
gcc -c ./comman.c
gcc -c ./file_io.c -I ../
ar crv file_io.a comman.o file_io.o

cd ../custom/
gcc -c ./initialize_memory.c
gcc -c ./Gauss_elimination.c
ar crv custom.a initialize_memory.o Gauss_elimination.o

cd ../Riemann_solver/
gcc -c ./GLACE.c -I ../
gcc -c ./EUCCLHYD.c -I ../
ar crv Riemann_solver.a GLACE.o EUCCLHYD.o

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
gcc -c ./first_order_solver.c -I ../
ar crv cell_centered_scheme.a first_order_solver.o

cd ../
gcc -c ./LAG_source.c -g
gcc -o LAG_source.out ./LAG_source.o ../../lib/file_io/file_io.a ./meshing/meshing.a ./cell_centered_scheme/cell_centered_scheme.a ../../lib/Riemann_solver/Riemann_solver.a ../../lib/custom/custom.a -lm

#Initilize

cd ../../../data_in/two-dim/
./data_initilize.sh
cd $INITIAL_PATH


#RUN

./LAG_source.out CB CB_EUC Sedov_mesh EUCCLHYD 90 -0.01
./LAG_source.out CB CB_GLA Sedov_mesh GLACE 3 -0.01

./LAG_source.out Sedov Sedov_GLA Sedov_mesh GLACE -1 -0.0001
./LAG_source.out Sedov Sedov_EUC Sedov_mesh EUCCLHYD -1 0.25

./LAG_source.out Saltzman Saltzman_0.93 Saltzman_mesh EUCCLHYD 9300 -0.0001

./LAG_source.out Sod Sod Sod_mesh EUCCLHYD -1 0.25
./LAG_source.out Sod_10 Sod_10 Sod_mesh EUCCLHYD -1 0.25
./LAG_source.out Sod Sod_2material Sod_2material_mesh EUCCLHYD -1 0.25
./LAG_source.out Sod_10 Sod_10_2material Sod_2material_mesh EUCCLHYD -1 0.25

./LAG_source.out Riemann_2D3_Tria Riemann_2D3_Tria Riemann_2D3_Tria_mesh EUCCLHYD -1 0.25
./LAG_source.out Riemann_2D3_Quad Riemann_2D3_Quad Riemann_2D3_Quad_mesh EUCCLHYD -1 0.25

./LAG_source.out Blunt Blunt Blunt_mesh EUCCLHYD -1 0.1



exit 0
