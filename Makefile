CC = gcc									#C compiler
CFLAGS += 
AR = ar										#Creat static libs
ARFLAGS = crv								#
INC = ./src									#Include
DIC = `ls -d ./src/*/`						#Subdirectories

all : modules
	$(AR)


.PHONY all

modules:						
	for n in $(DIC) 						#Enter each subdirectory
	do
		exit = make --directory=$$n			#Call the Makefile in the subdirectory

		if[ $exit = "1" ];then				#Check the status of the exit code in the Makefile.
			exit
		fi
	done

./PHONY modules

clean:
	for n in $(DIC) 					
	do
		make --directory=$$n clean			#Call the Makefile in the subdirectory, and run clean.
	done

./PHONY clean

ar crv cell_centered_scheme.a linear_GRP_solver_Edir.o linear_GRP_solver_Edir_2D.o slope_limiter_Ven.o first_order_solver.o first_order_two_species_solver.o second_order_solver.o second_order_two_species_solver.o
ranlib cell_centered_scheme.a


$(CC) -c ./EUL_source.c
$(CC) -o EUL_source.out ./EUL_source.o ../../lib/file_io/file_io.a ./meshing/meshing.a ./cell_centered_scheme/cell_centered_scheme.a ../../lib/Riemann_solver/Riemann_solver.a ../../lib/custom/custom.a -lmi


clean:
