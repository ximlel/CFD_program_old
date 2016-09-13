/* 
 * This is a implementation of 1-D or 2-D hydrocodes in Cartesian geometry
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <dirent.h>

#include "./include/var_struc.h"
#include "./include/tools.h"
//#include "./include/Riemann_solver.h"
#include "./include/file_io.h"
#include "./include/meshing.h"
#include "./include/finite_volume.h"



/*The global primitive variable and configuration array.
 */

double config[N_CONF];



int main(int argc, char *argv[])
{
	printf("\n%s\n", argv[2]);

	for (int i = 0; i < N_CONF; i++)
			config[i] = 1.0/0.0;

	// Set dimension.
	config[0] = (double)atoi(argv[3]);

	char * scheme;
	config[9] = (double)strtol(argv[4], &scheme, 10);	
	if (* scheme == '_')
		scheme++;
	else
		{
			printf("No order!\n");
			exit(2);
		}

	struct flu_var FV = flu_conf_load(argv[1]);

	struct mesh_var mv= mesh_load(argv[1], argv[5]);
	
	/*  
	double cpu_time[N];
	cpu_time[0] = 0.0;

	double X[NUM_POINT];
	double Y[NUM_POINT];
  
                                            
	int * CELL_POINT[NUM_CELL];
  
	int * BOUNDARY_POINT[2];

	double gamma[NUM_CELL];
	*/

	printf("Output: %s.\n", argv[2]);

	/*
	if(strcmp(argv[3],"Sod_mesh")==0)
		NUM_BOUNDARY = Sod_mesh(CELL_POINT, X, Y, BOUNDARY_POINT, gamma, config, m, n);
	else if(strcmp(argv[3],"odd_even_mesh")==0)
		NUM_BOUNDARY = odd_even_mesh(CELL_POINT, X, Y, BOUNDARY_POINT, gamma, config, m, n);
	else if(strcmp(argv[3],"odd_even_EW_upstream_mesh")==0)
		NUM_BOUNDARY = odd_even_EW_upstream_mesh(CELL_POINT, X, Y, BOUNDARY_POINT, gamma, config, m, n);
	else if(strcmp(argv[3],"Free_mesh")==0)
		NUM_BOUNDARY = Free_mesh(CELL_POINT, X, Y, BOUNDARY_POINT, gamma, config, m, n);
	else if(strcmp(argv[3],"Shear_mesh")==0)
		NUM_BOUNDARY = Shear_mesh(CELL_POINT, X, Y, BOUNDARY_POINT, gamma, config, m, n);
	else if(strcmp(argv[3],"Cylinder_mesh")==0)
		NUM_BOUNDARY = Cylinder_mesh(CELL_POINT, X, Y, BOUNDARY_POINT, gamma, config, m, n);
	else if(strcmp(argv[3],"RMI_mesh")==0)
		NUM_BOUNDARY = RMI_mesh(CELL_POINT, X, Y, BOUNDARY_POINT, gamma, config, m, n);
	else
		printf("No mesh setting!\n");

	int STEP;	
	if(atoi(argv[5])>=0)
		STEP = atoi(argv[5]);
	else
		STEP = N;
	


	if(strcmp(argv[argc-1],"second_order")==0||strcmp(argv[argc-2],"second_order")==0)
		{			
			printf("second order\n");
			if(strcmp(argv[argc-1],"Two_species")==0)
				{
					printf("Two species\n");		  
					second_order_two_species_solver(&STEP, config, NUM_CELL, NUM_POINT, NUM_BOUNDARY, CELL_POINT, BOUNDARY_POINT, m, n, RHO, U, V, P, CC, X, Y, gamma, cpu_time, argv[4], atof(argv[6]), argv[2]);
				}
			else        
				second_order_solver(&STEP, config, NUM_CELL, NUM_POINT, NUM_BOUNDARY, CELL_POINT, BOUNDARY_POINT, m, n, RHO, U, V, P, X, Y, gamma, cpu_time, argv[4], atof(argv[6]));
		}
	else       
		{
			printf("first order\n");
			if(strcmp(argv[argc-1],"Two_species")==0)
				{
					printf("Two species\n");		  
					first_order_two_species_solver(&STEP, config, NUM_CELL, NUM_POINT, NUM_BOUNDARY, CELL_POINT, BOUNDARY_POINT, m, n, RHO, U, V, P, CC, X, Y, gamma, cpu_time, argv[4], atof(argv[6]));
				}
			else
				first_order_solver(&STEP, config, NUM_CELL, NUM_POINT, NUM_BOUNDARY, CELL_POINT, BOUNDARY_POINT, m, n, RHO, U, V, P, X, Y, gamma, cpu_time, argv[4], atof(argv[6]), argv[2]);
		}
	
	//write the final data down.

	char address[50];
	if(strcmp(argv[argc-1],"second_order")==0||strcmp(argv[argc-2],"second_order")==0)
		strcpy(address,"2D_EUL_second_order\0" );
	else
		strcpy(address,"2D_EUL_first_order\0" );	

	if(strcmp(argv[argc-1],"Two_species")==0)
		file_two_species_write_TEC(NUM_POINT, X,x Y, NUM_CELL, CELL_POINT, RHO_t, U_t, V_t, P_t, CC_t, cpu_time, config, argv[2], address);
	else
		{		  
			//file_write_VTK(NUM_POINT, X, Y, NUM_CELL, CELL_POINT, RHO_t, U_t, V_t, P_t, cpu_time, config, argv[2], address); 
			file_write_TEC(NUM_POINT, X, Y, NUM_CELL, CELL_POINT, RHO_t, U_t, V_t, P_t, cpu_time, config, argv[2], address);  
		}
	*/

	

	
	Euler_scheme(&FV, &mv, scheme);


	file_write_TEC(FV, mv, argv[2], 0.0, (int)config[0]);

	if((int)config[0] > 1)
		file_write_VTK_3D(FV, mv, argv[2]);
	
/*
	for (int i = 0; i < 20; i++)
			printf("%d,%lf\n",i, config[i]);

	for (int i = 0; i < 1000; i++)
			printf("%d,%lf\n",i, FV.P[i]);
*/	
	return 0;	
}
