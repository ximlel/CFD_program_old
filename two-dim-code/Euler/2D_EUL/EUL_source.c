/* 
 * This is a implementation of 2-D Roe's scheme in Cartesian geometry
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "../../lib/file_io.h"
#include "meshing.h"
#include "cell_centered_scheme.h"
#include "../../lib/Riemann_solver.h"
#include "../../lib/custom.h"



#ifndef N_CONF
#define N_CONF 6
#endif /* N_CONF */

double * U0 = NULL;
double * P0 = NULL;
double * V0 = NULL;
double * RHO0 = NULL;
double * CC0 = NULL;

double SSTTAARR[3];

int main(int argc, char *argv[])
{
	int stat_mkdir = 0, len;
	char add_mkdir[200] = "../../../data_in/two-dim/\0";
	DIR * dir_test = NULL;
	strcat(add_mkdir, argv[1]);
	strcat(add_mkdir, "\0");

	dir_test = opendir(add_mkdir);
	if(dir_test == NULL)
		{
			printf("Input directory \"%s\"\n is not exist.\n", add_mkdir);
			return 1;
		}
	closedir(dir_test);

	char addRHO[200];
	strcpy(addRHO, add_mkdir);
	strcat(addRHO, "/RHO.txt\0");
	char addU[200];
	strcpy(addU, add_mkdir);
	strcat(addU, "/U.txt\0");
	char addV[200];
	strcpy(addV, add_mkdir);
	strcat(addV, "/V.txt\0");
	char addP[200];
	strcpy(addP, add_mkdir);
	strcat(addP, "/P.txt\0");
	char addCC[200];
	strcpy(addCC, add_mkdir);
	strcat(addCC, "/CC.txt\0");
	char addconfig[200];
	strcpy(addconfig, add_mkdir);
	strcat(addconfig, "/config.txt\0");


	if(strcmp(argv[argc-1],"Two_species")==0)
		initialize_CC(argv[0], addCC);

	initialize(argv[0], addRHO, addU, addV, addP);  /* Firstly we read the initial
													 * data file. The function 
													 * initialize return a point
													 * pointing to the position
													 * of a block of memory
													 * consisting (m*n+2) variables
													 * of type double.
													 * The first two value of these
													 * variables is m(line) and n(column). The
													 * following m*n variables
													 * are the initial value.
													 */

	int m = (int)RHO0[0], n = (int)RHO0[1];  
	int NUM_CELL;
	NUM_CELL = m*n;			/* m*n is the number of initial value
							 * as well as the number of grids.
							 * As m*n is frequently use to
							 * represent the number of grids,
							 * we do not use the name such as
							 * num_grid here to correspond to
							 * notation in the math theory.
							 */

	int NUM_POINT;
	NUM_POINT = (m+1)*(n+1);

	int NUM_BOUNDARY;
	double config[N_CONF];  /* config[0] is the polytropic index
							 * config[1] is the total time
							 * config[2] is the spatial grid size in x direction
							 * config[3] is the spatial grid size in y direction
							 * config[4] is the largest value can be seen as zero
							 * config[5] is the maximum number of time steps
							 */
	configurate(config, argv[0], addconfig); /* Read the config-
											  * uration data.
											  * The detail could
											  * be seen in the
											  * definition of
											  * array config.
											  */
	int i = 0, k = 0, N = (int)(config[5]);
  
	double cpu_time[N];
	cpu_time[0] = 0.0;


	double RHO_t[NUM_CELL];
	for(k=0; k<NUM_CELL; k++)
		RHO_t[k] = (RHO0+2)[k]; 

	double U_t[NUM_CELL];
	for(k=0; k<NUM_CELL; k++)
		U_t[k] = (U0+2)[k];
	
	double V_t[NUM_CELL];
	for(k=0; k<NUM_CELL; k++)
		V_t[k] = (V0+2)[k];
		
  	double P_t[NUM_CELL];
	for(k=0; k<NUM_CELL; k++)
		P_t[k] = (P0+2)[k];
		
  	double CC_t[NUM_CELL];
	if(strcmp(argv[argc-1],"Two_species")==0)
		for(k=0; k<NUM_CELL; k++)
			CC_t[k] = (CC0+2)[k];

	
	double *RHO[2];
	double *U[2];
	double *V[2];
	double *P[2];
	double *CC[2];
	
	RHO[0] = RHO0+2;
	RHO[1] = RHO_t;
	U[0] = U0+2;
	U[1] = U_t;
	V[0] = V0+2;
	V[1] = V_t;
	P[0] = P0+2;
	P[1] = P_t;	
	CC[0] = CC0+2;
	CC[1] = CC_t;
	

	double X[NUM_POINT];
	double Y[NUM_POINT];
  
                                            
	int * CELL_POINT[NUM_CELL];
  
	int * BOUNDARY_POINT[2];

	double gamma[NUM_CELL];

	printf("%s\n", argv[2]);

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
			file_two_species_write_TEC(NUM_POINT, X, Y, NUM_CELL, CELL_POINT, RHO_t, U_t, V_t, P_t, CC_t, cpu_time, config, argv[2], address);
		else
			{		  
				//file_write_VTK(NUM_POINT, X, Y, NUM_CELL, CELL_POINT, RHO_t, U_t, V_t, P_t, cpu_time, config, argv[2], address); 
				file_write_TEC(NUM_POINT, X, Y, NUM_CELL, CELL_POINT, RHO_t, U_t, V_t, P_t, cpu_time, config, argv[2], address);  
			}

		

	free(RHO0);
	free(U0);
	free(V0);
	free(P0);
	free(CC0);
	RHO0 = NULL;
	U0 = NULL;
	V0 = NULL;
	P0 = NULL;
	CC0 = NULL;
	RHO[0] = NULL;
	RHO[1] = NULL;
	U[0] = NULL;
	U[1] = NULL;
	V[0] = NULL;
	V[1] = NULL;
	P[0] = NULL;
	P[1] = NULL;	
	CC[0] = NULL;
	CC[1] = NULL;



	for(k = 0; k < NUM_CELL; ++k)
		{
			free(CELL_POINT[k]);
			CELL_POINT[k] = NULL;
		}
	free(BOUNDARY_POINT[0]);  
	BOUNDARY_POINT[0] = NULL;
	free(BOUNDARY_POINT[1]);  
	BOUNDARY_POINT[1] = NULL;

	printf("\n");
	return 0;
}
