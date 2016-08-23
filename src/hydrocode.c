/* 
 * This is a implementation of 1-D or 2-D hydrocodes in Cartesian geometry
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "file_io.h"
#include "meshing.h"
#include "finite_volume.h"
#include "Riemann_solver.h"
#include "tools.h"


/*The global primitive variable and configuration array.
 */

double config[200] = {1.0/0.0};

int DIM;



int main(int argc, char *argv[])
{
	/* Set dimension.
	 */
	DIM = atoi(argv[5]);
	
	/* Find the directory of input data.
	 */
	char add_mkdir[PATH_MAX];
	switch(DIM)
		{
		case 1 :
			strcpy(add_mkdir, "../data_in/one-dim/");
			break;
		case 2 :
			strcpy(add_mkdir, "../data_in/two-dim/");
			break;
		default :
			printf("Strange computational dimension.\n");
		}	
	DIR * dir_test = NULL;
	strcat(add_mkdir, argv[1]);

	dir_test = opendir(add_mkdir);
	if(dir_test == NULL)
		{
			printf("Input directory is not exist.\n");
			exit(1);
		}
	closedir(dir_test);


	/* Firstly we read the initial data file. 
	 */	
	printf("%s is configurated:\n", argv[1]);

	char addconfig[PATH_MAX];
	strcpy(addconfig, add_mkdir);
	strcat(addconfig, "/config.txt");
	configurate(addconfig);

	char addRHO[PATH_MAX];
	strcpy(addRHO, add_mkdir);
	strcat(addRHO, "/RHO.txt");
	initialize(addRHO, RHO);
	char addU[PATH_MAX];
	strcpy(addU, add_mkdir);
	strcat(addU, "/U.txt");
	initialize(addU, U);
	char addP[PATH_MAX];
	strcpy(addP, add_mkdir);
	strcat(addP, "/P.txt");
	initialize(addP, P);		
	if(atoi(argv[5])>1)
		{		
			char addV[PATH_MAX];
			strcpy(addV, add_mkdir);
			strcat(addV, "/V.txt");
			initialize(addV, V);
		}
	switch((int)config[2])
		{
		case 2 :
			{							
				char addPHI[PATH_MAX];
				strcpy(addPHI, add_mkdir);
				strcat(addPHI, "/PHI.txt");
				initialize(addPHI, PHI);
				break;
			}
		}

	printf("%s is initialized, grid number = %d.\n", argv[1], (int)config[3]);


	int NUM_POINT;
	NUM_POINT = (m+1)*(n+1);


	int i = 0, k = 0, N = (int)(config[5]);
  
	double cpu_time[N];
	cpu_time[0] = 0.0;


	

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

		

	free(RHO);
	free(U);
	free(P);
	free(V);
	free(PHI);
	RHO = NULL;
	U = NULL;
	P = NULL;
	V = NULL;
	PHI = NULL;



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
