#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"


static int tri_mesh(struct mesh_var * mv)
{
	int k;
	const int num_cell = (int)config[3];

	if(isinf(config[13]) || isinf(config[14]))
		{
			fprintf(stderr, "The initial data is not mentioned in a structural mesh!\n");
			exit(2);
		}
	if(isinf(config[13]) || isinf(config[14]) || config[13] < config[4] || config[14] < config[4])
		{
			fprintf(stderr, "Without a proper spatial grid size!\n");
			exit(2);
		}
	
	const int n_x = (int)config[13], n_y = (int)config[14];

	mv->num_pt = (n_x+1)*(n_y+1);
	mv->X = malloc(mv->num_pt * sizeof(double));
	mv->Y = malloc(mv->num_pt * sizeof(double));
	if(mv->X == NULL || mv->Y == NULL)
		{
			printf("Not enough memory in square mesh constructed!\n");
			goto return_0;
		}
	for(k = 0; k < mv->num_pt; ++k)
		{
			mv->X[k] = (k%(n_x+1))*config[10];
			mv->Y[k] = (k/(n_x+1))*config[11];
		}	
	
	mv->cell_pt = malloc(num_cell * sizeof(void *));
	if(mv->cell_pt == NULL)
		{
			fprintf(stderr, "Not enough memory in square mesh constructed!\n");
			goto return_0;
		}
	for(k = 0; k < num_cell; ++k)
		{
			mv->cell_pt[k] = malloc(5 * sizeof(int));
			if(mv->cell_pt[k] == NULL)
				{
					for(int i = 0; i < k; ++i)
						{
							free(mv->cell_pt[i]);
							mv->cell_pt[i] = NULL;
						}
					fprintf(stderr, "Not enough memory in CELL_POINT[%d]!\n", k);
					goto return_0;
				}
			
			mv->cell_pt[k][0] = 4;
			mv->cell_pt[k][1] = k + k/n_x;
			mv->cell_pt[k][2] = mv->cell_pt[k][1] + 1;
			mv->cell_pt[k][3] = k + k/n_x + n_x + 2;
			mv->cell_pt[k][4] = mv->cell_pt[k][3] - 1;
		}

	mv->num_border[0] = 1;	
	const int num_border = 2*n_x + 2*n_y;
	mv->num_border[1] = num_border;	
	mv->border_pt = malloc((num_border+1) * sizeof(int));	
	mv->border_cond = malloc(num_border * sizeof(int));
	if(mv->border_pt == NULL || mv->border_cond == NULL)
		{
			printf("Not enough memory in square mesh constructed!\n");
			goto return_0;
		}
	
	for(k = 0; k < n_x+1; ++k)	
		mv->border_pt[k] = k; 
	for(k = n_x+1; k < n_x+n_y+1; ++k)	
		mv->border_pt[k] = mv->border_pt[k-1] + n_x + 1;
	for(k = n_x+n_y+1; k < n_x*2 + n_y + 1; ++k)	
		mv->border_pt[k] = mv->border_pt[k-1] - 1;
	for(k = n_x*2 + n_y + 1; k < num_border; ++k)	
		mv->border_pt[k] = mv->border_pt[k-1] - n_x - 1;
	mv->border_pt[num_border] = 0;

	return 1;

 return_0:
	free(mv->X);
	mv->X = NULL;
	free(mv->Y);
	mv->Y = NULL;	
	free(mv->border_pt);
	mv->border_pt = NULL;
	free(mv->border_cond);
	mv->border_cond = NULL;	
	if (mv->cell_pt != NULL)
		{
			for(int i = 0; i < num_cell; ++i)
				{
					free(mv->cell_pt[i]);
					mv->cell_pt[i] = NULL;
				}
			free(mv->cell_pt);
			mv->cell_pt = NULL;
		}
	return 0;	
}


int Tria_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * gamma, double * config, int m, int n)
{
	int num_boundary;
	int i, k;
	num_boundary = 2*m+2*n;
	BOUNDARY_POINT[0] = (int *)malloc(num_boundary * sizeof(int));
	BOUNDARY_POINT[1] = (int *)malloc(num_boundary * sizeof(int));
	if((BOUNDARY_POINT[0] == NULL)|(BOUNDARY_POINT[1] == NULL))
				{
					printf("NOT enough memory! BOUNDARY_POINT\n");
					exit(5);
				}

	for(k = 0; k < m*n; ++k)
		{
			CELL_POINT[k] = (int *)malloc(5 * sizeof(int));
			if(CELL_POINT[k] == NULL)
				{
					for(i = 0; i < k; ++i)
						{
							free(CELL_POINT[i]);
							CELL_POINT[i] = NULL;
						}
					printf("NOT enough memory! CELL_POINT[%d]\n", k);
					exit(5);
				}
		}
	//  CELL_POINT
	for(k = 0; k < m*n; ++k)
		{
		if(k%2)
			{
				CELL_POINT[k][0] = 3;
				CELL_POINT[k][1] = ((k-1)/2) + ((k-1)/2)/(n/2);
				CELL_POINT[k][2] = CELL_POINT[k][1] + 1;
				CELL_POINT[k][3] = CELL_POINT[k][1] + n/2 + 2;
			}
		else
			{
				CELL_POINT[k][0] = 3;
				CELL_POINT[k][1] = (k/2) + (k/2)/(n/2);
				CELL_POINT[k][2] = CELL_POINT[k][1] + n/2 + 2;
				CELL_POINT[k][3] = CELL_POINT[k][1] + n/2 + 1;
			}
		}
    // X, Y
	for(k = 0; k < (m+1)*(n/2+1); ++k)
	{
		X[k] = (k%(n/2+1))*config[2];
		Y[k] = (k/(n/2+1))*config[3];
	}
	// BOUNDARY
		for(k = 0; k < n/2+1; ++k)	
		{
			BOUNDARY_POINT[0][k] = k; 
		}
	for(k = n/2+1; k < n/2+m+1; ++k)	
		{
			BOUNDARY_POINT[0][k] = BOUNDARY_POINT[0][k-1] + n/2 + 1;
		}
	for(k = n/2+m+1; k < n + m + 1; ++k)	
		{
			BOUNDARY_POINT[0][k] = BOUNDARY_POINT[0][k-1] - 1;
		}
	for(k = n + m + 1; k < num_boundary; ++k)	
		{
			BOUNDARY_POINT[0][k] = BOUNDARY_POINT[0][k-1] - n/2 - 1;
		}


		for(k = 0; k < n; ++k)	
		{
			BOUNDARY_POINT[1][k] = -6; 
		}
	for(k = n; k < n+m; ++k)	
		{
			BOUNDARY_POINT[1][k] = -3; //prescribed boundary condition.
		}
	for(k = n+m; k < n*2 + m; ++k)	
		{
			BOUNDARY_POINT[1][k] = -6;
		}
	for(k = n*2 + m; k < num_boundary; ++k)	
		{
			BOUNDARY_POINT[1][k] = -3;
		}


		//gamma
		for(k = 0; k < m*n; ++k)
			{
                                gamma[k] = config[0];
			}
		
		return 	num_boundary;
}
