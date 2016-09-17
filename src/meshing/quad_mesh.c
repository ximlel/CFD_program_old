#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "../include/var_struc.h"
#include "../include/meshing.h"

#ifndef M_PI
#define M_PI acos(-1.0)
#endif



static int quad_mesh(struct mesh_var * mv, int n_x, int n_y)
{
	int k;
	const int num_cell = n_x * n_y;
	if (num_cell > (int)config[3])
		printf("There are %d ghost cell!\n", mv->num_ghost = num_cell - (int)config[3]);
	else if (num_cell < (int)config[3])
		{
			fprintf(stderr,"There are not enough cell in .msh file!\n");
			exit(2);
		}
	
	if(isinf(config[13]) || isinf(config[14]))
		{
			fprintf(stderr, "The initial data is not mentioned in a structural mesh!\n");
			exit(2);
		}
	if(isinf(config[10]) || isinf(config[11]) || config[10] < 0.0 || config[11] < 0.0)
		{
			fprintf(stderr, "Without a proper spatial grid size!\n");
			exit(2);
		}

	mv->num_pt = (n_x+1)*(n_y+1);
	mv->X = malloc(mv->num_pt * sizeof(double));
	mv->Y = malloc(mv->num_pt * sizeof(double));
	if(mv->X == NULL || mv->Y == NULL)
		{
			printf("Not enough memory in quadrilateral mesh constructed!\n");
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
			fprintf(stderr, "Not enough memory in quadrilateral mesh constructed!\n");
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
			printf("Not enough memory in quadrilateral mesh constructed!\n");
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

static void quad_border_cond(struct mesh_var * mv, int n_x, int n_y, int down, int right, int up, int left)
{
	if (down >= 0 || right >= 0 || up >= 0|| left >= 0)
		{
			fprintf(stderr, "Input wrong boundary condition in quadrilateral mesh!\n");
			exit(2);
		}
	else if ((up == -7) - (down == -7) != 0 || (left == -7) - (right == -7) != 0)
		{
			fprintf(stderr, "Periodic boundary condition error!\n");
			exit(2);
		}

	int k;
	const int num_cell = n_x * n_y;
	const int num_border = mv->num_border[1];
	
	mv->peri_cell = malloc(num_cell * sizeof(int));
	if(mv->peri_cell == NULL)
		{
			printf("Not enough memory in quad periodic boundary constructed!\n");
			exit(5);
		}	
	for(k = 0; k < num_cell; k++)
		mv->peri_cell[k] = -1;
	
	for(k = 0; k < n_x; ++k)
		{
			if (down == -7)
				mv->peri_cell[k] = k + n_x * (n_y-2);				
			// mv->border_cond[k] = k + n_x * (n_y-1);
			mv->border_cond[k] = down;
		}
	for(k = n_x; k < n_x+n_y; ++k)
		{
			if (right == -7)				
				mv->peri_cell[n_x * (k - n_x + 1) - 1] = n_x * (k - n_x) + 1;			
			// mv->border_cond[k] = n_x * (k - n_x);
			mv->border_cond[k] = right;
		}
	for(k = n_x+n_y; k < n_x*2 + n_y; ++k)
		{
			if (up == -7)
				mv->peri_cell[n_x*(n_y + 1) + n_y - k - 1] = n_x*3 + n_y - k - 1;					
			// mv->border_cond[k] = n_x*2 + n_y - k - 1;
			mv->border_cond[k] = up;
		}
	for(k = n_x*2 + n_y; k < num_border; ++k)
		{
			if (left == -7)
				mv->peri_cell[(num_border - k -1) * n_x] = (num_border - k) * n_x - 2;				
			// mv->border_cond[k] = (mv->num_border[1] - k) * n_x - 1;
			mv->border_cond[k] = left;
		}

	if(mv->num_ghost > 0)
		period_cell_modi(mv);
}

void Sod_mesh(struct mesh_var * mv)
{
	const int n_x = (int)config[13], n_y = (int)config[14];
	if (quad_mesh(mv, n_x, n_y) == 0)
		exit(5);

	quad_border_cond(mv, n_x, n_y, -2, -3, -2, -3);
}

void Shear_mesh(struct mesh_var * mv)
{
	const int n_x = (int)config[13], n_y = (int)config[14];		
	if (quad_mesh(mv, n_x, n_y) == 0)
		exit(5);

	quad_border_cond(mv, n_x, n_y, -1, -3, -3, -3);
}

void free_mesh(struct mesh_var * mv)
{
	const int n_x = (int)config[13], n_y = (int)config[14];				
	if (quad_mesh(mv, n_x, n_y) == 0)
		exit(5);

	quad_border_cond(mv, n_x, n_y, -3, -3, -3, -3);
}

void RMI_mesh(struct mesh_var * mv)
{
	const int n_x = (int)config[13] + 2, n_y = (int)config[14];				
	if (quad_mesh(mv, n_x, n_y) == 0)
		exit(5);

	quad_border_cond(mv, n_x, n_y, -3, -7, -3, -7);
}

void cylinder_mesh(struct mesh_var * mv)
{
	const int n_x = (int)config[13], n_y = (int)config[14];
	if (quad_mesh(mv, n_x, n_y) == 0)
		exit(5);

	const double R = config[11]*n_y/M_PI*9.0/8.0;
	
	for(int k = 0; k < (n_x+1)*(n_y+1); ++k)
		{
			mv->X[k] = (R+(n_x-k%(n_x+1))*config[10])*cos(13.0/9.0*M_PI-(k/(n_x+1))*8.0/9.0*M_PI/n_y);
			mv->Y[k] = (R+(n_x-k%(n_x+1))*config[10])*sin(13.0/9.0*M_PI-(k/(n_x+1))*8.0/9.0*M_PI/n_y);
		}
	
	quad_border_cond(mv, n_x, n_y, -3, -2, -3, -1);
}

void odd_even_mesh(struct mesh_var * mv)
{
	const int n_x = (int)config[13], n_y = (int)config[14];

	if (quad_mesh(mv, n_x, n_y) == 0)
		exit(5);

	for(int k = 0; k < (n_y+1); ++k)
		mv->Y[(n_x/2)*(n_y+1)+k] += ((k%2)-0.5)*0.002*config[11];

	quad_border_cond(mv, n_x, n_y, -2, -3, -2, -1);
}

void odd_even_periodic_mesh(struct mesh_var * mv)
{
	const int n_x = (int)config[13] + 2, n_y = (int)config[14];
	if (quad_mesh(mv, n_x, n_y) == 0)
		exit(5);

	for(int k = 0; k < (n_y+1); ++k)	
		mv->Y[(n_x/2)*(n_y+1)+k] += ((k%2)-0.5)*0.002*config[11];

	quad_border_cond(mv, n_x, n_y, -2, -7, -2, -7);
}

void odd_even_inflow_mesh(struct mesh_var * mv)
{
	const int n_x = (int)config[13], n_y = (int)config[14] + 2;
	if (quad_mesh(mv, n_x, n_y) == 0)
		exit(5);

	for(int k = 0; k < (n_y+1); ++k)
		mv->Y[(n_x/2)*(n_y+1)+k] += ((k%2)-0.5)*0.002*config[11];

	quad_border_cond(mv, n_x, n_y, -7, -3, -7, -1);
}

void rand_disturb_inflow_mesh(struct mesh_var * mv)
{
	const int n_x = (int)config[13], n_y = (int)config[14] + 2;
	if (quad_mesh(mv, n_x, n_y) == 0)
		exit(5);

	srand((unsigned) time(NULL)); //seed--time.

	for(int k = n_x+1; k < n_y*(n_x+1); ++k)
		{
			if(k%(n_x+1)&&(k%(n_x+1))!=n_x+1)
				{							
					mv->X[k] += (0.5-(rand()%10001)/10000.0)*0.000001*config[11];
					mv->Y[k] += (0.5-(rand()%10001)/10000.0)*0.000001*config[11];
				}
		}

	quad_border_cond(mv, n_x, n_y, -7, -3, -7, -1);
}
