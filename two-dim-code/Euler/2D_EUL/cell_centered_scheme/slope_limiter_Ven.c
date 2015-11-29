#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "../../../../lib/custom.h"

void slope_limiter_Ven
(double * X_c, double * Y_c, double  * X, double * Y,
 double * grad_W_x, double * grad_W_y, 
 double *W[],  int NUM_CELL, double * config,
 int * CELL_CELL[], int * CELL_POINT[],
 int i, int m, int n)
{
	double const eps = config[4];

	int k,j;
	double X_c_n, Y_c_n;
	double M_c[2][2], M_c_i[2][2];
	double temp_grad_W_x, temp_grad_W_y;
	double W_c_min, W_c_max, W_c_x_p;
	double FAI_W;

	double const h_x = config[2];      // the length of the spatial grids
	double const h_y = config[3];      // the length of the spatial grids
	
	int CELL_RIGHT, STEP_RIGHT;     //For boundary condition
	int judge_x;

	

	for(k = 0; k < NUM_CELL; ++k)
		{
			M_c[0][0] = 0.0;
			M_c[0][1] = 0.0;
			M_c[1][0] = 0.0;
			M_c[1][1] = 0.0;
			grad_W_x[k] = 0.0;
			grad_W_y[k] = 0.0;

			judge_x = 1;
			for(j = 0; j < CELL_POINT[k][0]; ++j)
				{							
					if (CELL_CELL[k][j]>=0)
						{											
							X_c_n = X_c[CELL_CELL[k][j]];	
							Y_c_n = Y_c[CELL_CELL[k][j]];
							STEP_RIGHT = i;						
							CELL_RIGHT = CELL_CELL[k][j];
						}
					else if (CELL_CELL[k][j]==-4)//periodic boundary condition.
						{
							STEP_RIGHT = i;
							if(!(k%n))
								{  										   								
									CELL_RIGHT = k+n-1;
									X_c_n = X_c[k]-h_x;	
									Y_c_n = Y_c[k];
								}
							else if(k%n==n-1)
								{																					
									CELL_RIGHT = k-n+1;
									X_c_n = X_c[k]+h_x;	
									Y_c_n = Y_c[k];
								}
							else
								{																									
									printf("Something wrong as we construct boundary condition foe reconstruction.\n");
									exit(8);
								}
						}
					else
						{
							if((!(k%n))&&judge_x)
								{  										   								
									X_c_n = X_c[k]-h_x;	
									Y_c_n = Y_c[k];
									judge_x = 0;
								}
							else if(k%n==n-1&&judge_x)
								{																					
									X_c_n = X_c[k]+h_x;	
									Y_c_n = Y_c[k];
									judge_x = 0;
								}
							else if(k<n)
								{  										   								
									X_c_n = X_c[k];	
									Y_c_n = Y_c[k]-h_y;
								}
							else if(k>=(m-1)*n)
								{																					
									X_c_n = X_c[k];	
									Y_c_n = Y_c[k]+h_y;
								}
							else
								{																									
									printf("Something wrong as we construct boundary condition for reconstruction.\n");
									exit(8);
								}
									
							if (CELL_CELL[k][j]==-1)//initial boundary condition.
								{
									STEP_RIGHT = 0;
									CELL_RIGHT = k;
								}
							else if (CELL_CELL[k][j]==-3)//prescribed boundary condition.
								{
									STEP_RIGHT = i;
									CELL_RIGHT = k;
								}
							else if(CELL_CELL[k][j]==-2)//reflecting boundary condition.
								break;
							else
								{
									printf("No suitable boundary!\n");
									exit(7);
								}
						}
							
					M_c[0][0] = M_c[0][0] + (X_c_n - X_c[k]) * (X_c_n - X_c[k]);
					M_c[0][1] = M_c[0][1] + (X_c_n - X_c[k]) * (Y_c_n - Y_c[k]);
					M_c[1][0] = M_c[1][0] + (Y_c_n - Y_c[k]) * (X_c_n - X_c[k]);
					M_c[1][1] = M_c[1][1] + (Y_c_n - Y_c[k]) * (Y_c_n - Y_c[k]);
					grad_W_x[k] = grad_W_x[k] +  (W[STEP_RIGHT][CELL_RIGHT] - W[i][k]) * (X_c_n - X_c[k]);
					grad_W_y[k] = grad_W_y[k] +  (W[STEP_RIGHT][CELL_RIGHT] - W[i][k]) * (Y_c_n - Y_c[k]);
				}
					
			M_c_i[0][0] =   M_c[1][1]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
			M_c_i[0][1] = - M_c[0][1]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
			M_c_i[1][0] = - M_c[1][0]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
			M_c_i[1][1] =   M_c[0][0]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
			temp_grad_W_x = M_c_i[0][0] * grad_W_x[k] + M_c_i[0][1] * grad_W_y[k];
			temp_grad_W_y = M_c_i[1][0] * grad_W_x[k] + M_c_i[1][1] * grad_W_y[k];
			grad_W_x[k] = temp_grad_W_x;
			grad_W_y[k] = temp_grad_W_y;						
		}

	for(k = 0; k < NUM_CELL; ++k)
		{
			W_c_min = W[i][k];
			W_c_max = W[i][k];
			for(j = 0; j < CELL_POINT[k][0]; ++j)
				{
					if (CELL_CELL[k][j]>=0)
						{											
							STEP_RIGHT = i;							
							CELL_RIGHT = CELL_CELL[k][j];
						}
					else if (CELL_CELL[k][j]==-1)//initial boundary condition.
						{
							STEP_RIGHT = 0;
							CELL_RIGHT = k;
						}
					else if (CELL_CELL[k][j]==-3)//prescribed boundary condition.
						{
							STEP_RIGHT = i;
							CELL_RIGHT = k;
						}
					else if (CELL_CELL[k][j]==-4)//periodic boundary condition in x-direction.
						{
							STEP_RIGHT = i;
							if(!(k%n))
								CELL_RIGHT = k+n-1;
							else if(k%n==n-1)
								CELL_RIGHT = k-n+1;
							else
								{																									
									printf("Something wrong as we construct periodic boundary condition in x-direction.\n");
									exit(8);
								}
						}
					else if (CELL_CELL[k][j]==-2)//reflecting boundary condition.
						break;
					else
						{
							printf("No suitable boundary!\n");
							exit(7);
						}	
							
					if(W[STEP_RIGHT][CELL_RIGHT] < W_c_min)
						W_c_min = W[STEP_RIGHT][CELL_RIGHT];
					else if(W[STEP_RIGHT][CELL_RIGHT] > W_c_max)
						W_c_max = W[STEP_RIGHT][CELL_RIGHT];
				}
			FAI_W = DBL_MAX;
			for(j = 0; j < CELL_POINT[k][0]; ++j)
				{
					for(j = 0; j < CELL_POINT[k][0]; ++j)
						{
							if (CELL_CELL[k][j]>=0)
								{											
									STEP_RIGHT = i;							
									CELL_RIGHT = CELL_CELL[k][j];
								}
							else if (CELL_CELL[k][j]==-1)//initial boundary condition.
								{
									STEP_RIGHT = 0;
									CELL_RIGHT = k;
								}
							else if (CELL_CELL[k][j]==-3)//prescribed boundary condition.
								{
									STEP_RIGHT = i;
									CELL_RIGHT = k;
								}
							else if (CELL_CELL[k][j]==-4)//periodic boundary condition in x-direction.
								{
									STEP_RIGHT = i;
									if(!(k%n))
										CELL_RIGHT = k+n-1;
									else if(k%n==n-1)
										CELL_RIGHT = k-n+1;
									else
										{																									
											printf("Something wrong as we construct periodic boundary condition in x-direction.\n");
											exit(8);
										}
								}
							else if (CELL_CELL[k][j]==-2)//reflecting boundary condition.
								break;
							else
								{
									printf("No suitable boundary!\n");
									exit(7);
								}
							
							W_c_x_p = W[i][k] + grad_W_x[k] * (X[CELL_POINT[k][j+1]] - X_c[k]) + grad_W_y[k] * (Y[CELL_POINT[k][j+1]] - Y_c[k]);
							if (fabs(W_c_x_p - W[i][k]) <eps)
								FAI_W = (FAI_W < 1.0) ? FAI_W : 1.0;	
							else if((W_c_x_p - W[i][k]) > 0.0)
								FAI_W = (FAI_W < miu_Ven((W_c_max - W[i][k])/(W_c_x_p - W[i][k]))) ? FAI_W : miu_Ven((W_c_max - W[i][k])/(W_c_x_p - W[i][k]));
							else
								FAI_W = (FAI_W < miu_Ven((W_c_min - W[i][k])/(W_c_x_p - W[i][k]))) ? FAI_W : miu_Ven((W_c_min - W[i][k])/(W_c_x_p - W[i][k]));
						}
					grad_W_x[k] = grad_W_x[k] * FAI_W;
					grad_W_y[k] = grad_W_y[k] * FAI_W;
				}	   	
		}
}
