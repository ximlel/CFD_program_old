#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#define min(x,y)  ( x<y?x:y )

#include <time.h>


#include "../../../lib/Riemann_solver.h"
#include "../../../lib/custom.h"

/* This function use second order scheme to solve 2-D
 * equations of motion by Lagrangian method.
 *
 *config is the array of configuration data, the detail
 *         could be seen in the comments of the main function.
 *NUM_CELL is the number of the grids.
 */

//
int second_order_solver_GLACE
(double * config, int NUM_CELL, int NUM_POINT, int NUM_BOUNDARY, int * CELL_POINT[], int * BOUNDARY_POINT[],
 double * RHO[], double * U[], double * V[], double * P[], double * X[], double * Y[],
 double * NORMAL_VELOCITY[], double * gamma, double * cpu_time, char * limiter, double CFL)
{


	double (*miu)(double);
	
	if(strcmp(limiter,"BJ")==0)	
		miu = miu_BJ;
	else if(strcmp(limiter,"Ven")==0)
		miu = miu_Ven;
        else
        printf("No limiter!\n");
	
	int i, j, k, l;
	int NUM_CC;
	double S, S_tri;

  clock_t tic, toc;
  double sum = 0.0;

  int const N = (int)(config[5]);  // the number of time steps
  double const eps = config[4];    // the largest value could be
                                   // seen as zero
  double tau;
  tau = DBL_MAX;
  double t_all;
  t_all = 0;

  
  double lambda[NUM_CELL]; //minimal distance
  
  double dt_volume_temp;
  
  int stop_step=0;


  int p_p,p_n;

  int * CELL_CELL[NUM_CELL];
  int CELL_record[NUM_CELL];
  int CELL_save[NUM_CELL];

  	for(k = 0; k < NUM_CELL; ++k)
		{ 
  			NUM_CC=0;
			
			for(i = 0; i < NUM_CELL; ++i)				
					CELL_record[i] = 0;				

			for(j = 0; j < CELL_POINT[k][0]; ++j)
				{					
					for(i = 0; i < NUM_CELL; ++i)
						{
							if(i!=k)
								{								
									for(l = 0; l < CELL_POINT[i][0]; ++l)
										{
											if((CELL_POINT[i][l+1]==CELL_POINT[k][j+1])&&(CELL_record[i]==0))
												{
													CELL_record[i] = 1;
													CELL_save[NUM_CC] = i;
													NUM_CC = NUM_CC + 1;
												}
										}						  
								}
						}
					
				}

			CELL_CELL[k] = (int *)malloc((NUM_CC+1) * sizeof(int));
			if(CELL_CELL[k] == NULL)
				{
					for(i = 0; i < k; ++i)
						{
							free(CELL_CELL[i]);
							CELL_CELL[i] = NULL;
						}
					printf("NOT enough memory! CELL_CELL[%d]\n", k);
					exit(5);
				}
                        CELL_CELL[k][0] = NUM_CC;
						
						for(i = 0; i < NUM_CC; ++i)
							CELL_CELL[k][i+1] = CELL_save[i];						
		}

	double M_c[2][2], M_c_i[2][2];
	
	double grad_P_x[NUM_CELL], grad_P_y[NUM_CELL];
	double temp_grad_P_x, temp_grad_P_y;	
	double P_c_min, P_c_max, P_c_x_p;
	double FAI_P[NUM_CELL];

	double grad_U_x[NUM_CELL], grad_U_y[NUM_CELL];
	double temU_grad_U_x, temU_grad_U_y;	
	double U_c_min, U_c_max, U_c_x_U;
	double FAI_U[NUM_CELL];

	double grad_V_x[NUM_CELL], grad_V_y[NUM_CELL];
	double temV_grad_V_x, temV_grad_V_y;	
	double V_c_min, V_c_max, V_c_x_V;
	double FAI_V[NUM_CELL];

	
  double MASS[NUM_CELL];
  double volume_temp;
  for(k = 0; k < NUM_CELL; ++k)
	  {
		  volume_temp = 0;
		  for(j = 0; j < CELL_POINT[k][0]; ++j)
			  {
				  if(j == CELL_POINT[k][0]-1) 
					  {
						  p_p=CELL_POINT[k][1];
						  p_n=CELL_POINT[k][j+1];
					  }				  
				  else
					  {
						  p_p=CELL_POINT[k][j+2];
						  p_n=CELL_POINT[k][j+1];
					  } 
				  volume_temp = volume_temp + 0.5 * (X[0][p_n] * Y[0][p_p] - Y[0][p_n] * X[0][p_p]);
			  }
		  MASS[k] = volume_temp * RHO[0][k];
	  }

  double X_c[NUM_CELL], Y_c[NUM_CELL];

  double * l_pc[NUM_CELL];
  initialize_memory(l_pc,NUM_CELL,CELL_POINT);
  double * n_pc_x[NUM_CELL];
  initialize_memory(n_pc_x,NUM_CELL,CELL_POINT);
  double * n_pc_y[NUM_CELL];
  initialize_memory(n_pc_y,NUM_CELL,CELL_POINT);
  double * dt_l_n_pc_x[NUM_CELL];
  initialize_memory(dt_l_n_pc_x,NUM_CELL,CELL_POINT);
  double * dt_l_n_pc_y[NUM_CELL];
  initialize_memory(dt_l_n_pc_y,NUM_CELL,CELL_POINT);
  double * dt_n_pc_x[NUM_CELL];
  initialize_memory(dt_n_pc_x,NUM_CELL,CELL_POINT);
  double * dt_n_pc_y[NUM_CELL];
  initialize_memory(dt_n_pc_y,NUM_CELL,CELL_POINT);

  double * posi_l_pc[NUM_CELL];
  initialize_memory(posi_l_pc,NUM_CELL,CELL_POINT);
  double * posi_n_pc_x[NUM_CELL];
  initialize_memory(posi_n_pc_x,NUM_CELL,CELL_POINT);
  double * posi_n_pc_y[NUM_CELL];
  initialize_memory(posi_n_pc_y,NUM_CELL,CELL_POINT);
  double * dt_posi_l_n_pc_x[NUM_CELL];
  initialize_memory(dt_posi_l_n_pc_x,NUM_CELL,CELL_POINT);
  double * dt_posi_l_n_pc_y[NUM_CELL];
  initialize_memory(dt_posi_l_n_pc_y,NUM_CELL,CELL_POINT);
  double * dt_posi_n_pc_x[NUM_CELL];
  initialize_memory(dt_posi_n_pc_x,NUM_CELL,CELL_POINT);
  double * dt_posi_n_pc_y[NUM_CELL];
  initialize_memory(dt_posi_n_pc_y,NUM_CELL,CELL_POINT);

  double * nega_l_pc[NUM_CELL];
  initialize_memory(nega_l_pc,NUM_CELL,CELL_POINT);
  double * nega_n_pc_x[NUM_CELL];
  initialize_memory(nega_n_pc_x,NUM_CELL,CELL_POINT);
  double * nega_n_pc_y[NUM_CELL];
  initialize_memory(nega_n_pc_y,NUM_CELL,CELL_POINT);
  double * dt_nega_l_n_pc_x[NUM_CELL];
  initialize_memory(dt_nega_l_n_pc_x,NUM_CELL,CELL_POINT);
  double * dt_nega_l_n_pc_y[NUM_CELL];
  initialize_memory(dt_nega_l_n_pc_y,NUM_CELL,CELL_POINT);
  double * dt_nega_n_pc_x[NUM_CELL];
  initialize_memory(dt_nega_n_pc_x,NUM_CELL,CELL_POINT);
  double * dt_nega_n_pc_y[NUM_CELL];
  initialize_memory(dt_nega_n_pc_y,NUM_CELL,CELL_POINT);


  double M_pc[2][2];

  double * M_pc_00[NUM_CELL];
  initialize_memory(M_pc_00,NUM_CELL,CELL_POINT);
  double * M_pc_01[NUM_CELL];
  initialize_memory(M_pc_01,NUM_CELL,CELL_POINT);
  double * M_pc_10[NUM_CELL];
  initialize_memory(M_pc_10,NUM_CELL,CELL_POINT);
  double * M_pc_11[NUM_CELL];
  initialize_memory(M_pc_11,NUM_CELL,CELL_POINT);
  double * dt_M_pc_00[NUM_CELL];
  initialize_memory(dt_M_pc_00,NUM_CELL,CELL_POINT);
  double * dt_M_pc_01[NUM_CELL];
  initialize_memory(dt_M_pc_01,NUM_CELL,CELL_POINT);
  double * dt_M_pc_10[NUM_CELL];
  initialize_memory(dt_M_pc_10,NUM_CELL,CELL_POINT);
  double * dt_M_pc_11[NUM_CELL];
  initialize_memory(dt_M_pc_11,NUM_CELL,CELL_POINT);

  double U_p_x[NUM_POINT];
  double U_p_y[NUM_POINT];
  double PAI_p[NUM_POINT];
  double dt_U_p_x[NUM_POINT];
  double dt_U_p_y[NUM_POINT];

  
  double temp_U_p_x, temp_U_p_y;

  int index[NUM_POINT];
  	for(l = 0; l < NUM_POINT; ++l)
		{
			index[l] = -1;
			for(k = 0; k < NUM_BOUNDARY; ++k)
				{
					if(BOUNDARY_POINT[0][k] == l) index[l] = k;
				}
		}

	double M_p[2][2], M_p_i[2][2];

	double A_p[3][3], A_p_i[3][3];
	double A_p_det;

	double b_temp_x, b_temp_y, b_temp_z;

	double posi_l_p[NUM_BOUNDARY], posi_n_p_x[NUM_BOUNDARY], posi_n_p_y[NUM_BOUNDARY];
	double nega_l_p[NUM_BOUNDARY], nega_n_p_x[NUM_BOUNDARY], nega_n_p_y[NUM_BOUNDARY];
	double dt_posi_n_p_x[NUM_BOUNDARY], dt_posi_n_p_y[NUM_BOUNDARY];
	double dt_nega_n_p_x[NUM_BOUNDARY], dt_nega_n_p_y[NUM_BOUNDARY];

	double * F_pc_x[NUM_CELL];
	initialize_memory(F_pc_x,NUM_CELL,CELL_POINT);
	double * F_pc_y[NUM_CELL];
	initialize_memory(F_pc_y,NUM_CELL,CELL_POINT);
	double * dt_F_pc_x[NUM_CELL];
	initialize_memory(dt_F_pc_x,NUM_CELL,CELL_POINT);
	double * dt_F_pc_y[NUM_CELL];
	initialize_memory(dt_F_pc_y,NUM_CELL,CELL_POINT);


	double * P_c_lim[NUM_CELL];
	initialize_memory(P_c_lim,NUM_CELL,CELL_POINT);
	double * U_c_lim[NUM_CELL];
	initialize_memory(U_c_lim,NUM_CELL,CELL_POINT);
	double * V_c_lim[NUM_CELL];
	initialize_memory(V_c_lim,NUM_CELL,CELL_POINT);

	double a[NUM_CELL];   //acoustic

	double * POINT_POINT[NUM_POINT];
	double POINT_save[NUM_CELL*2];
	int NUM_PP;
	for(l = 0; l < NUM_POINT; ++l)
		{
		NUM_PP=0;
			for(k = 0; k < NUM_CELL; ++k)
					{
						for(j = 0; j < CELL_POINT[k][0]; ++j)
							{							   
								if (CELL_POINT[k][j+1]==l)
									{
										POINT_save[NUM_PP*2]=k;
										POINT_save[NUM_PP*2+1]=j;
										NUM_PP++;
									}
							}
					}
		POINT_POINT[l] = (double *)malloc((NUM_PP*2+1) * sizeof(double));
			if(POINT_POINT[l] == NULL)
				{
					for(i = 0; i < l; ++i)
						{
							free(POINT_POINT[i]);
							POINT_POINT[i] = NULL;
						}
					printf("NOT enough memory! POINT_POINT[%d]\n", l);
					exit(5);
				}
                        POINT_POINT[l][0] = NUM_PP;
						
						for(i = 0; i < NUM_PP*2; ++i)
							POINT_POINT[l][i+1] = POINT_save[i];

		}


	printf("Grid has been constructed.\n");

  
//------------THE MAIN LOOP-------------
	for(i = 0; i < N; ++i)
		{	

			tic = clock();

			for(k = 0; k < NUM_CELL; ++k)
				{
					a[k] = sqrt(gamma[k] * P[i][k] / RHO[i][k]);
				}

			for(k = 0; k < NUM_CELL; ++k)
				{
					S = 0;
					X_c[k] = 0;
					Y_c[k] = 0;
					for(j = 2; j < CELL_POINT[k][0]; ++j)
						{
							S_tri =X[i][CELL_POINT[k][1]]*Y[i][CELL_POINT[k][j]] + X[i][CELL_POINT[k][j+1]]*Y[i][CELL_POINT[k][1]] + X[i][CELL_POINT[k][j]]*Y[i][CELL_POINT[k][j+1]] - X[i][CELL_POINT[k][j+1]]*Y[i][CELL_POINT[k][j]] - X[i][CELL_POINT[k][1]]*Y[i][CELL_POINT[k][j+1]] - X[i][CELL_POINT[k][j]]*Y[i][CELL_POINT[k][1]];
							X_c[k] = X_c[k] + (X[i][CELL_POINT[k][1]] + X[i][CELL_POINT[k][j]] + X[i][CELL_POINT[k][j+1]]) * S_tri;
							Y_c[k] = Y_c[k] + (Y[i][CELL_POINT[k][1]] + Y[i][CELL_POINT[k][j]] + Y[i][CELL_POINT[k][j+1]]) * S_tri;
							S = S + S_tri;
						}			
					X_c[k] = X_c[k]/S/3;
					Y_c[k] = Y_c[k]/S/3;
				}



			for(k = 0; k < NUM_CELL; ++k)
				{
					M_c[0][0] = 0;
					M_c[0][1] = 0;
					M_c[1][0] = 0;
					M_c[1][1] = 0;
					grad_P_x[k] = 0;
					grad_P_y[k] = 0;
					
					for(j = 0; j < CELL_CELL[k][0]; ++j)
						{
							M_c[0][0] = M_c[0][0] + (X_c[CELL_CELL[k][j+1]] - X_c[k]) * (X_c[CELL_CELL[k][j+1]] - X_c[k]);
							M_c[0][1] = M_c[0][1] + (X_c[CELL_CELL[k][j+1]] - X_c[k]) * (Y_c[CELL_CELL[k][j+1]] - Y_c[k]);
							M_c[1][0] = M_c[1][0] + (Y_c[CELL_CELL[k][j+1]] - Y_c[k]) * (X_c[CELL_CELL[k][j+1]] - X_c[k]);
							M_c[1][1] = M_c[1][1] + (Y_c[CELL_CELL[k][j+1]] - Y_c[k]) * (Y_c[CELL_CELL[k][j+1]] - Y_c[k]);
							grad_P_x[k] = grad_P_x[k] +  (P[i][CELL_CELL[k][j+1]] - P[i][k]) * (X_c[CELL_CELL[k][j+1]] - X_c[k]);
							grad_P_y[k] = grad_P_y[k] +  (P[i][CELL_CELL[k][j+1]] - P[i][k]) * (Y_c[CELL_CELL[k][j+1]] - Y_c[k]);
						}
					
					M_c_i[0][0] =   M_c[1][1]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
					M_c_i[0][1] = - M_c[0][1]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
					M_c_i[1][0] = - M_c[1][0]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
					M_c_i[1][1] =   M_c[0][0]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
					temp_grad_P_x = M_c_i[0][0] * grad_P_x[k] + M_c_i[0][1] * grad_P_y[k];
					temp_grad_P_y = M_c_i[1][0] * grad_P_x[k] + M_c_i[1][1] * grad_P_y[k];
					grad_P_x[k] = temp_grad_P_x;
					grad_P_y[k] = temp_grad_P_y;						
				}

			for(k = 0; k < NUM_CELL; ++k)
				{
					P_c_min = P[i][k];
					P_c_max = P[i][k];
					for(j = 0; j < CELL_CELL[k][0]; ++j)
						{
							if(P[i][CELL_CELL[k][j+1]] < P_c_min)
								P_c_min = P[i][CELL_CELL[k][j+1]];
							else if(P[i][CELL_CELL[k][j+1]] > P_c_max)
								P_c_max = P[i][CELL_CELL[k][j+1]];
						}
					FAI_P[k] = DBL_MAX;
					for(j = 0; j < CELL_POINT[k][0]; ++j)
						{
							P_c_x_p = P[i][k] + grad_P_x[k] * (X[i][CELL_POINT[k][j+1]] - X_c[k]) + grad_P_y[k] * (Y[i][CELL_POINT[k][j+1]] - Y_c[k]);
							if (fabs(P_c_x_p - P[i][k]) <eps)
								FAI_P[k] = (FAI_P[k] < 1) ? FAI_P[k] : 1;	
							else if((P_c_x_p - P[i][k]) > 0)
								FAI_P[k] = (FAI_P[k] < (*miu)((P_c_max - P[i][k])/(P_c_x_p - P[i][k]))) ? FAI_P[k] : (*miu)((P_c_max - P[i][k])/(P_c_x_p - P[i][k]));
							else
								FAI_P[k] = (FAI_P[k] < (*miu)((P_c_min - P[i][k])/(P_c_x_p - P[i][k]))) ? FAI_P[k] : (*miu)((P_c_min - P[i][k])/(P_c_x_p - P[i][k]));
						}
				}

			for(k = 0; k < NUM_CELL; ++k)
				{
					M_c[0][0] = 0;
					M_c[0][1] = 0;
					M_c[1][0] = 0;
					M_c[1][1] = 0;
					grad_U_x[k] = 0;
					grad_U_y[k] = 0;
					
					for(j = 0; j < CELL_CELL[k][0]; ++j)
						{
							M_c[0][0] = M_c[0][0] + (X_c[CELL_CELL[k][j+1]] - X_c[k]) * (X_c[CELL_CELL[k][j+1]] - X_c[k]);
							M_c[0][1] = M_c[0][1] + (X_c[CELL_CELL[k][j+1]] - X_c[k]) * (Y_c[CELL_CELL[k][j+1]] - Y_c[k]);
							M_c[1][0] = M_c[1][0] + (Y_c[CELL_CELL[k][j+1]] - Y_c[k]) * (X_c[CELL_CELL[k][j+1]] - X_c[k]);
							M_c[1][1] = M_c[1][1] + (Y_c[CELL_CELL[k][j+1]] - Y_c[k]) * (Y_c[CELL_CELL[k][j+1]] - Y_c[k]);
							grad_U_x[k] = grad_U_x[k] +  (U[i][CELL_CELL[k][j+1]] - U[i][k]) * (X_c[CELL_CELL[k][j+1]] - X_c[k]);
							grad_U_y[k] = grad_U_y[k] +  (U[i][CELL_CELL[k][j+1]] - U[i][k]) * (Y_c[CELL_CELL[k][j+1]] - Y_c[k]);
						}
					
					M_c_i[0][0] =   M_c[1][1]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
					M_c_i[0][1] = - M_c[0][1]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
					M_c_i[1][0] = - M_c[1][0]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
					M_c_i[1][1] =   M_c[0][0]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
					temU_grad_U_x = M_c_i[0][0] * grad_U_x[k] + M_c_i[0][1] * grad_U_y[k];
					temU_grad_U_y = M_c_i[1][0] * grad_U_x[k] + M_c_i[1][1] * grad_U_y[k];
					grad_U_x[k] = temU_grad_U_x;
					grad_U_y[k] = temU_grad_U_y;						
				}

			for(k = 0; k < NUM_CELL; ++k)
				{
					U_c_min = U[i][k];
					U_c_max = U[i][k];
					for(j = 0; j < CELL_CELL[k][0]; ++j)
						{
							if(U[i][CELL_CELL[k][j+1]] < U_c_min)
								U_c_min = U[i][CELL_CELL[k][j+1]];
							else if(U[i][CELL_CELL[k][j+1]] > U_c_max)
								U_c_max = U[i][CELL_CELL[k][j+1]];
						}
					FAI_U[k] = DBL_MAX;
					for(j = 0; j < CELL_POINT[k][0]; ++j)
						{
							U_c_x_U = U[i][k] + grad_U_x[k] * (X[i][CELL_POINT[k][j+1]] - X_c[k]) + grad_U_y[k] * (Y[i][CELL_POINT[k][j+1]] - Y_c[k]);
							if (fabs(U_c_x_U - U[i][k]) <eps)
								FAI_U[k] = (FAI_U[k] < 1) ? FAI_U[k] : 1;	
							else if((U_c_x_U - U[i][k]) > 0)
								FAI_U[k] = (FAI_U[k] < (*miu)((U_c_max - U[i][k])/(U_c_x_U - U[i][k]))) ? FAI_U[k] : (*miu)((U_c_max - U[i][k])/(U_c_x_U - U[i][k]));
							else
								FAI_U[k] = (FAI_U[k] < (*miu)((U_c_min - U[i][k])/(U_c_x_U - U[i][k]))) ? FAI_U[k] : (*miu)((U_c_min - U[i][k])/(U_c_x_U - U[i][k]));
						}
				}


			for(k = 0; k < NUM_CELL; ++k)
				{
					M_c[0][0] = 0;
					M_c[0][1] = 0;
					M_c[1][0] = 0;
					M_c[1][1] = 0;
					grad_V_x[k] = 0;
					grad_V_y[k] = 0;
					
					for(j = 0; j < CELL_CELL[k][0]; ++j)
						{
							M_c[0][0] = M_c[0][0] + (X_c[CELL_CELL[k][j+1]] - X_c[k]) * (X_c[CELL_CELL[k][j+1]] - X_c[k]);
							M_c[0][1] = M_c[0][1] + (X_c[CELL_CELL[k][j+1]] - X_c[k]) * (Y_c[CELL_CELL[k][j+1]] - Y_c[k]);
							M_c[1][0] = M_c[1][0] + (Y_c[CELL_CELL[k][j+1]] - Y_c[k]) * (X_c[CELL_CELL[k][j+1]] - X_c[k]);
							M_c[1][1] = M_c[1][1] + (Y_c[CELL_CELL[k][j+1]] - Y_c[k]) * (Y_c[CELL_CELL[k][j+1]] - Y_c[k]);
							grad_V_x[k] = grad_V_x[k] +  (V[i][CELL_CELL[k][j+1]] - V[i][k]) * (X_c[CELL_CELL[k][j+1]] - X_c[k]);
							grad_V_y[k] = grad_V_y[k] +  (V[i][CELL_CELL[k][j+1]] - V[i][k]) * (Y_c[CELL_CELL[k][j+1]] - Y_c[k]);
						}
					
					M_c_i[0][0] =   M_c[1][1]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
					M_c_i[0][1] = - M_c[0][1]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
					M_c_i[1][0] = - M_c[1][0]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
					M_c_i[1][1] =   M_c[0][0]/(M_c[0][0]*M_c[1][1] - M_c[0][1]*M_c[1][0]);
					temV_grad_V_x = M_c_i[0][0] * grad_V_x[k] + M_c_i[0][1] * grad_V_y[k];
					temV_grad_V_y = M_c_i[1][0] * grad_V_x[k] + M_c_i[1][1] * grad_V_y[k];
					grad_V_x[k] = temV_grad_V_x;
					grad_V_y[k] = temV_grad_V_y;						
				}

			for(k = 0; k < NUM_CELL; ++k)
				{
					V_c_min = V[i][k];
					V_c_max = V[i][k];
					for(j = 0; j < CELL_CELL[k][0]; ++j)
						{
							if(V[i][CELL_CELL[k][j+1]] < V_c_min)
								V_c_min = V[i][CELL_CELL[k][j+1]];
							else if(V[i][CELL_CELL[k][j+1]] > V_c_max)
								V_c_max = V[i][CELL_CELL[k][j+1]];
						}
					FAI_V[k] = DBL_MAX;
					for(j = 0; j < CELL_POINT[k][0]; ++j)
						{
							V_c_x_V = V[i][k] + grad_V_x[k] * (X[i][CELL_POINT[k][j+1]] - X_c[k]) + grad_V_y[k] * (Y[i][CELL_POINT[k][j+1]] - Y_c[k]);
							if (fabs(V_c_x_V - V[i][k]) <eps)
								FAI_V[k] = (FAI_V[k] < 1) ? FAI_V[k] : 1;	
							else if((V_c_x_V - V[i][k]) > 0)
								FAI_V[k] = (FAI_V[k] < (*miu)((V_c_max - V[i][k])/(V_c_x_V - V[i][k]))) ? FAI_V[k] : (*miu)((V_c_max - V[i][k])/(V_c_x_V - V[i][k]));
							else
								FAI_V[k] = (FAI_V[k] < (*miu)((V_c_min - V[i][k])/(V_c_x_V - V[i][k]))) ? FAI_V[k] : (*miu)((V_c_min - V[i][k])/(V_c_x_V - V[i][k]));
						}
				}


			for(k = 0; k < NUM_CELL; ++k)
				{
					for(j = 0; j < CELL_POINT[k][0]; ++j)
						{
							P_c_lim[k][j] = P[i][k] + FAI_P[k] * (grad_P_x[k] * (X[i][CELL_POINT[k][j+1]] - X_c[k]) + grad_P_y[k] * (Y[i][CELL_POINT[k][j+1]] - Y_c[k]));
							U_c_lim[k][j] = U[i][k] + FAI_U[k] * (grad_U_x[k] * (X[i][CELL_POINT[k][j+1]] - X_c[k]) + grad_U_y[k] * (Y[i][CELL_POINT[k][j+1]] - Y_c[k]));
							V_c_lim[k][j] = V[i][k] + FAI_V[k] * (grad_V_x[k] * (X[i][CELL_POINT[k][j+1]] - X_c[k]) + grad_V_y[k] * (Y[i][CELL_POINT[k][j+1]] - Y_c[k]));	
						}						
				}





							
			for(k = 0; k < NUM_CELL; ++k)
				{
					for(j = 0; j < CELL_POINT[k][0]; ++j)
						{						
						 if(j == CELL_POINT[k][0]-1) 
							 {
								 p_p=CELL_POINT[k][1];
								 p_n=CELL_POINT[k][j];
							 }						 
						 else if(j)
							 {
								 p_p=CELL_POINT[k][j+2];
								 p_n=CELL_POINT[k][j];
							 } 
						 else 
							 {
								 p_p=CELL_POINT[k][j+2];
								 p_n=CELL_POINT[k][CELL_POINT[k][0]];
							 }
						 							 
						 l_pc[k][j] = 0.5 * sqrt((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p]));
						 n_pc_x[k][j] = 0.5 * (Y[i][p_p] - Y[i][p_n]) / l_pc[k][j];						
						 n_pc_y[k][j] = 0.5 * (X[i][p_n] - X[i][p_p]) / l_pc[k][j];						 
						}
				}
		
			for(k = 0; k < NUM_CELL; ++k)
				{
					for(j = 0; j < CELL_POINT[k][0]; ++j)
						{				
							if(j == CELL_POINT[k][0]-1) 
								{
									p_p=CELL_POINT[k][1];
									p_n=CELL_POINT[k][j+1];
								}							
							else
								{
									p_p=CELL_POINT[k][j+2];
									p_n=CELL_POINT[k][j+1];
								} 
						
							posi_l_pc[k][j] = 0.5 * sqrt((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p]));
							posi_n_pc_x[k][j] = 0.5 * (Y[i][p_p] - Y[i][p_n]) / posi_l_pc[k][j];
							posi_n_pc_y[k][j] = 0.5 * (X[i][p_n] - X[i][p_p]) / posi_l_pc[k][j];
						}	
				}
		
			for(k = 0; k < NUM_CELL; ++k)
				{
					for(j = 0; j < CELL_POINT[k][0]; ++j)
						{																
							if(j)
								{
									p_p=CELL_POINT[k][j+1];
									p_n=CELL_POINT[k][j];
								}
							else
								{
									p_p=CELL_POINT[k][j+1];
									p_n=CELL_POINT[k][CELL_POINT[k][0]];
								}
						
							nega_l_pc[k][j] = 0.5 * sqrt((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p]));
							nega_n_pc_x[k][j] = 0.5 * (Y[i][p_p] - Y[i][p_n]) / nega_l_pc[k][j];
							nega_n_pc_y[k][j] = 0.5 * (X[i][p_n] - X[i][p_p]) / nega_l_pc[k][j];
						}
				}
		
			for(k = 0; k < NUM_BOUNDARY; ++k)
				{			
					if(k == NUM_BOUNDARY-1) 
						{
							p_p=BOUNDARY_POINT[0][0];
							p_n=BOUNDARY_POINT[0][k];
						}					
					else
						{
							p_p=BOUNDARY_POINT[0][k+1];
							p_n=BOUNDARY_POINT[0][k];
						}							  
					posi_l_p[k] = 0.5 * sqrt((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p]));
					posi_n_p_x[k] = 0.5 * (Y[i][p_p] - Y[i][p_n]) / posi_l_p[k];
					posi_n_p_y[k] = 0.5 * (X[i][p_n] - X[i][p_p]) / posi_l_p[k];
				}

			for(k = 0; k < NUM_BOUNDARY; ++k)
				{			
					if(k)
						{
							p_p=BOUNDARY_POINT[0][k];
							p_n=BOUNDARY_POINT[0][k-1];
						} 
					else 
						{
							p_p=BOUNDARY_POINT[0][k];
							p_n=BOUNDARY_POINT[0][NUM_BOUNDARY-1];
						}							  
					nega_l_p[k] = 0.5 * sqrt((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p]));
					nega_n_p_x[k] = 0.5 * (Y[i][p_p] - Y[i][p_n]) / nega_l_p[k];
					nega_n_p_y[k] = 0.5 * (X[i][p_n] - X[i][p_p]) / nega_l_p[k];
				}

			
			
    for(k = 0; k < NUM_CELL; ++k)
	  {
		  for(j = 0; j < CELL_POINT[k][0]; ++j)
			  {

				  GLACE(M_pc, gamma[k], P[i][k], RHO[i][k], l_pc[k][j], n_pc_x[k][j], n_pc_y[k][j]);
                                 M_pc_00[k][j] = M_pc[0][0];
                                 M_pc_01[k][j] = M_pc[0][1];
                                 M_pc_10[k][j] = M_pc[1][0];
                                 M_pc_11[k][j] = M_pc[1][1];
			  }
	  }	  	  
	
	
	for(l = 0; l < NUM_POINT; ++l)
		{
			U_p_x[l]=0.0;
			U_p_y[l]=0.0;
			M_p[0][0]=0.0;
			M_p[0][1]=0.0;
			M_p[1][0]=0.0;
			M_p[1][1]=0.0;
			A_p[0][0]=0.0;
			A_p[0][1]=0.0;
			A_p[0][2]=0.0;
			A_p[1][0]=0.0;
			A_p[1][1]=0.0;
			A_p[1][2]=0.0;
			A_p[2][0]=0.0;
			A_p[2][1]=0.0;
			A_p[2][2]=0.0;
			PAI_p[l]=0;

			for(NUM_PP = 0; NUM_PP < POINT_POINT[l][0]; ++NUM_PP)
					{
						k=POINT_POINT[l][NUM_PP*2+1];
						j=POINT_POINT[l][NUM_PP*2+2];
										U_p_x[l]=U_p_x[l] + l_pc[k][j]*n_pc_x[k][j]*P_c_lim[k][j] + M_pc_00[k][j]*U_c_lim[k][j] + M_pc_01[k][j]*V_c_lim[k][j];
										U_p_y[l]=U_p_y[l] + l_pc[k][j]*n_pc_y[k][j]*P_c_lim[k][j] + M_pc_10[k][j]*U_c_lim[k][j] + M_pc_11[k][j]*V_c_lim[k][j];
										M_p[0][0]=M_p[0][0] +  M_pc_00[k][j];
										M_p[0][1]=M_p[0][1] +  M_pc_01[k][j];
										M_p[1][0]=M_p[1][0] +  M_pc_10[k][j];
										M_p[1][1]=M_p[1][1] +  M_pc_11[k][j];
										A_p[0][2]=A_p[0][2] + l_pc[k][j]*n_pc_x[k][j];
										A_p[1][2]=A_p[1][2] + l_pc[k][j]*n_pc_y[k][j];
					}
								
			if(index[l]==-1)
				{
					M_p_i[0][0] =   M_p[1][1]/(M_p[0][0]*M_p[1][1] - M_p[0][1]*M_p[1][0]);
					M_p_i[0][1] = - M_p[0][1]/(M_p[0][0]*M_p[1][1] - M_p[0][1]*M_p[1][0]);
					M_p_i[1][0] = - M_p[1][0]/(M_p[0][0]*M_p[1][1] - M_p[0][1]*M_p[1][0]);
					M_p_i[1][1] =   M_p[0][0]/(M_p[0][0]*M_p[1][1] - M_p[0][1]*M_p[1][0]);
					temp_U_p_x = M_p_i[0][0]*U_p_x[l] + M_p_i[0][1]*U_p_y[l];
					temp_U_p_y = M_p_i[1][0]*U_p_x[l] + M_p_i[1][1]*U_p_y[l];
					U_p_x[l] = temp_U_p_x;
					U_p_y[l] = temp_U_p_y;
				}
			else if(fabs(fabs(posi_n_p_x[index[l]]*nega_n_p_y[index[l]]) - fabs(posi_n_p_y[index[l]]*nega_n_p_x[index[l]])) > eps)
				{

					M_p[0][0] = posi_n_p_x[index[l]];
					M_p[0][1] = posi_n_p_y[index[l]];
					M_p[1][0] = nega_n_p_x[index[l]];
					M_p[1][1] = nega_n_p_y[index[l]];
					M_p_i[0][0] =   M_p[1][1]/(M_p[0][0]*M_p[1][1] - M_p[0][1]*M_p[1][0]);
					M_p_i[0][1] = - M_p[0][1]/(M_p[0][0]*M_p[1][1] - M_p[0][1]*M_p[1][0]);
					M_p_i[1][0] = - M_p[1][0]/(M_p[0][0]*M_p[1][1] - M_p[0][1]*M_p[1][0]);
					M_p_i[1][1] =   M_p[0][0]/(M_p[0][0]*M_p[1][1] - M_p[0][1]*M_p[1][0]);
					U_p_x[l] = M_p_i[0][0]*NORMAL_VELOCITY[0][index[l]] + M_p_i[0][1]*NORMAL_VELOCITY[1][index[l]];
					U_p_y[l] = M_p_i[1][0]*NORMAL_VELOCITY[0][index[l]] + M_p_i[1][1]*NORMAL_VELOCITY[1][index[l]];
				}
			else
				{

					A_p[0][0]=M_p[0][0];
					A_p[0][1]=M_p[0][1];
					A_p[1][0]=M_p[1][0];
					A_p[1][1]=M_p[1][1];
					A_p[2][0]=A_p[0][2];
					A_p[2][1]=A_p[1][2];
					A_p_det = A_p[2][0]*A_p[0][1]*A_p[1][2]-A_p[2][0]*A_p[0][2]*A_p[1][1]-A_p[1][0]*A_p[0][1]*A_p[2][2]+A_p[1][0]*A_p[0][2]*A_p[2][1]+A_p[0][0]*A_p[1][1]*A_p[2][2]-A_p[0][0]*A_p[1][2]*A_p[2][1];
					A_p_i[0][0] =   (A_p[1][1]*A_p[2][2]-A_p[1][2]*A_p[2][1]) / A_p_det;
					A_p_i[0][1] = - (A_p[0][1]*A_p[2][2]-A_p[0][2]*A_p[2][1]) / A_p_det;
					A_p_i[0][2] =   (A_p[0][1]*A_p[1][2]-A_p[0][2]*A_p[1][1]) / A_p_det;
					A_p_i[1][0] = - (A_p[1][0]*A_p[2][2]-A_p[2][0]*A_p[1][2]) / A_p_det;
					A_p_i[1][1] =   (A_p[0][0]*A_p[2][2]-A_p[2][0]*A_p[0][2]) / A_p_det;
					A_p_i[1][2] = - (A_p[0][0]*A_p[1][2]-A_p[1][0]*A_p[0][2]) / A_p_det;
					A_p_i[2][0] =   (A_p[1][0]*A_p[2][1]-A_p[2][0]*A_p[1][1]) / A_p_det;
					A_p_i[2][1] = - (A_p[0][0]*A_p[2][1]-A_p[2][0]*A_p[0][1]) / A_p_det;
					A_p_i[2][2] =   (A_p[0][0]*A_p[1][1]-A_p[1][0]*A_p[0][1]) / A_p_det;
					temp_U_p_x = A_p_i[0][0]*U_p_x[l] + A_p_i[0][1]*U_p_y[l] + A_p_i[0][2]*(NORMAL_VELOCITY[0][index[l]]*posi_l_p[index[l]] + NORMAL_VELOCITY[1][index[l]]*nega_l_p[index[l]]);
					temp_U_p_y = A_p_i[1][0]*U_p_x[l] + A_p_i[1][1]*U_p_y[l] + A_p_i[1][2]*(NORMAL_VELOCITY[0][index[l]]*posi_l_p[index[l]] + NORMAL_VELOCITY[1][index[l]]*nega_l_p[index[l]]);
					  PAI_p[l] = A_p_i[2][0]*U_p_x[l] + A_p_i[2][1]*U_p_y[l] + A_p_i[2][2]*(NORMAL_VELOCITY[0][index[l]]*posi_l_p[index[l]] + NORMAL_VELOCITY[1][index[l]]*nega_l_p[index[l]]);
					U_p_x[l] = temp_U_p_x;
					U_p_y[l] = temp_U_p_y;
					
				}
			
		}

	for(k = 0; k < NUM_CELL; ++k)
		{
			for(j = 0; j < CELL_POINT[k][0]; ++j)
				{
					F_pc_x[k][j]=l_pc[k][j]*n_pc_x[k][j]*P_c_lim[k][j] - M_pc_00[k][j]*(U_p_x[CELL_POINT[k][j+1]] - U_c_lim[k][j]) - M_pc_01[k][j]*(U_p_y[CELL_POINT[k][j+1]] - V_c_lim[k][j]);
					F_pc_y[k][j]=l_pc[k][j]*n_pc_y[k][j]*P_c_lim[k][j] - M_pc_10[k][j]*(U_p_x[CELL_POINT[k][j+1]] - U_c_lim[k][j]) - M_pc_11[k][j]*(U_p_y[CELL_POINT[k][j+1]] - V_c_lim[k][j]);

				}
		}

			for(k = 0; k < NUM_BOUNDARY; ++k)
				{			
					if(k == NUM_BOUNDARY-1) 
						{
							p_p=BOUNDARY_POINT[0][0];
							p_n=BOUNDARY_POINT[0][k];
						}					
					else
						{
							p_p=BOUNDARY_POINT[0][k+1];
							p_n=BOUNDARY_POINT[0][k];
						}							  					
					dt_posi_n_p_x[k] = ((U_p_y[p_p] - U_p_y[p_n])*((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) - (Y[i][p_p] - Y[i][p_n]) * (X[i][p_p] - X[i][p_n]) * (U_p_x[p_p] - U_p_x[p_n]) - (Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) * (U_p_y[p_p] - U_p_y[p_n])) / sqrt((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) / ((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p]));
					dt_posi_n_p_y[k] = ((U_p_x[p_n] - U_p_x[p_p])*((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) + (X[i][p_p] - X[i][p_n]) * (X[i][p_p] - X[i][p_n]) * (U_p_x[p_p] - U_p_x[p_n]) + (X[i][p_p] - X[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) * (U_p_y[p_p] - U_p_y[p_n])) / sqrt((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) / ((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p]));
				}

			for(k = 0; k < NUM_BOUNDARY; ++k)
				{			
					if(k)
						{
							p_p=BOUNDARY_POINT[0][k];
							p_n=BOUNDARY_POINT[0][k-1];
						} 
					else 
						{
							p_p=BOUNDARY_POINT[0][k];
							p_n=BOUNDARY_POINT[0][NUM_BOUNDARY-1];
						}							  
					
					dt_nega_n_p_x[k] = ((U_p_y[p_p] - U_p_y[p_n])*((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) - (Y[i][p_p] - Y[i][p_n]) * (X[i][p_p] - X[i][p_n]) * (U_p_x[p_p] - U_p_x[p_n]) - (Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) * (U_p_y[p_p] - U_p_y[p_n])) / sqrt((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) / ((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p]));
					dt_nega_n_p_y[k] = ((U_p_x[p_n] - U_p_x[p_p])*((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) + (X[i][p_p] - X[i][p_n]) * (X[i][p_p] - X[i][p_n]) * (U_p_x[p_p] - U_p_x[p_n]) + (X[i][p_p] - X[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) * (U_p_y[p_p] - U_p_y[p_n])) / sqrt((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) / ((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p]));
				}

	for(k = 0; k < NUM_CELL; ++k)
				{
					for(j = 0; j < CELL_POINT[k][0]; ++j)
						{						
						 if(j == CELL_POINT[k][0]-1) 
							 {
								 p_p=CELL_POINT[k][1];
								 p_n=CELL_POINT[k][j];
							 }						 
						 else if(j)
							 {
								 p_p=CELL_POINT[k][j+2];
								 p_n=CELL_POINT[k][j];
							 } 
						 else 
							 {
								 p_p=CELL_POINT[k][j+2];
								 p_n=CELL_POINT[k][CELL_POINT[k][0]];
							 }
						 							 						 
						 dt_l_n_pc_x[k][j] = 0.5 * (U_p_y[p_p] - U_p_y[p_n]);						
						 dt_l_n_pc_y[k][j] = 0.5 * (U_p_x[p_n] - U_p_x[p_p]);
						 dt_n_pc_x[k][j]   = ((U_p_y[p_p] - U_p_y[p_n])*((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) - (Y[i][p_p] - Y[i][p_n]) * (X[i][p_p] - X[i][p_n]) * (U_p_x[p_p] - U_p_x[p_n]) - (Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) * (U_p_y[p_p] - U_p_y[p_n])) / sqrt((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) / ((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p]));
						 dt_n_pc_y[k][j]   = ((U_p_x[p_n] - U_p_x[p_p])*((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) + (X[i][p_p] - X[i][p_n]) * (X[i][p_p] - X[i][p_n]) * (U_p_x[p_p] - U_p_x[p_n]) + (X[i][p_p] - X[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) * (U_p_y[p_p] - U_p_y[p_n])) / sqrt((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) / ((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p]));
						}
				}
		
			for(k = 0; k < NUM_CELL; ++k)
				{
					for(j = 0; j < CELL_POINT[k][0]; ++j)
						{				
							if(j == CELL_POINT[k][0]-1) 
								{
									p_p=CELL_POINT[k][1];
									p_n=CELL_POINT[k][j+1];
								}							
							else
								{
									p_p=CELL_POINT[k][j+2];
									p_n=CELL_POINT[k][j+1];
								} 
												
							dt_posi_l_n_pc_x[k][j] = 0.5 * (U_p_y[p_p] - U_p_y[p_n]);
							dt_posi_l_n_pc_y[k][j] = 0.5 * (U_p_x[p_n] - U_p_x[p_p]);
							dt_posi_n_pc_x[k][j]   = ((U_p_y[p_p] - U_p_y[p_n])*((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) - (Y[i][p_p] - Y[i][p_n]) * (X[i][p_p] - X[i][p_n]) * (U_p_x[p_p] - U_p_x[p_n]) - (Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) * (U_p_y[p_p] - U_p_y[p_n])) / sqrt((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) / ((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p]));
							dt_posi_n_pc_y[k][j]   = ((U_p_x[p_n] - U_p_x[p_p])*((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) + (X[i][p_p] - X[i][p_n]) * (X[i][p_p] - X[i][p_n]) * (U_p_x[p_p] - U_p_x[p_n]) + (X[i][p_p] - X[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) * (U_p_y[p_p] - U_p_y[p_n])) / sqrt((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) / ((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p]));
						}	
				}
		
			for(k = 0; k < NUM_CELL; ++k)
				{
					for(j = 0; j < CELL_POINT[k][0]; ++j)
						{																
							if(j)
								{
									p_p=CELL_POINT[k][j+1];
									p_n=CELL_POINT[k][j];
								}
							else
								{
									p_p=CELL_POINT[k][j+1];
									p_n=CELL_POINT[k][CELL_POINT[k][0]];
								}
						
							dt_nega_l_n_pc_x[k][j] = 0.5 * (U_p_y[p_p] - U_p_y[p_n]);
							dt_nega_l_n_pc_y[k][j] = 0.5 * (U_p_x[p_n] - U_p_x[p_p]);
							dt_nega_n_pc_x[k][j]   = ((U_p_y[p_p] - U_p_y[p_n])*((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) - (Y[i][p_p] - Y[i][p_n]) * (X[i][p_p] - X[i][p_n]) * (U_p_x[p_p] - U_p_x[p_n]) - (Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) * (U_p_y[p_p] - U_p_y[p_n])) / sqrt((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) / ((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p]));
							dt_nega_n_pc_y[k][j]   = ((U_p_x[p_n] - U_p_x[p_p])*((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) + (X[i][p_p] - X[i][p_n]) * (X[i][p_p] - X[i][p_n]) * (U_p_x[p_p] - U_p_x[p_n]) + (X[i][p_p] - X[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) * (U_p_y[p_p] - U_p_y[p_n])) / sqrt((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p])) / ((Y[i][p_p] - Y[i][p_n]) * (Y[i][p_p] - Y[i][p_n]) + (X[i][p_n] - X[i][p_p]) * (X[i][p_n] - X[i][p_p]));
						}
				}

			for(k = 0; k < NUM_CELL; ++k)
				{
					for(j = 0; j < CELL_POINT[k][0]; ++j)
						{
							
							GLACE_dt(M_pc, gamma[k], P[i][k], RHO[i][k], l_pc[k][j], n_pc_x[k][j], n_pc_y[k][j],  dt_l_n_pc_x[k][j], dt_l_n_pc_y[k][j], dt_n_pc_x[k][j], dt_n_pc_y[k][j]);
							dt_M_pc_00[k][j] = M_pc[0][0];
							dt_M_pc_01[k][j] = M_pc[0][1];
							dt_M_pc_10[k][j] = M_pc[1][0];
							dt_M_pc_11[k][j] = M_pc[1][1];
						}
				}	  
			
	for(l = 0; l < NUM_POINT; ++l)
		{
			dt_U_p_x[l]=0.0;
			dt_U_p_y[l]=0.0;
			M_p[0][0]=0.0;
			M_p[0][1]=0.0;
			M_p[1][0]=0.0;
			M_p[1][1]=0.0;
			A_p[0][0]=0.0;
			A_p[0][1]=0.0;
			A_p[0][2]=0.0;
			A_p[1][0]=0.0;
			A_p[1][1]=0.0;
			A_p[1][2]=0.0;
			A_p[2][0]=0.0;
			A_p[2][1]=0.0;
			A_p[2][2]=0.0;
			b_temp_x=0.0;
			b_temp_y=0.0;
			b_temp_z=0.0;

			for(NUM_PP = 0; NUM_PP < POINT_POINT[l][0]; ++NUM_PP)
					{
						k=POINT_POINT[l][NUM_PP*2+1];
						j=POINT_POINT[l][NUM_PP*2+2];
										dt_U_p_x[l]=dt_U_p_x[l] - 1.0/RHO[i][k]*(M_pc_00[k][j]*grad_P_x[k] + M_pc_01[k][j]*grad_P_y[k] + RHO[i][k]*a[k]*RHO[i][k]*a[k]*l_pc[k][j]*n_pc_x[k][j]*(grad_U_x[k]+grad_V_y[k])) + P[i][k]*dt_l_n_pc_x[k][j] - dt_M_pc_00[k][j]*(U_p_x[l]-U[i][k]) - dt_M_pc_01[k][j]*(U_p_y[l]-V[i][k]);
										dt_U_p_y[l]=dt_U_p_y[l] - 1.0/RHO[i][k]*(M_pc_10[k][j]*grad_P_x[k] + M_pc_11[k][j]*grad_P_y[k] + RHO[i][k]*a[k]*RHO[i][k]*a[k]*l_pc[k][j]*n_pc_y[k][j]*(grad_U_x[k]+grad_V_y[k])) + P[i][k]*dt_l_n_pc_y[k][j] - dt_M_pc_10[k][j]*(U_p_x[l]-U[i][k]) - dt_M_pc_11[k][j]*(U_p_y[l]-V[i][k]);
										b_temp_x=b_temp_x - dt_l_n_pc_x[k][j]*PAI_p[l];
										b_temp_y=b_temp_y - dt_l_n_pc_y[k][j]*PAI_p[l];
										b_temp_z=b_temp_z - dt_l_n_pc_x[k][j]*U_p_x[l] - dt_l_n_pc_y[k][j]*U_p_y[l];
										M_p[0][0]=M_p[0][0] +  M_pc_00[k][j];
										M_p[0][1]=M_p[0][1] +  M_pc_01[k][j];
										M_p[1][0]=M_p[1][0] +  M_pc_10[k][j];
										M_p[1][1]=M_p[1][1] +  M_pc_11[k][j];
										A_p[0][2]=A_p[0][2] + l_pc[k][j]*n_pc_x[k][j];
										A_p[1][2]=A_p[1][2] + l_pc[k][j]*n_pc_y[k][j];
								
					}
								
			if(index[l]==-1)
				{
					M_p_i[0][0] =   M_p[1][1]/(M_p[0][0]*M_p[1][1] - M_p[0][1]*M_p[1][0]);
					M_p_i[0][1] = - M_p[0][1]/(M_p[0][0]*M_p[1][1] - M_p[0][1]*M_p[1][0]);
					M_p_i[1][0] = - M_p[1][0]/(M_p[0][0]*M_p[1][1] - M_p[0][1]*M_p[1][0]);
					M_p_i[1][1] =   M_p[0][0]/(M_p[0][0]*M_p[1][1] - M_p[0][1]*M_p[1][0]);
					temp_U_p_x = M_p_i[0][0]*dt_U_p_x[l] + M_p_i[0][1]*dt_U_p_y[l];
					temp_U_p_y = M_p_i[1][0]*dt_U_p_x[l] + M_p_i[1][1]*dt_U_p_y[l];
					dt_U_p_x[l] = temp_U_p_x;
					dt_U_p_y[l] = temp_U_p_y;
				}
			else if(fabs(fabs(posi_n_p_x[index[l]]*nega_n_p_y[index[l]]) - fabs(posi_n_p_y[index[l]]*nega_n_p_x[index[l]])) > eps)
				{

					M_p[0][0] = posi_n_p_x[index[l]];
					M_p[0][1] = posi_n_p_y[index[l]];
					M_p[1][0] = nega_n_p_x[index[l]];
					M_p[1][1] = nega_n_p_y[index[l]];
					M_p_i[0][0] =   M_p[1][1]/(M_p[0][0]*M_p[1][1] - M_p[0][1]*M_p[1][0]);
					M_p_i[0][1] = - M_p[0][1]/(M_p[0][0]*M_p[1][1] - M_p[0][1]*M_p[1][0]);
					M_p_i[1][0] = - M_p[1][0]/(M_p[0][0]*M_p[1][1] - M_p[0][1]*M_p[1][0]);
					M_p_i[1][1] =   M_p[0][0]/(M_p[0][0]*M_p[1][1] - M_p[0][1]*M_p[1][0]);
					dt_U_p_x[l] = M_p_i[0][0]*(- U_p_x[l]*dt_posi_n_p_x[index[l]] - U_p_y[l]*dt_posi_n_p_y[index[l]]) + M_p_i[0][1]*(- U_p_x[l]*dt_nega_n_p_x[index[l]] - U_p_y[l]*dt_nega_n_p_y[index[l]]);
					dt_U_p_y[l] = M_p_i[1][0]*(- U_p_x[l]*dt_posi_n_p_x[index[l]] - U_p_y[l]*dt_posi_n_p_y[index[l]]) + M_p_i[1][1]*(- U_p_x[l]*dt_nega_n_p_x[index[l]] - U_p_y[l]*dt_nega_n_p_y[index[l]]);
				}
			else
				{

					A_p[0][0]=M_p[0][0];
					A_p[0][1]=M_p[0][1];
					A_p[1][0]=M_p[1][0];
					A_p[1][1]=M_p[1][1];
					A_p[2][0]=A_p[0][2];
					A_p[2][1]=A_p[1][2];
					A_p_det = A_p[2][0]*A_p[0][1]*A_p[1][2]-A_p[2][0]*A_p[0][2]*A_p[1][1]-A_p[1][0]*A_p[0][1]*A_p[2][2]+A_p[1][0]*A_p[0][2]*A_p[2][1]+A_p[0][0]*A_p[1][1]*A_p[2][2]-A_p[0][0]*A_p[1][2]*A_p[2][1];
					A_p_i[0][0] =   (A_p[1][1]*A_p[2][2]-A_p[1][2]*A_p[2][1]) / A_p_det;
					A_p_i[0][1] = - (A_p[0][1]*A_p[2][2]-A_p[0][2]*A_p[2][1]) / A_p_det;
					A_p_i[0][2] =   (A_p[0][1]*A_p[1][2]-A_p[0][2]*A_p[1][1]) / A_p_det;
					A_p_i[1][0] = - (A_p[1][0]*A_p[2][2]-A_p[2][0]*A_p[1][2]) / A_p_det;
					A_p_i[1][1] =   (A_p[0][0]*A_p[2][2]-A_p[2][0]*A_p[0][2]) / A_p_det;
					A_p_i[1][2] = - (A_p[0][0]*A_p[1][2]-A_p[1][0]*A_p[0][2]) / A_p_det;
					A_p_i[2][0] =   (A_p[1][0]*A_p[2][1]-A_p[2][0]*A_p[1][1]) / A_p_det;
					A_p_i[2][1] = - (A_p[0][0]*A_p[2][1]-A_p[2][0]*A_p[0][1]) / A_p_det;
					A_p_i[2][2] =   (A_p[0][0]*A_p[1][1]-A_p[1][0]*A_p[0][1]) / A_p_det;
					temp_U_p_x = A_p_i[0][0]*(dt_U_p_x[l]+b_temp_x) + A_p_i[0][1]*(dt_U_p_y[l]+b_temp_y) + A_p_i[0][2]*b_temp_z;
					temp_U_p_y = A_p_i[1][0]*(dt_U_p_x[l]+b_temp_x) + A_p_i[1][1]*(dt_U_p_y[l]+b_temp_y) + A_p_i[1][2]*b_temp_z;
					dt_U_p_x[l] = temp_U_p_x;
					dt_U_p_y[l] = temp_U_p_y;
					}
			
		}

			

	for(k = 0; k < NUM_CELL; ++k)
		{
			for(j = 0; j < CELL_POINT[k][0]; ++j)
				{
					dt_F_pc_x[k][j]= - a[k]*RHO[i][k]*a[k]*l_pc[k][j]*n_pc_x[k][j]*(grad_U_x[k]+grad_V_y[k]) - M_pc_00[k][j]*(grad_P_x[k]/RHO[i][k]+dt_U_p_x[CELL_POINT[k][j+1]]) - M_pc_01[k][j]*(grad_P_y[k]/RHO[i][k]+dt_U_p_y[CELL_POINT[k][j+1]]) + P[i][k]*dt_l_n_pc_x[k][j] - dt_M_pc_00[k][j]*(U_p_x[l]-U[i][k]) - dt_M_pc_01[k][j]*(U_p_y[i]-V[i][k]);						
					dt_F_pc_y[k][j]= - a[k]*RHO[i][k]*a[k]*l_pc[k][j]*n_pc_y[k][j]*(grad_U_x[k]+grad_V_y[k]) - M_pc_10[k][j]*(grad_P_x[k]/RHO[i][k]+dt_U_p_x[CELL_POINT[k][j+1]]) - M_pc_11[k][j]*(grad_P_y[k]/RHO[i][k]+dt_U_p_y[CELL_POINT[k][j+1]]) + P[i][k]*dt_l_n_pc_y[k][j] - dt_M_pc_10[k][j]*(U_p_x[l]-U[i][k]) - dt_M_pc_11[k][j]*(U_p_y[i]-V[i][k]);		
				}
		}

tau = 1.01*tau;
for(k = 0; k < NUM_CELL; ++k)
	{
		lambda[k] = DBL_MAX;
		for(j = 1; j <= CELL_POINT[k][0]; ++j)
			{
				for(l = j+1; l <= CELL_POINT[k][0]; ++l)
					{
						lambda[k] = min(lambda[k],sqrt((X[i][CELL_POINT[k][j]]-X[i][CELL_POINT[k][l]])*(X[i][CELL_POINT[k][j]]-X[i][CELL_POINT[k][l]]) + (Y[i][CELL_POINT[k][j]]-Y[i][CELL_POINT[k][l]])*(Y[i][CELL_POINT[k][j]]-Y[i][CELL_POINT[k][l]])));
					}
			}
		tau = min(tau,CFL*lambda[k]/a[k]);
		volume_temp = 0.0;
		dt_volume_temp = 0.0;
		for(j = 0; j < CELL_POINT[k][0]; ++j)
			  {
				  if(j == CELL_POINT[k][0]-1) 
					  {
						  p_p=CELL_POINT[k][1];
						  p_n=CELL_POINT[k][j+1];
					  }				  
				  else
					  {
						  p_p=CELL_POINT[k][j+2];
						  p_n=CELL_POINT[k][j+1];
					  } 
				  volume_temp = volume_temp + 0.5 * (X[i][p_n] * Y[i][p_p] - Y[i][p_n] * X[i][p_p]);
				  dt_volume_temp += l_pc[k][j]*n_pc_x[k][j]*U_p_x[CELL_POINT[k][j+1]] + l_pc[k][j]*n_pc_y[k][j]*U_p_y[CELL_POINT[k][j+1]];
			  }
		tau = min(tau,0.1*volume_temp/fabs(dt_volume_temp));
	}

	if (CFL<0.0)
		{
			tau=-CFL;
		}
	t_all += tau;
	
	if(tau < config[1]/N/5)
		{
			printf("The length of the time step is so small at step %d, t_all=%lf, tau=%lf.\n",i,t_all,tau);
			stop_step=1;
		}
	
	if(t_all > config[1])
		{
			printf("The time is enough at step %d.\n",i);
			tau = tau - (t_all-config[1]);
			stop_step=1;
		}	

	for(k = 0; k < NUM_CELL; ++k)
		{
			RHO[i+1][k]    = 0.0;
			U[i+1][k]      = 0.0;
			V[i+1][k]      = 0.0;
			P[i+1][k]      = 0.0;				
			for(j = 0; j < CELL_POINT[k][0]; ++j)
				{
					RHO[i+1][k] = RHO[i+1][k] + (l_pc[k][j]*n_pc_x[k][j]+tau/2.0*dt_l_n_pc_x[k][j])*(U_p_x[CELL_POINT[k][j+1]]+tau/2.0*dt_U_p_x[CELL_POINT[k][j+1]]) + (l_pc[k][j]*n_pc_y[k][j]+tau/2.0*dt_l_n_pc_y[k][j])*(U_p_y[CELL_POINT[k][j+1]]+tau/2.0*dt_U_p_y[CELL_POINT[k][j+1]]);
					     U[i+1][k] = U[i+1][k] - F_pc_x[k][j] - tau/2.0*dt_F_pc_x[k][j];
						 V[i+1][k] = V[i+1][k] - F_pc_y[k][j] - tau/2.0*dt_F_pc_y[k][j];
						 P[i+1][k] = P[i+1][k] - (F_pc_x[k][j] + tau/2.0*dt_F_pc_x[k][j])*(U_p_x[CELL_POINT[k][j+1]] + tau/2.0*dt_U_p_x[CELL_POINT[k][j+1]]) - (F_pc_y[k][j] + tau/2.0*dt_F_pc_y[k][j])*(U_p_y[CELL_POINT[k][j+1]] + tau/2.0*dt_U_p_y[CELL_POINT[k][j+1]]);
				}
			
			RHO[i+1][k] = 1.0/RHO[i][k] + tau/MASS[k]*RHO[i+1][k];
			RHO[i+1][k] = 1.0/RHO[i+1][k];
			     U[i+1][k] = U[i][k] + tau/MASS[k]*U[i+1][k];
			     V[i+1][k] = V[i][k] + tau/MASS[k]*V[i+1][k];
				 P[i+1][k] = P[i][k]/(gamma[k]-1.0)/RHO[i][k] + 0.5*(U[i][k]*U[i][k]+V[i][k]*V[i][k]) + tau/MASS[k]*P[i+1][k];
				 P[i+1][k] = (P[i+1][k] - 0.5*(U[i+1][k]*U[i+1][k]+V[i+1][k]*V[i+1][k]))*(gamma[k]-1.0)*RHO[i+1][k];			 	 
		}
	

	for(k = 0; k < NUM_POINT; ++k)
		{				
			X[i+1][k] = X[i][k] + tau*(U_p_x[k]+tau/2*dt_U_p_x[k]);
			Y[i+1][k] = Y[i][k] + tau*(U_p_y[k]+tau/2*dt_U_p_y[k]);								
		}

//==================================================

    toc = clock();
    cpu_time[i] = ((double)toc - (double)tic) / (double)CLOCKS_PER_SEC;;
    sum += cpu_time[i];


	if(stop_step) break;
		}


  printf("The cost of CPU time for 2D equations of motion by Lagrangian method is %g seconds.\n", sum);

	if(!stop_step)
		printf("The maximum number of time steps is not enough for this calculation, t_all=%lf.\n",t_all);

//------------END OF THE MAIN LOOP-------------

  for(k = 0; k < NUM_CELL; ++k)
  {
    	free(l_pc[k]);
    	free(n_pc_x[k]);
    	free(n_pc_y[k]);
	free(dt_l_n_pc_x[k]);
	free(dt_l_n_pc_y[k]);
	free(dt_n_pc_x[k]);
	free(dt_n_pc_y[k]);
	free(posi_l_pc[k]);
	free(posi_n_pc_x[k]);
	free(posi_n_pc_y[k]);
	free(dt_posi_l_n_pc_x[k]);
	free(dt_posi_l_n_pc_y[k]);
	free(dt_posi_n_pc_x[k]);
	free(dt_posi_n_pc_y[k]);
	free(nega_l_pc[k]);
	free(nega_n_pc_x[k]);
	free(nega_n_pc_y[k]);
	free(dt_nega_l_n_pc_x[k]);
	free(dt_nega_l_n_pc_y[k]);
	free(dt_nega_n_pc_x[k]);
	free(dt_nega_n_pc_y[k]);
	free(M_pc_00[k]);
	free(M_pc_01[k]);
	free(M_pc_10[k]);
	free(M_pc_11[k]);
	free(dt_M_pc_00[k]);
	free(dt_M_pc_01[k]);
	free(dt_M_pc_10[k]);
	free(dt_M_pc_11[k]);
	free(F_pc_x[k]);
	free(F_pc_y[k]);
	free(dt_F_pc_x[k]);
	free(dt_F_pc_y[k]);
	free(P_c_lim[k]);
	free(U_c_lim[k]);
	free(V_c_lim[k]);
	free(CELL_CELL[k]);

	l_pc[k] = NULL;
	n_pc_x[k] = NULL;
	n_pc_y[k] = NULL;
	dt_l_n_pc_x[k] = NULL;
	dt_l_n_pc_y[k] = NULL;
	dt_n_pc_x[k] = NULL;
	dt_n_pc_y[k] = NULL;
	posi_l_pc[k] = NULL;
	posi_n_pc_x[k] = NULL;
	posi_n_pc_y[k] = NULL;
	dt_posi_l_n_pc_x[k] = NULL;
	dt_posi_l_n_pc_y[k] = NULL;
	dt_posi_n_pc_x[k] = NULL;
	dt_posi_n_pc_y[k] = NULL;
	nega_l_pc[k] = NULL;
	nega_n_pc_x[k] = NULL;
	nega_n_pc_y[k] = NULL;
	dt_nega_l_n_pc_x[k] = NULL;
	dt_nega_l_n_pc_y[k] = NULL;
	dt_nega_n_pc_x[k] = NULL;
	dt_nega_n_pc_y[k] = NULL;
	M_pc_00[k] = NULL;
	M_pc_01[k] = NULL;
	M_pc_10[k] = NULL;
	M_pc_11[k] = NULL;
	dt_M_pc_00[k] = NULL;
	dt_M_pc_01[k] = NULL;
	dt_M_pc_10[k] = NULL;
	dt_M_pc_11[k] = NULL;
	F_pc_x[k] = NULL;
	F_pc_y[k] = NULL;
	dt_F_pc_x[k] = NULL;
	dt_F_pc_y[k] = NULL;
	P_c_lim[k] = NULL;
	U_c_lim[k] = NULL;
	V_c_lim[k] = NULL;
	CELL_CELL[k] = NULL;
  }

  return i;
}
