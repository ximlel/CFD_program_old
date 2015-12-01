#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>


#include "../../../lib/Riemann_solver.h"
#include "../../../lib/custom.h"

#include <float.h>
#define min(x,y)  ( x<y?x:y )

/* This function use first order scheme to solve 2-D
 * equations of motion by Lagrangian method.
 *
 *config is the array of configuration data, the detail
 *         could be seen in the comments of the main function.
 *NUM_CELL is the number of the grids.
 */

//
int first_order_solver
(double * config, int NUM_CELL, int NUM_POINT, int NUM_BOUNDARY, int * CELL_POINT[], int * BOUNDARY_POINT[],
 double * RHO[], double * U[], double * V[], double * P[], double * X[], double * Y[],
 double * NORMAL_VELOCITY[], double * gamma, double * cpu_time, char * scheme, double CFL)
{
	int i, j, k, l; 

  clock_t tic, toc;
  double sum = 0.0;

int stop_step=0;

  int const N = (int)(config[5]);  // the number of time steps
  double const eps = config[4];    // the largest value could be
                                   // seen as zero
  double tau;
  tau = DBL_MAX;    // the length of the time step
  double t_all;
  t_all = 0.0;

  double a[NUM_CELL];   //acoustic
  double lambda[NUM_CELL]; //minimal distance

  double dt_volume_temp;

  int p_p,p_n;
  
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

  double * l_pc[NUM_CELL];
  initialize_memory(l_pc,NUM_CELL,CELL_POINT);
  double * n_pc_x[NUM_CELL];
  initialize_memory(n_pc_x,NUM_CELL,CELL_POINT);
  double * n_pc_y[NUM_CELL];
  initialize_memory(n_pc_y,NUM_CELL,CELL_POINT);

  double * posi_l_pc[NUM_CELL];
  initialize_memory(posi_l_pc,NUM_CELL,CELL_POINT);
  double * posi_n_pc_x[NUM_CELL];
  initialize_memory(posi_n_pc_x,NUM_CELL,CELL_POINT);
  double * posi_n_pc_y[NUM_CELL];
  initialize_memory(posi_n_pc_y,NUM_CELL,CELL_POINT);

  double * nega_l_pc[NUM_CELL];
  initialize_memory(nega_l_pc,NUM_CELL,CELL_POINT);
  double * nega_n_pc_x[NUM_CELL];
  initialize_memory(nega_n_pc_x,NUM_CELL,CELL_POINT);
  double * nega_n_pc_y[NUM_CELL];
  initialize_memory(nega_n_pc_y,NUM_CELL,CELL_POINT);

  double M_pc[2][2];

  double * M_pc_00[NUM_CELL];
  initialize_memory(M_pc_00,NUM_CELL,CELL_POINT);
  double * M_pc_01[NUM_CELL];
  initialize_memory(M_pc_01,NUM_CELL,CELL_POINT);
  double * M_pc_10[NUM_CELL];
  initialize_memory(M_pc_10,NUM_CELL,CELL_POINT);
  double * M_pc_11[NUM_CELL];
  initialize_memory(M_pc_11,NUM_CELL,CELL_POINT);

  double U_p_x[NUM_POINT];
  double U_p_y[NUM_POINT];

  int index[NUM_POINT];
  	for(l = 0; l < NUM_POINT; ++l)
		{
			index[l] = -1;
			for(k = 0; k < NUM_BOUNDARY; ++k)
				{
					if(BOUNDARY_POINT[0][k] == l) index[l] = k;
				}
		}

	double M_p[2][3];
	double temp_2[2];

	double A_p[3][4];
	double temp_3[3];

	double posi_l_p[NUM_BOUNDARY], posi_n_p_x[NUM_BOUNDARY], posi_n_p_y[NUM_BOUNDARY];
        double nega_l_p[NUM_BOUNDARY], nega_n_p_x[NUM_BOUNDARY], nega_n_p_y[NUM_BOUNDARY];

	double * F_pc_x[NUM_CELL];
	initialize_memory(F_pc_x,NUM_CELL,CELL_POINT);
	double * F_pc_y[NUM_CELL];
	initialize_memory(F_pc_y,NUM_CELL,CELL_POINT);


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

	
	double rho_i[NUM_CELL];
	double E[NUM_CELL];
 
	for(k = 0; k < NUM_CELL; ++k)
		{
			rho_i[k]=1.0/RHO[0][k];
			E[k]=P[0][k]/(gamma[k]-1.0)/RHO[0][k] + 0.5*(U[0][k]*U[0][k]+V[0][k]*V[0][k]);			 	 
		}



//------------THE MAIN LOOP-------------
	for(i = 0; i < N; ++i)
		{		
			tic = clock();		

			
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

				if(strcmp(scheme,"GLACE")==0)			  
				  GLACE(M_pc, gamma[k], P[i][k], RHO[i][k], l_pc[k][j], n_pc_x[k][j], n_pc_y[k][j]);
				else if(strcmp(scheme,"EUCCLHYD")==0)
				  EUCCLHYD(M_pc, gamma[k], P[i][k], RHO[i][k], posi_l_pc[k][j], posi_n_pc_x[k][j], posi_n_pc_y[k][j], nega_l_pc[k][j], nega_n_pc_x[k][j], nega_n_pc_y[k][j]);
                                else
                                 {
                                    printf("No Riemann solver!\n");
                                    break;
                                 }

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
			A_p[0][2]=0.0;
			A_p[1][2]=0.0;
			A_p[2][2]=0.0;

			for(NUM_PP = 0; NUM_PP < POINT_POINT[l][0]; ++NUM_PP)
					{
						k=POINT_POINT[l][NUM_PP*2+1];
						j=POINT_POINT[l][NUM_PP*2+2];
						U_p_x[l]=U_p_x[l] + l_pc[k][j]*n_pc_x[k][j]*P[i][k] + M_pc_00[k][j]*U[i][k] + M_pc_01[k][j]*V[i][k];
						U_p_y[l]=U_p_y[l] + l_pc[k][j]*n_pc_y[k][j]*P[i][k] + M_pc_10[k][j]*U[i][k] + M_pc_11[k][j]*V[i][k];
						M_p[0][0]=M_p[0][0] +  M_pc_00[k][j];
						M_p[0][1]=M_p[0][1] +  M_pc_01[k][j];
						M_p[1][0]=M_p[1][0] +  M_pc_10[k][j];
						M_p[1][1]=M_p[1][1] +  M_pc_11[k][j];
						A_p[0][2]=A_p[0][2] + l_pc[k][j]*n_pc_x[k][j];
						A_p[1][2]=A_p[1][2] + l_pc[k][j]*n_pc_y[k][j];

					}
								
			if(index[l]==-1)
				{
					M_p[0][2] = U_p_x[l];
					M_p[1][2] = U_p_y[l];					
					Gauss_elimination(2, M_p, temp_2);
					U_p_x[l] = temp_2[0];
					U_p_y[l] = temp_2[1];
				}
			else if(fabs(fabs(posi_n_p_x[index[l]]*nega_n_p_y[index[l]]) - fabs(posi_n_p_y[index[l]]*nega_n_p_x[index[l]])) > eps)
				{

				        M_p[0][0] = posi_n_p_x[index[l]];
					M_p[0][1] = posi_n_p_y[index[l]];
					M_p[1][0] = nega_n_p_x[index[l]];
					M_p[1][1] = nega_n_p_y[index[l]];
					M_p[0][2] = NORMAL_VELOCITY[0][index[l]];
					M_p[1][2] = NORMAL_VELOCITY[1][index[l]];
					Gauss_elimination(2, M_p, temp_2);
					U_p_x[l] = temp_2[0];
					U_p_y[l] = temp_2[1];			
				}
			else
				{

					A_p[0][0] = M_p[0][0];
					A_p[0][1] = M_p[0][1];
					A_p[1][0] = M_p[1][0];
					A_p[1][1] = M_p[1][1];
					A_p[2][0] = A_p[0][2];
					A_p[2][1] = A_p[1][2];
					A_p[0][3] = U_p_x[l];
					A_p[1][3] = U_p_y[l];
					A_p[2][3] = NORMAL_VELOCITY[0][index[l]]*posi_l_p[index[l]] + NORMAL_VELOCITY[1][index[l]]*nega_l_p[index[l]];
					Gauss_elimination(3, A_p, temp_3);
					U_p_x[l] = temp_3[0];
					U_p_y[l] = temp_3[1];
				}

		}

	for(k = 0; k < NUM_CELL; ++k)
		{
			for(j = 0; j < CELL_POINT[k][0]; ++j)
				{
					F_pc_x[k][j]=l_pc[k][j]*n_pc_x[k][j]*P[i][k] - M_pc_00[k][j]*(U_p_x[CELL_POINT[k][j+1]] - U[i][k]) - M_pc_01[k][j]*(U_p_y[CELL_POINT[k][j+1]] - V[i][k]);
					F_pc_y[k][j]=l_pc[k][j]*n_pc_y[k][j]*P[i][k] - M_pc_10[k][j]*(U_p_x[CELL_POINT[k][j+1]] - U[i][k]) - M_pc_11[k][j]*(U_p_y[CELL_POINT[k][j+1]] - V[i][k]);

				}
		}
	tau = 1.01*tau;
	for(k = 0; k < NUM_CELL; ++k)
		{
			a[k] = sqrt(gamma[k] * P[i][k] / RHO[i][k]);
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
				}
		dt_volume_temp = 0.0;
		for(j = 0; j < CELL_POINT[k][0]; ++j)
				{
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
					RHO[i+1][k] = RHO[i+1][k] + l_pc[k][j]*n_pc_x[k][j]*U_p_x[CELL_POINT[k][j+1]] +l_pc[k][j]*n_pc_y[k][j]*U_p_y[CELL_POINT[k][j+1]];
					     U[i+1][k] = U[i+1][k] - F_pc_x[k][j];
						 V[i+1][k] = V[i+1][k] - F_pc_y[k][j];
						 P[i+1][k] = P[i+1][k] - F_pc_x[k][j]*U_p_x[CELL_POINT[k][j+1]] - F_pc_y[k][j]*U_p_y[CELL_POINT[k][j+1]];
				}
			
			rho_i[k] = rho_i[k] + tau/MASS[k]*RHO[i+1][k];
			RHO[i+1][k] = 1.0/rho_i[k];
			     U[i+1][k] = U[i][k] + tau/MASS[k]*U[i+1][k];
			     V[i+1][k] = V[i][k] + tau/MASS[k]*V[i+1][k];
				 E[k] = E[k] + tau/MASS[k]*P[i+1][k];
				 P[i+1][k] = (E[k] - 0.5*(U[i+1][k]*U[i+1][k]+V[i+1][k]*V[i+1][k]))*(gamma[k]-1.0)*RHO[i+1][k];			 	 
		}
	

	for(k = 0; k < NUM_POINT; ++k)
		{				
			X[i+1][k] = X[i][k] + tau*U_p_x[k];
			Y[i+1][k] = Y[i][k] + tau*U_p_y[k];								
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
	free(posi_l_pc[k]);
	free(posi_n_pc_x[k]);
	free(posi_n_pc_y[k]);
	free(nega_l_pc[k]);
	free(nega_n_pc_x[k]);
	free(nega_n_pc_y[k]);
	free(M_pc_00[k]);
	free(M_pc_01[k]);
	free(M_pc_10[k]);
	free(M_pc_11[k]);
	free(F_pc_x[k]);
	free(F_pc_y[k]);
    l_pc[k] = NULL;
    n_pc_x[k] = NULL;
    n_pc_y[k] = NULL;
	posi_l_pc[k] = NULL;
	posi_n_pc_x[k] = NULL;
	posi_n_pc_y[k] = NULL;
	nega_l_pc[k] = NULL;
	nega_n_pc_x[k] = NULL;
	nega_n_pc_y[k] = NULL;
	M_pc_00[k] = NULL;
	M_pc_01[k] = NULL;
	M_pc_10[k] = NULL;
	M_pc_11[k] = NULL;
	F_pc_x[k] = NULL;
	F_pc_y[k] = NULL;
  }

  return i;  
}


