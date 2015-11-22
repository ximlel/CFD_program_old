#include <stdio.h>
#include <math.h>
#include <stdlib.h>


int Blunt_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * NORMAL_VELOCITY[], double * gamma, double * config, int m, int n)
{
	int num_boundary;
	int i, k;
	num_boundary = 2*m+2*n;
	BOUNDARY_POINT[0] = (int *)malloc(num_boundary * sizeof(int));
	if(BOUNDARY_POINT[0] == NULL)
				{
					printf("NOT enough memory! BOUNDARY_POINT\n");
					exit(5);
				}
	NORMAL_VELOCITY[0] = (double *)malloc(num_boundary * sizeof(double));//posi
	NORMAL_VELOCITY[1] = (double *)malloc(num_boundary * sizeof(double));//nega
	if((NORMAL_VELOCITY[0] == NULL)|(NORMAL_VELOCITY[1] == NULL))
				{
					printf("NOT enough memory! NORMAL_VELOCITY\n");
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
			CELL_POINT[k][0] = 4;
			CELL_POINT[k][1] = k + k/n;
			CELL_POINT[k][2] = CELL_POINT[k][1] + 1;
			CELL_POINT[k][3] = k + k/n + n + 2;
			CELL_POINT[k][4] = CELL_POINT[k][3] - 1;
		}
    	// X, Y
	double hyp;
	hyp = (m/2*config[3])*(m/2*config[3]) / ((1.0+n/4*config[2])*(1.0+n/4*config[2]) - 1.0);

	for(i = 0; i <= m; i++)
		{
		for(k = 0; k <= n/2; ++k)
			{
		
				X[k+i*(n+1)] = k*config[2];
				Y[k+i*(n+1)] = i*config[3];
			}
				X[n+i*(n+1)] = n*config[2] - (1.0 + n/4*config[2] - sqrt(1.0+(i*config[3]-m/2*config[3])*(i*config[3]-m/2*config[3])/hyp));
				Y[n+i*(n+1)] = i*config[3];
		for(k = 1; k <= n/2; ++k)
			{
		
				X[k+n/2+i*(n+1)] = X[n/2+i*(n+1)] + (X[n+i*(n+1)] - X[n/2+i*(n+1)])/(n/2)*k;
				Y[k+n/2+i*(n+1)] = i*config[3];
			}
		}
	// BOUNDARY
		for(k = 0; k < n+1; ++k)	
		{
			BOUNDARY_POINT[0][k] = k; 
		}
	for(k = n+1; k < n+m+1; ++k)	
		{
			BOUNDARY_POINT[0][k] = BOUNDARY_POINT[0][k-1] + n + 1;
		}
	for(k = n+m+1; k < n*2 + m + 1; ++k)	
		{
			BOUNDARY_POINT[0][k] = BOUNDARY_POINT[0][k-1] - 1;
		}
	for(k = n*2 + m + 1; k < num_boundary; ++k)	
		{
			BOUNDARY_POINT[0][k] = BOUNDARY_POINT[0][k-1] - n - 1;
		}

	// NORMAL_VELOCIT
		for(k = 0; k < num_boundary; ++k)	
			{
				NORMAL_VELOCITY[0][k] = 0.0;
				NORMAL_VELOCITY[1][k] = 0.0;
			}
                for(k = n*2 + m + 1; k < num_boundary; ++k)	
		       {
			        NORMAL_VELOCITY[0][k] = -0.0118321595662000;
			        NORMAL_VELOCITY[1][k] = -0.0118321595662000;
		       }
                NORMAL_VELOCITY[0][n*2+m] = -0.0118321595662000;
                NORMAL_VELOCITY[1][0] = -0.0118321595662000;
		//gamma
		for(k = 0; k < m*n; ++k)
			{
				gamma[k] = config[0];
			}
		
		return 	num_boundary;
}
