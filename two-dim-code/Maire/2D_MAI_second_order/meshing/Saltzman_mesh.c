#include <stdio.h>
#include <math.h>
#include <stdlib.h>


int Saltzman_mesh
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
	for(k = 0; k < (m+1)*(n+1); ++k)
	{
		X[k] = (k%(n+1))*config[2];
		Y[k] = (k/(n+1))*config[3];
                X[k] = X[k] + (0.1 - Y[k]) * sin(M_PI*X[k]);
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
			        NORMAL_VELOCITY[0][k] = -1.0;
			        NORMAL_VELOCITY[1][k] = -1.0;
		       }
                NORMAL_VELOCITY[0][n*2+m] = -1.0;
                NORMAL_VELOCITY[1][0] = -1.0;

		//gamma
		for(k = 0; k < m*n; ++k)
			{
				gamma[k] = config[0];
			}
		
		return 	num_boundary;
}
