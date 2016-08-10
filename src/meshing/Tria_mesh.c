#include <stdio.h>
#include <math.h>
#include <stdlib.h>


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
