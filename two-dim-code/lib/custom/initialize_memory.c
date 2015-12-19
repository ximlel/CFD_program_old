#include <stdio.h>
#include <math.h>
#include <stdlib.h>


void initialize_memory(double * p[],int N,int * CELL_POINT[])
{
int k,i;
for(k = 0; k < N; ++k)
		{
			p[k] = (double *)malloc(CELL_POINT[k][0] * sizeof(double));
			if(p[k] == NULL)
				{
					for(i = 0; i < k; ++i)
						{
							free(p[i]);
							p[i] = NULL;
						}
					printf("Initialize_memory fail.\n");
					exit(5);
				}
		}
}

void initialize_memory_int(int * p[],int N,int * CELL_POINT[])
{
int k,i;
for(k = 0; k < N; ++k)
		{
			p[k] = (int *)malloc(CELL_POINT[k][0] * sizeof(int));
			if(p[k] == NULL)
				{
					for(i = 0; i < k; ++i)
						{
							free(p[i]);
							p[i] = NULL;
						}
					printf("Initialize_memory fail.\n");
					exit(5);
				}
		}
}
