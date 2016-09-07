#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "../include/var_struc.h"

struct mesh_var mesh_load(char *mesh_name)
{
	struct mesh_var mv = {0, 0, NULL, NULL,	1, NULL, NULL, NULL, NULL, NULL};
	
}




int msh_read(FILE * fp, struct mesh_var * mv)
{
	int s_now = -1, s_max = 0; //section index
	int num_of = 1, num_tag; // number of section or tags.
	int *idx_N = NULL; // order of NODE;
	int temp[30]; // store cell_point data
	int i, type;
	const int dim = (int)config[0];
	
	const char *section[] = {
		"MeshFormat", "PhysicalNames", "Nodes",
		"Elements" , "Periodic", "NodeData", "ElementData",
		"ElementNodeData", "InterpolationScheme" };	
	char one_line[1000]; // input one line
	char *headptr, *endptr; // PTR of read data
	char Endsection[25];

	
	while (fgets(one_line, sizeof(one_line), fp) != NULL)
		{
			for (headptr = one_line; isspace(headptr); headptr++) ;

			if (*headptr == '$')
				{
					for (endptr = (++headptr); !isspace(endptr); endptr++) ;

					*endptr = '\0';
					if (s_now >= 0)
						{
							strcpy(Endsection, "End");
							strcpy(Endsection, section[s_now]);
							if (strcmp(headptr, Endsection) != 0)
								{
									fprintf(stderr, "End of the section in .msh file doesn't match!");
									return 0;										
								}
							else if (num_of > 0)
								{
									fprintf(stderr, "The count is not complete in the section in .msh file!");
									return 0;
								}
							else
								s_now = -1;
								
							break;
						}
				
					for (int s = s_max ; s < 9; s++)
							if (strcmp(headptr, section[s]) == 0)
								{
									s_now = s;
									s_max = s;
									break;
								}
				}
			else if(strlen(headptr) <= 0)
				;
			else if (s_now == 0 && num_of == 1)
				{
					if (strtod(headptr, &headptr) - 2.2 > EPS)
						{
							fprintf(stderr, "Version-number isn't not equal to 2.2 in .msh file!");
							return 0;
						}
					if (strtol(headptr, &headptr, 10) != 0)
						{
							fprintf(stderr, "The .msh file isn't ASCII file format!");
							return 0;
						}
					if (strtol(headptr, NULL, 10) != 8)
						{
							fprintf(stderr, "Currently only data-size = sizeof(double) is supported in .msh file!");
							return 0;
						}
					num_of--;					
				}
			else if (s_now == 2)
				{
					if(num_of <= 0)
						{											
							num_of = strtol(headptr, NULL, 10);
							mv->num_pt = num_of;
							if (num_of > 0)
								{
									idx_N = calloc(num_of, sizeof(int));
									mv->X = malloc(num_of * sizeof(double));
									mv->Y = malloc(num_of * sizeof(double));
									mv->Z = malloc(num_of * sizeof(double));
								}
						}
					else if (idx_N != NULL)
						{							
							idx_N[--num_of] = strtol(headptr, &headptr, 10);
							mv->X[num_of] = strtod(headptr, &headptr);
							mv->Y[num_of] = strtod(headptr, &headptr);
							mv->Z[num_of] = strtod(headptr, NULL);
						}
				}
			else if (s_now == 3)
				{
					if(num_of <= 0)
						{											
							num_of = strtol(headptr, NULL, 10);							
							if (num_of <= 0)
								break;
							else if (num_of > (int)config[3])
								printf("There are %d ghost cell!", mv->num_ghost = num_of - (int)config[3]);
							else if (num_of < (int)config[3])
								{
									fprintf(stderr,"There are not enough cell in .msh file!");
									return 0;
								}
									
							mv->cell_type = malloc(num_of * sizeof(int)); 
							mv->cell_pt = malloc(num_of * sizeof(void *));
						}
					else
						{
							strtol(headptr, &headptr, 10);

							type = strtol(headptr, &headptr, 10);
							if (dim == 1 && type != 1)
								break;
							else if (dim == 2 && (type < 2 || type > 3))
								break;
							else if (dim == 3 && (type < 4 || type > 7))
								break;
							mv->cell_type[--num_of] = type;				

							num_tag = strtol(headptr, &headptr, 10);
							while(num_tag-- != 0)
									strtol(headptr, &headptr, 10);

							for(i = 1; (temp[i] = strtol(headptr, &headptr, 10)) > 0; )
								i++;

							mv->cell_pt[num_of] = malloc((i+1) * sizeof(int));
							mv->cell_pt[num_of][0] = i;
							for(int j = 1; j <= i; j++)
								for(int k = 0; k < mv->num_pt; k++)
									if(temp[j] == idx_N[k])
										mv->cell_pt[num_of][j] = k;							
						}
				}
				
		}

	if (ferror(fp))
		{
			fprintf(stderr, "Read error occurrs in .msh file!\n");
			return 0;
		}
	free(idx_N);
	idx_N = NULL;
	
	return 1;
}



int Sod_mesh
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


		for(k = 0; k < n; ++k)	
		{
			BOUNDARY_POINT[1][k] = -2; //reflecting boundary condition.
		}
	for(k = n; k < n+m; ++k)	
		{
			BOUNDARY_POINT[1][k] = -1; //initial boundary condition.
		}
	for(k = n+m; k < n*2 + m; ++k)	
		{
			BOUNDARY_POINT[1][k] = -2;
		}
	for(k = n*2 + m; k < num_boundary; ++k)	
		{
			BOUNDARY_POINT[1][k] = -1;
		}


		//gamma
		for(k = 0; k < m*n; ++k)
			{
                                gamma[k] = config[0];
			}
		
		return 	num_boundary;
}
