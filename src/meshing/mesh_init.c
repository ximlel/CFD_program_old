#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"
#include "../include/file_io.h"
#include "../include/meshing.h"




struct mesh_var mesh_load(const char *example, const char *mesh_name)
{
	struct mesh_var mv = {0, 0, NULL, NULL, {1}, NULL, NULL, NULL, NULL, NULL, NULL};

	char add_mkdir[FILENAME_MAX];
	example_io(example, add_mkdir, 1);
	char add[FILENAME_MAX];
	strcpy(add, add_mkdir);
	strcat(add, mesh_name);
	strcat(add, ".msh");

	FILE * fp;

	if ((fp = fopen(add, "r")) != NULL)
		{
			if(msh_read(fp, &mv) == 0)
				{
					free(mv.X);
					mv.X = NULL;
					free(mv.Y);
					mv.Y = NULL;
					free(mv.Z);
					mv.Z = NULL;
					free(mv.cell_type);
					mv.cell_type = NULL;
					free(mv.cell_pt);
					mv.cell_pt = NULL; // Dosen't free cell_pt[*], it is so difficult.
					// Dosen't free idx_N.
					fclose(fp);
					exit(2);
				}
			else
				{									
					fclose(fp);
					printf("Mesh file %s.msh has been read!\n", mesh_name);
					return mv;
				}
		}


	
	return mv;
}

