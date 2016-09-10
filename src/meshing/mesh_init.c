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
					fclose(fp);
					exit(2);
				}
			else
				{									
					fclose(fp);
					printf("Mesh file(%s.msh) has been read!\n", mesh_name);
					return mv;
				}
		}

	if (strcmp(mesh_name,"Sod_mesh") == 0)
		Sod_mesh(&mv);
	else if (strcmp(mesh_name,"Shear_mesh") == 0)
		Shear_mesh(&mv);
	else if (strcmp(mesh_name,"freee_mesh") == 0)
		free_mesh(&mv);
	else if (strcmp(mesh_name,"RMI_mesh") == 0)
		RMI_mesh(&mv);
	else if (strcmp(mesh_name,"cylinder_mesh") == 0)
		cylinder_mesh(&mv);
	else if (strcmp(mesh_name,"odd_even_mesh") == 0)
		odd_even_mesh(&mv);
	else if (strcmp(mesh_name,"odd_even_EW_mesh") == 0)
		odd_even_EW_mesh(&mv);
	else if (strcmp(mesh_name,"odd_even_EW_upstream_mesh") == 0)
		odd_even_EW_upstream_mesh(&mv);
	else if (strcmp(mesh_name,"odd_even_all_mesh") == 0)
		odd_even_all_mesh(&mv);
	else if (strcmp(mesh_name,"free_1D_mesh") == 0)
		free_1D_mesh(&mv);
	else
		{
			fprintf(stderr, "No mesh setting!\n");
			exit(2);
		}
	
	return mv;
}

