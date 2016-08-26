#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../include/var_struc.h"
#include "../include/file_io.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))



#define PRINT_LINE(v)														\
	do {																\
		for(k = 0; k < num_cell; ++k)									\
			fprintf(fp, "\t%.15g", v[k]);								\
		fprintf(fp,"\n");												\
	} while (0)

void file_write_TEC(const struct flu_var FV, const struct mesh_var mv, const char * problem)
{
	int k;
	
	const int dim = (int)config[0], num_cell = (int)config[3];

	int cell_type = 0;
	for (k = 0; k < num_cell; k++)
		cell_type= MAX(mv.cell_pt[0][0], cell_type);
	
	FILE * fp;
  
	char file_data[FILENAME_MAX];	
	example_io(problem, file_data, 0);	

	strcat(file_data, "/FLU_VAR.tec");

	//===================Write solution File=========================
	
	if ((fp = fopen(file_data, "w")) == NULL)
		{
			fprintf(stderr, "Cannot open solution output file!\n");
			exit(1);
		}
  
	fprintf(fp, "TITLE = \"FE-Volume Brick Data\"\n");
	fprintf(fp, "VARIABLES = \"X\"");
	if (dim > 1)
		fprintf(fp, ", \"Y\"");
	if (dim > 2)
		fprintf(fp, ", \"Z\"");
	fprintf(fp, ", \"P\", \"RHO\", \"U\"");
	if (dim > 1)
		fprintf(fp, ", \"V\"");
	if (dim > 2)
		fprintf(fp, ", \"W\"");
	if (!isinf(config[2]))
		switch((int)config[2])
			{
			case 2 :
				fprintf(fp, ", \"PHI\"");
				break;				
			}
	fprintf(fp, "\n");
					
	fprintf(fp, "ZONE  NODES=%d , ELEMENTS=%d , DATAPACKING=BLOCK, ", mv.num_pt, num_cell);
	
	if (cell_type < 2)
		{
			printf("NON ZONETYPE!");
			fclose(fp);
			remove(file_data);
			exit(2);
		}
	else if (cell_type <= 2 && dim == 1)
		fprintf(fp, "ZONETYPE=FELINESEG\n");
	else if (cell_type <= 3 && dim == 2)
		fprintf(fp, "ZONETYPE=FETETRAHEDRON\n");
	else if (cell_type <= 4 && dim == 2)
		fprintf(fp, "ZONETYPE=FEQUADRILATERAL\n");
	else if (cell_type <= 4 && dim == 3)
		fprintf(fp, "ZONETYPE=FETETRAHEDRON\n");
	else if (cell_type <= 8 && dim == 3)
		fprintf(fp, "ZONETYPE=FEBRICK\n");
	else
		{
			printf("NON ZONETYPE!");
			fclose(fp);
			remove(file_data);
			exit(2);
		}	

	PRINT_LINE(mv.X);
	if (dim > 1)
		PRINT_LINE(mv.Y);
	if (dim > 2)
		PRINT_LINE(mv.Z);
	PRINT_LINE(FV.P);
	PRINT_LINE(FV.RHO);
	PRINT_LINE(FV.U);
	if (dim > 1)
		PRINT_LINE(FV.V);
	if (dim > 2)
		PRINT_LINE(FV.W);
	if (!isinf(config[2]))
		switch((int)config[2])
			{
			case 2 :
				PRINT_LINE(FV.PHI);
				break;				
			}									
	
	for(k = 0; k < num_cell; ++k)
		{
			if (mv.cell_pt[k][0] == cell_type)
				for(int i = 1; i <= cell_type; i++)
					{
						if (i <= mv.cell_pt[k][0])
							fprintf(fp, "\t%d", mv.cell_pt[k][i]+1);
						else
							fprintf(fp, "\t%d", mv.cell_pt[k][mv.cell_pt[k][0]]+1);
					}
			fprintf(fp, "\n");
		}

	fclose(fp);
}


