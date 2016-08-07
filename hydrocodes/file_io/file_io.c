#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>


#include "../file_io.h"
#include "../tools.h"



extern double * RHO;
extern double * U;
extern double * P;
extern double * V;
extern double * PHI;

extern double * config;



/* this function counts how many numbers are there
 * the initial data file.
 */
int file_pre_read_line(FILE * fp)
{
	int line = 0, column = 0, M = 0;
	int flag = 0;  /* We need to know how many numbers are there in
					* the initial data. "flag" helps us to count.
					* We read characters one by one from the data
					* file. The value of "flag" is 1 when read a
					* number-using character (i.e. 0, 1, 2, and so
					* on and the dot), while is 0 when read a 
					* non-number-using character. 
					*/
	char ch;

	while((ch = getc(fp)) != EOF)
		{
			//if(((ch < 45) || ((ch > 57) && (ch != 69) && (ch != 101)) || (ch == 47)) && (flag))
			if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (flag))
				{
					++column;
					flag = 0;
					if(ch == '\n')
						{
							if(!line)
								M = column;
							else if(column != M)
								{
									printf("Input data error, line=%d, M=%d, column=%d\n", line, M, column);
									if(RHO)
										free(RHO);
									if(U)
										free(U);
									if(V)
										free(V);
									if(P)
										free(P);
									exit(1);
								}

							++line;
							column = 0;
						}
				}
			else if( ((ch == 46)||(ch == 45)||(ch == 69)||(ch == 101)||((ch > 47) && (ch < 58))) && (!flag) )
				flag = 1;
			else if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (!flag))
				{
					if(ch == '\n')
						{
							if(!line)
								M = column;
							else if(column != M)
								{
									printf("Input data error, line=%d, M=%d, column=%d\n", line, M, column);
									if(RHO0)
										free(RHO0);
									if(U0)
										free(U0);
									if(V0)
										free(V0);
									if(P0)
										free(P0);
									exit(1);
								}

							++line;
							column = 0;
						}
					continue;
				}
			else if( ((ch == 46)||(ch == 45)||(ch == 69)||(ch == 101)||((ch > 47) && (ch < 58))) && (flag) )
				continue;
			else
				{
					printf("Input contains illigal character(ASCII=%d, flag=%d)!\n", (int)ch, flag);
					if(RHO0)
						free(RHO0);
					if(U0)
						free(U0);
					if(V0)
						free(V0);
					if(P0)
						free(P0);
					exit(1);
				}
		}

	if(flag)
		++column;
	if(column)
		{
			if(column != M)
				{
					if(ch == EOF)
						printf("hhh\n");
					printf("Input data error, line=%d, M=%d, column=%d\n", line, M, column);
					if(RHO0)
						free(RHO0);
					if(U0)
						free(U0);
					if(V0)
						free(V0);
					if(P0)
						free(P0);
					exit(1);
				}
			++line;
			column = 0;
		}

	return line;
}



int file_pre_read_column(FILE * fp)
{
	int column = 0;
	int flag = 0;  /* We need to know how many numbers are there in
					* the initial data. "flag" helps us to count.
					* We read characters one by one from the data
					* file. The value of "flag" is 1 when read a
					* number-using character (i.e. 0, 1, 2, and so
					* on and the dot), while is 0 when read a 
					* non-number-using character. 
					*/
	char ch;

	while((ch = getc(fp)) != EOF)
		{
			//if(((ch < 45) || ((ch > 57) && (ch != 69) && (ch != 101)) || (ch == 47)) && (flag))
			if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (flag))
				{
					++column;
					flag = 0;
					if(ch == '\n')
						break;
				}
			else if( ((ch == 46)||(ch == 45)||(ch == 69)||(ch == 101)||((ch > 47) && (ch < 58))) && (!flag) )
				flag = 1;
			else if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (!flag))
				{
					if(ch == '\n')
						break;
					continue;
				}
			else if( ((ch == 46)||(ch == 45)||(ch == 69)||(ch == 101)||((ch > 47) && (ch < 58))) && (flag) )
				continue;
			else
				{
					printf("Input contains illigal character(ASCII=%d, flag=%d)!\n", (int)ch, flag);
					if(RHO0)
						free(RHO0);
					if(U0)
						free(U0);
					if(V0)
						free(V0);
					if(P0)
						free(P0);
					exit(1);
				}
		}

	if(flag)
		++column;

	return column;
}


/* This function reads the initial data file. The function 
 * initialize return a pointer pointing to the position of
 * a block of memory consisting (m+1) variables* of type
 * double. The value of first of these variables is m.
 * The following m variables are the initial value.
 */
void initialize(char * name, char * addrho, char * addu, char * addv, char * addp)
{
	FILE * fp_data;
	int num_rho = 0, line_rho, column_rho;  
	int num_u = 0, line_u, column_u;  
	int num_v = 0, line_v, column_v;  
	int num_p = 0, line_p, column_p;  
	char  ch;
	int file_read_state;
	//double * U0;


	//open the initial data file
	//printf("%s will open the initial data file: ", name);
	//printf("%s\n", add);
	if((fp_data = fopen(addrho, "r")) == 0)
		{
			printf("Cannot open initial data file rho!\n");
			exit(1);
		}
	//read the initial data file, RHO0 is the array of initial data
	line_rho = file_pre_read_line(fp_data);
	fseek(fp_data, 0L, SEEK_SET);
	column_rho = file_pre_read_column(fp_data);
	fseek(fp_data, 0L, SEEK_SET);
	num_rho = line_rho * column_rho;
	//read
	RHO0 = (double *)malloc((num_rho + 2) * sizeof(double));
	RHO0[0] = (double)line_rho;
	RHO0[1] = (double)column_rho;
	file_read_state = file_read(fp_data, RHO0+2, num_rho);
	fclose(fp_data);
	//check
	if(file_read_state)
		{
			free(RHO0);
			if(file_read_state == num_rho)
				printf("Error on file reading! rho\n");
			else
				printf("\nThe %dth entry in the file rho is not a number.\n", file_read_state);
			exit(2);
		}


	if((fp_data = fopen(addu, "r")) == 0)
		{
			printf("Cannot open initial data file u!\n");
			exit(1);
		}
	//read the initial data file, U0 is the array of initial data
	line_u = file_pre_read_line(fp_data);
	fseek(fp_data, 0L, SEEK_SET);
	column_u = file_pre_read_column(fp_data);
	fseek(fp_data, 0L, SEEK_SET);
	num_u = line_u * column_u;
	//check
	if(line_u != line_rho)
		{
			printf("Unequal! line_u=%d\tline_rho=%d\n", line_u, line_rho);
			free(RHO0);
			exit(5);
		}
	if(column_u != column_rho)
		{
			printf("Unequal! column_u=%d\tcolumn_rho=%d\n", column_u, column_rho);
			free(RHO0);
			exit(5);
		}
	//read
	U0 = (double *)malloc((num_u + 2) * sizeof(double));
	U0[0] = (double)line_u;
	U0[1] = (double)column_u;
	file_read_state = file_read(fp_data, U0+2, num_u);
	fclose(fp_data);
	//check
	if(file_read_state)
		{
			free(RHO0);
			free(U0);
			if(file_read_state == num_u)
				printf("Error on file reading! u\n");
			else
				printf("\nThe %dth entry in the file u is not a number.\n", file_read_state);
			exit(2);
		}



	if((fp_data = fopen(addv, "r")) == 0)
		{
			printf("Cannot open initial data file v!\n");
			exit(1);
		}
	//read the initial data file, U0 is the array of initial data
	line_v = file_pre_read_line(fp_data);
	fseek(fp_data, 0L, SEEK_SET);
	column_v = file_pre_read_column(fp_data);
	fseek(fp_data, 0L, SEEK_SET);
	num_v = line_v * column_v;
	//check
	if(line_v != line_rho)
		{
			printf("Unequal! line_u=%d\tline_rho=%d\n", line_v, line_rho);
			free(RHO0);
			free(U0);
			exit(5);
		}
	if(column_u != column_rho)
		{
			printf("Unequal! column_u=%d\tcolumn_rho=%d\n", column_v, column_rho);
			free(RHO0);
			free(U0);
			exit(5);
		}
	//read
	V0 = (double *)malloc((num_v + 2) * sizeof(double));
	V0[0] = (double)line_v;
	V0[1] = (double)column_v;
	file_read_state = file_read(fp_data, V0+2, num_v);
	fclose(fp_data);
	//check
	if(file_read_state)
		{
			free(RHO0);
			free(U0);
			free(V0);
			if(file_read_state == num_u)
				printf("Error on file reading! v\n");
			else
				printf("\nThe %dth entry in the file v is not a number.\n", file_read_state);
			exit(2);
		}





	if((fp_data = fopen(addp, "r")) == 0)
		{
			printf("Cannot open initial data file p!\n");
			exit(1);
		}
	//read the initial data file, U0 is the array of initial data
	line_p = file_pre_read_line(fp_data);
	fseek(fp_data, 0L, SEEK_SET);
	column_p = file_pre_read_column(fp_data);
	fseek(fp_data, 0L, SEEK_SET);
	num_p = line_p * column_p;
	//check
	if(line_p != line_rho)
		{
			printf("Unequal! line_p=%d\tline_rho=%d\n", line_p, line_rho);
			free(RHO0);
			free(U0);
			free(V0);
			exit(5);
		}
	if(column_p != column_rho)
		{
			printf("Unequal! column_p=%d\tcolumn_rho=%d\n", column_p, column_rho);
			free(RHO0);
			free(U0);
			free(V0);
			exit(5);
		}
	//read
	P0 = (double *)malloc((num_p + 2) * sizeof(double));
	P0[0] = (double)line_p;
	P0[1] = (double)column_p;
	file_read_state = file_read(fp_data, P0+2, num_p);
	fclose(fp_data);
	//check
	if(file_read_state)
		{
			free(RHO0);
			free(U0);
			free(V0);
			free(P0);
			if(file_read_state == num_u)
				printf("Error on file reading! p\n");
			else
				printf("\nThe %dth entry in the file p is not a number.\n", file_read_state);
			exit(2);
		}

	printf("%s initialized, line=%d, column=%d.\n", name, line_rho, column_rho);
}


void initialize_CC(char * name, char * addcc)
{
	FILE * fp_data;
	int num_cc = 0, line_cc, column_cc;    

	int file_read_state;

	if((fp_data = fopen(addcc, "r")) == 0)
		{
			printf("Cannot open initial data file cc!\n");
			exit(1);
		}
	//read the initial data file, CC0 is the array of initial data
	line_cc = file_pre_read_line(fp_data);
	fseek(fp_data, 0L, SEEK_SET);
	column_cc = file_pre_read_column(fp_data);
	fseek(fp_data, 0L, SEEK_SET);
	num_cc = line_cc * column_cc;
	//read
	CC0 = (double *)malloc((num_cc + 2) * sizeof(double));
	CC0[0] = (double)line_cc;
	CC0[1] = (double)column_cc;
	file_read_state = file_read(fp_data, CC0+2, num_cc);
	fclose(fp_data);
	//check
	if(file_read_state)
		{
			free(CC0);
			if(file_read_state == num_cc)
				printf("Error on file reading! cc\n");
			else
				printf("\nThe %dth entry in the file cc is not a number.\n", file_read_state);
			exit(2);
		}
}



/* This function read the configuration data file,
 * and store the configuration data in the array
 * "config".
 * config[0] is the polytropic index
 * config[1] is the total time
 * config[2] is the spatial grid size in x direction
 * config[3] is the spatial grid size in y direction
 * config[4] is the largest value can be seen as zero
 * config[5] is the maximum number of time steps
 */
void configurate(double * config, char * name, char * add)
{
	FILE * fp_data;

	//open the configuration data file
	//printf("%s will open the cinfiguration data file: ", name);
	//printf("%s\n\n", add);
	if((fp_data = fopen(add, "r")) == 0)
		{
			printf("Cannot open configuration data file!\n");
			free(RHO0);
			free(U0);
			free(V0);
			free(P0);
			exit(1);
		}

	int n_conf, state;

	//read the configuration data file
	if((n_conf = file_pre_read_column(fp_data)) != N_CONF)
		{
			printf("Configuration data file error, n_config=%d.\n", n_conf);
			free(RHO0);
			free(U0);
			free(V0);
			free(P0);
			exit(2);
		}

	fseek(fp_data, 0L, SEEK_SET);
	state = file_read(fp_data, config, n_conf);
	fclose(fp_data);

	if(state)
		{
			printf("Configuration data file error, at %d.\n", state);
			free(RHO0);
			free(U0);
			free(P0);
			exit(2);
		}
	if(config[0] < 1)
		{
			printf("The polytropic index(%lf) should be larger than 1.\n", config[0]);
			free(RHO0);
			free(U0);
			free(V0);
			free(P0);
			exit(2);
		}
	if(config[3] < 0)
		{
			printf("eps(%lf) should be positive.\n", config[3]);
			free(RHO0);
			free(U0);
			free(V0);
			free(P0);
			exit(2);
		}

	printf("%s configurated:\n", name);
	printf("gamma     = %g\n", config[0]);
	printf("t_all     = %g\n", config[1]);
	printf("h_x       = %g\n", config[2]);
	printf("h_y       = %g\n", config[3]);
	printf("eps       = %g\n", config[4]);
	printf("tim_step       = %d\n", (int)config[5]);
}



void file_write_VTK(int NUM_POINT, double * X, double * Y, int NUM_CELL, int * CELL_POINT[], double * RHO, double * U, double * V, double * P, double * cpu_time, double * config, char * example, char * problem)
{

	int i, j, k;
	int SIZE=0;
	for(k = 0; k < NUM_CELL; ++k)
		{
			SIZE=SIZE+CELL_POINT[k][0]+1;
		}

	FILE * fp_write;
	char file_data[200] = "";



	//===================Write solution File=========================

	strcat(file_data, "../../../data_out/two-dim/");
	strcat(file_data, problem);

	int stat_mkdir = 0;
	DIR * dir_test = NULL;
	strcat(file_data, example);
	strcat(file_data, "\0");

	dir_test = opendir(file_data);
	if(dir_test != NULL)
		printf("Output directory \"%s\" already exists.\n", file_data);
	else
		{
			stat_mkdir = mkdir(file_data, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			if(stat_mkdir)
				{
					printf("Output directory \"%s\" construction failed.\n", file_data);
					exit(9);
				}
			else
				printf("Output directory \"%s\" constructed.\n", file_data);
		}
	closedir(dir_test);

	char rho_data[200] = "";
	strcpy(rho_data, file_data);
	strcat(rho_data, "/\0");
	strcat(rho_data, "RHO.vtk\0");
	char uv_data[200] = "";
	strcpy(uv_data, file_data);
	strcat(uv_data, "/\0");
	strcat(uv_data, "UV.vtk\0");
	char p_data[200] = "";
	strcpy(p_data, file_data);
	strcat(p_data, "/\0");
	strcat(p_data, "P.vtk\0");



	if((fp_write = fopen(rho_data, "w")) == 0)
		{
			printf("Cannot open solution output file!\n");
			exit(1);
		}
	fprintf(fp_write, "# vtk DataFile Version 1.0\n");
	fprintf(fp_write, "RHO\n");
	fprintf(fp_write, "ASCII\n");
	fprintf(fp_write, "\n");
	fprintf(fp_write, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp_write, "POINTS %d double\n",NUM_POINT);
	for(k = 0; k < NUM_POINT; ++k)
		{
			fprintf(fp_write, "\t%.16lf", X[k]);
			fprintf(fp_write, "\t%.16lf", Y[k]);
			fprintf(fp_write, "\t0.0\n");
		}
    fprintf(fp_write, "\n");
    fprintf(fp_write, "CELLS %d %d\n",NUM_CELL ,SIZE);
	for(k = 0; k < NUM_CELL; ++k)
		{
			for(i = 0; i <= CELL_POINT[k][0]; ++i)
				{
					fprintf(fp_write, "\t%d", CELL_POINT[k][i]);
				}
			fprintf(fp_write, "\n");
		}
	fprintf(fp_write, "\n");

	fprintf(fp_write, "CELL_TYPES %d\n",NUM_CELL);
	for(k = 0; k < NUM_CELL; ++k)
		{
			fprintf(fp_write, "\t7");
			fprintf(fp_write, "\n");
		}
	fprintf(fp_write, "\n");

	fprintf(fp_write, "CELL_DATA %d\n",NUM_CELL);
	fprintf(fp_write, "SCALARS rho double\n");
	fprintf(fp_write, "LOOKUP_TABLE default\n");
	for(k = 0; k < NUM_CELL; ++k)
		{
			fprintf(fp_write, "\t%.16lf", RHO[k]);
			fprintf(fp_write, "\n");
		}

	fclose(fp_write);



	if((fp_write = fopen(uv_data, "w")) == 0)
		{
			printf("Cannot open solution output file!\n");
			exit(1);
		}
    fprintf(fp_write, "# vtk DataFile Version 1.0\n");
	fprintf(fp_write, "U\n");
	fprintf(fp_write, "ASCII\n");
	fprintf(fp_write, "\n");
	fprintf(fp_write, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp_write, "POINTS %d double\n",NUM_POINT);
	for(k = 0; k < NUM_POINT; ++k)
		{
			fprintf(fp_write, "\t%.16lf", X[k]);
			fprintf(fp_write, "\t%.16lf", Y[k]);
			fprintf(fp_write, "\t0.0\n");
		}
    fprintf(fp_write, "\n");
    fprintf(fp_write, "CELLS %d %d\n",NUM_CELL ,SIZE);
	for(k = 0; k < NUM_CELL; ++k)
		{
			for(i = 0; i <= CELL_POINT[k][0]; ++i)
				{
					fprintf(fp_write, "\t%d", CELL_POINT[k][i]);
				}
			fprintf(fp_write, "\n");
		}
	fprintf(fp_write, "\n");

	fprintf(fp_write, "CELL_TYPES %d\n",NUM_CELL);
	for(k = 0; k < NUM_CELL; ++k)
		{
			fprintf(fp_write, "\t7");
			fprintf(fp_write, "\n");
		}
	fprintf(fp_write, "\n");

	fprintf(fp_write, "CELL_DATA %d\n",NUM_CELL);
	fprintf(fp_write, "VECTORS UV float\n");
	for(k = 0; k < NUM_CELL; ++k)
		{
			fprintf(fp_write, "\t%.16lf", U[k]);
			fprintf(fp_write, "\t%.16lf", V[k]);
			fprintf(fp_write, "\t0.0\n");
		}
	fprintf(fp_write, "\n");
	fprintf(fp_write, "SCALARS size double\n");
	fprintf(fp_write, "LOOKUP_TABLE default\n");
	for(k = 0; k < NUM_CELL; ++k)
		{
			fprintf(fp_write, "\t%.16lf", sqrt(U[k]*U[k]+V[k]*V[k]));
			fprintf(fp_write, "\n");
		}
	fclose(fp_write);



	if((fp_write = fopen(p_data, "w")) == 0)
		{
			printf("Cannot open solution output file!\n");
			exit(1);
		}
	fprintf(fp_write, "# vtk DataFile Version 1.0\n");
	fprintf(fp_write, "P\n");
	fprintf(fp_write, "ASCII\n");
	fprintf(fp_write, "\n");
	fprintf(fp_write, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp_write, "POINTS %d double\n",NUM_POINT);
	for(k = 0; k < NUM_POINT; ++k)
		{
			fprintf(fp_write, "\t%.16lf", X[k]);
			fprintf(fp_write, "\t%.16lf", Y[k]);
			fprintf(fp_write, "\t0.0\n");
		}
    fprintf(fp_write, "\n");
    fprintf(fp_write, "CELLS %d %d\n",NUM_CELL ,SIZE);
	for(k = 0; k < NUM_CELL; ++k)
		{
			for(i = 0; i <= CELL_POINT[k][0]; ++i)
				{
					fprintf(fp_write, "\t%d", CELL_POINT[k][i]);
				}
			fprintf(fp_write, "\n");
		}
	fprintf(fp_write, "\n");

	fprintf(fp_write, "CELL_TYPES %d\n",NUM_CELL);
	for(k = 0; k < NUM_CELL; ++k)
		{
			fprintf(fp_write, "\t7");
			fprintf(fp_write, "\n");
		}
	fprintf(fp_write, "\n");

	fprintf(fp_write, "CELL_DATA %d\n",NUM_CELL);
	fprintf(fp_write, "SCALARS p double\n");
	fprintf(fp_write, "LOOKUP_TABLE default\n");
	for(k = 0; k < NUM_CELL; ++k)
		{
			fprintf(fp_write, "\t%.16lf", P[k]);
			fprintf(fp_write, "\n");
		}
	fclose(fp_write);
}



void file_write_TEC(int NUM_POINT, double * X, double * Y, int NUM_CELL, int * CELL_POINT[], double * RHO, double * U, double * V, double * P, double * cpu_time, double * config, char * example, char * problem)
{

	int i, j, k;
	int SIZE=0;
	for(k = 0; k < NUM_CELL; ++k)
		{
			SIZE=SIZE+CELL_POINT[k][0]+1;
		}

	FILE * fp_write;
	char file_data[200] = "";

	//===================Write solution File=========================

	strcat(file_data, "../../../data_out/two-dim/");
	strcat(file_data, problem);
	strcat(file_data, "/\0");

	int stat_mkdir = 0;
	DIR * dir_test = NULL;
	strcat(file_data, example);


	dir_test = opendir(file_data);
	if(dir_test != NULL)
		printf("Output directory \"%s\" already exists.\n", file_data);
	else
		{
			stat_mkdir = CreateDir(file_data);
			if(stat_mkdir)
				{
					printf("Output directory \"%s\" construction failed.\n", file_data);
					exit(9);
				}
			else
				printf("Output directory \"%s\" constructed.\n", file_data);
		}
	closedir(dir_test);

	char rho_data[200] = "";
	strcpy(rho_data, file_data);
	strcat(rho_data, "/RHO.tec\0");
	char uv_data[200] = "";
	strcpy(uv_data, file_data);
	strcat(uv_data, "/UV.tec\0");
	char p_data[200] = "";
	strcpy(p_data, file_data);
	strcat(p_data, "/P.tec\0");



	if((fp_write = fopen(rho_data, "w")) == 0)
		{
			printf("Cannot open solution output file!\n");
			exit(1);
		}
	fprintf(fp_write, "VARIABLES = \"X\", \"Y\", \"RHO\"\n");
	fprintf(fp_write, "ZONE  N= %d , E= %d , ",NUM_POINT,NUM_CELL);
	if(CELL_POINT[0][0]==4)
		fprintf(fp_write, "ZONETYPE=FEQUADRILATERAL\n");
	else if(CELL_POINT[0][0]==3)
		fprintf(fp_write, "ZONETYPE=FETRIANGLE\n");
	else
		{
			printf("NON ZONETYPE!");
			exit(1);
		}	
	fprintf(fp_write, "DATAPACKING=BLOCK\n");
	fprintf(fp_write, "VARLOCATION=([1-2]=NODAL, [3]=CELLCENTERED)\n");
	for(k = 0; k < NUM_POINT; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", X[k]);
		}
	for(k = 0; k < NUM_POINT; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", Y[k]);
		}
	for(k = 0; k < NUM_CELL; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", RHO[k]);
		}
	for(k = 0; k < NUM_CELL; ++k)
		{
			for(i = 1; i <= CELL_POINT[k][0]; ++i)
				{
					fprintf(fp_write, "\t%d", CELL_POINT[k][i]+1);
				}
			fprintf(fp_write, "\n");
		}
	fclose(fp_write);



	if((fp_write = fopen(uv_data, "w")) == 0)
		{
			printf("Cannot open solution output file!\n");
			exit(1);
		}
	fprintf(fp_write, "VARIABLES = \"X\", \"Y\", \"U\", \"V\"\n");
	fprintf(fp_write, "ZONE  N= %d , E= %d , ",NUM_POINT,NUM_CELL);
	if(CELL_POINT[0][0]==4)
		fprintf(fp_write, "ZONETYPE=FEQUADRILATERAL\n");
	else if(CELL_POINT[0][0]==3)
		fprintf(fp_write, "ZONETYPE=FETRIANGLE\n");
	else
		{
			printf("NON ZONETYPE!");
			exit(1);
		}	
	fprintf(fp_write, "DATAPACKING=BLOCK\n");
	fprintf(fp_write, "VARLOCATION=([1-2]=NODAL, [3-4]=CELLCENTERED)\n");
	for(k = 0; k < NUM_POINT; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", X[k]);
		}
	for(k = 0; k < NUM_POINT; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", Y[k]);
		}
	for(k = 0; k < NUM_CELL; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", U[k]);
		}
	for(k = 0; k < NUM_CELL; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", V[k]);
		}
	for(k = 0; k < NUM_CELL; ++k)
		{
			for(i = 1; i <= CELL_POINT[k][0]; ++i)
				{
					fprintf(fp_write, "\t%d", CELL_POINT[k][i]+1);
				}
			fprintf(fp_write, "\n");
		}
	fclose(fp_write);



	if((fp_write = fopen(p_data, "w")) == 0)
		{
			printf("Cannot open solution output file!\n");
			exit(1);
		}
	fprintf(fp_write, "VARIABLES = \"X\", \"Y\", \"P\"\n");
	fprintf(fp_write, "ZONE  N= %d , E= %d , ",NUM_POINT,NUM_CELL);
	if(CELL_POINT[0][0]==4)
		fprintf(fp_write, "ZONETYPE=FEQUADRILATERAL\n");
	else if(CELL_POINT[0][0]==3)
		fprintf(fp_write, "ZONETYPE=FETRIANGLE\n");
	else
		{
			printf("NON ZONETYPE!");
			exit(1);
		}	
	fprintf(fp_write, "DATAPACKING=BLOCK\n");
	fprintf(fp_write, "VARLOCATION=([1-2]=NODAL, [3]=CELLCENTERED)\n");
	for(k = 0; k < NUM_POINT; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", X[k]);
		}
	for(k = 0; k < NUM_POINT; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", Y[k]);
		}
	for(k = 0; k < NUM_CELL; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", P[k]);
		}
	for(k = 0; k < NUM_CELL; ++k)
		{
			for(i = 1; i <= CELL_POINT[k][0]; ++i)
				{
					fprintf(fp_write, "\t%d", CELL_POINT[k][i]+1);
				}
			fprintf(fp_write, "\n");
		}
	fclose(fp_write);
}



void file_two_species_write_TEC(int NUM_POINT, double * X, double * Y, int NUM_CELL, int * CELL_POINT[], double * RHO, double * U, double * V, double * P, double * CC, double * cpu_time, double * config, char * example, char * problem)
{

	int i, j, k;
	int SIZE=0;
	for(k = 0; k < NUM_CELL; ++k)
		{
			SIZE=SIZE+CELL_POINT[k][0]+1;
		}

	FILE * fp_write;
	char file_data[200] = "";

	//===================Write solution File=========================

	strcat(file_data, "../../../data_out/two-dim/\0");
	strcat(file_data, problem);
	strcat(file_data, "/\0");

	int stat_mkdir = 0;
	DIR * dir_test = NULL;
	strcat(file_data, example);


	dir_test = opendir(file_data);
	if(dir_test != NULL)
		printf("Output directory \"%s\" already exists.\n", file_data);
	else
		{
			stat_mkdir = CreateDir(file_data);
			if(stat_mkdir)
				{
					printf("Output directory \"%s\" construction failed.\n", file_data);
					exit(9);
				}
			else
				printf("Output directory \"%s\" constructed.\n", file_data);
		}
	closedir(dir_test);

	char rho_data[200] = "";
	strcpy(rho_data, file_data);
	strcat(rho_data, "/RHO.tec\0");
	char uv_data[200] = "";
	strcpy(uv_data, file_data);
	strcat(uv_data, "/UV.tec\0");
	char p_data[200] = "";
	strcpy(p_data, file_data);
	strcat(p_data, "/P.tec\0");
	char cc_data[200] = "";
	strcpy(cc_data, file_data);
	strcat(cc_data, "/CC.tec\0");



	if((fp_write = fopen(rho_data, "w")) == 0)
		{
			printf("Cannot open solution output file!\n");
			exit(1);
		}
	fprintf(fp_write, "VARIABLES = \"X\", \"Y\", \"RHO\"\n");
	fprintf(fp_write, "ZONE  N= %d , E= %d , ",NUM_POINT,NUM_CELL);
	if(CELL_POINT[0][0]==4)
		fprintf(fp_write, "ZONETYPE=FEQUADRILATERAL\n");
	else if(CELL_POINT[0][0]==3)
		fprintf(fp_write, "ZONETYPE=FETRIANGLE\n");
	else
		{
			printf("NON ZONETYPE!");
			exit(1);
		}	
	fprintf(fp_write, "DATAPACKING=BLOCK\n");
	fprintf(fp_write, "VARLOCATION=([1-2]=NODAL, [3]=CELLCENTERED)\n");
	for(k = 0; k < NUM_POINT; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", X[k]);
		}
	for(k = 0; k < NUM_POINT; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", Y[k]);
		}
	for(k = 0; k < NUM_CELL; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", RHO[k]);
		}
	for(k = 0; k < NUM_CELL; ++k)
		{
			for(i = 1; i <= CELL_POINT[k][0]; ++i)
				{
					fprintf(fp_write, "\t%d", CELL_POINT[k][i]+1);
				}
			fprintf(fp_write, "\n");
		}
	fclose(fp_write);



	if((fp_write = fopen(uv_data, "w")) == 0)
		{
			printf("Cannot open solution output file!\n");
			exit(1);
		}
	fprintf(fp_write, "VARIABLES = \"X\", \"Y\", \"U\", \"V\"\n");
	fprintf(fp_write, "ZONE  N= %d , E= %d , ",NUM_POINT,NUM_CELL);
	if(CELL_POINT[0][0]==4)
		fprintf(fp_write, "ZONETYPE=FEQUADRILATERAL\n");
	else if(CELL_POINT[0][0]==3)
		fprintf(fp_write, "ZONETYPE=FETRIANGLE\n");
	else
		{
			printf("NON ZONETYPE!");
			exit(1);
		}	
	fprintf(fp_write, "DATAPACKING=BLOCK\n");
	fprintf(fp_write, "VARLOCATION=([1-2]=NODAL, [3-4]=CELLCENTERED)\n");
	for(k = 0; k < NUM_POINT; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", X[k]);
		}
	for(k = 0; k < NUM_POINT; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", Y[k]);
		}
	for(k = 0; k < NUM_CELL; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", U[k]);
		}
	for(k = 0; k < NUM_CELL; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", V[k]);
		}
	for(k = 0; k < NUM_CELL; ++k)
		{
			for(i = 1; i <= CELL_POINT[k][0]; ++i)
				{
					fprintf(fp_write, "\t%d", CELL_POINT[k][i]+1);
				}
			fprintf(fp_write, "\n");
		}
	fclose(fp_write);



	if((fp_write = fopen(p_data, "w")) == 0)
		{
			printf("Cannot open solution output file!\n");
			exit(1);
		}
	fprintf(fp_write, "VARIABLES = \"X\", \"Y\", \"P\"\n");
	fprintf(fp_write, "ZONE  N= %d , E= %d , ",NUM_POINT,NUM_CELL);
	if(CELL_POINT[0][0]==4)
		fprintf(fp_write, "ZONETYPE=FEQUADRILATERAL\n");
	else if(CELL_POINT[0][0]==3)
		fprintf(fp_write, "ZONETYPE=FETRIANGLE\n");
	else
		{
			printf("NON ZONETYPE!");
			exit(1);
		}	
	fprintf(fp_write, "DATAPACKING=BLOCK\n");
	fprintf(fp_write, "VARLOCATION=([1-2]=NODAL, [3]=CELLCENTERED)\n");
	for(k = 0; k < NUM_POINT; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", X[k]);
		}
	for(k = 0; k < NUM_POINT; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", Y[k]);
		}
	for(k = 0; k < NUM_CELL; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", P[k]);
		}
	for(k = 0; k < NUM_CELL; ++k)
		{
			for(i = 1; i <= CELL_POINT[k][0]; ++i)
				{
					fprintf(fp_write, "\t%d", CELL_POINT[k][i]+1);
				}
			fprintf(fp_write, "\n");
		}
	fclose(fp_write);


	if((fp_write = fopen(cc_data, "w")) == 0)
		{
			printf("Cannot open solution output file!\n");
			exit(1);
		}
	fprintf(fp_write, "VARIABLES = \"X\", \"Y\", \"CC\"\n");
	fprintf(fp_write, "ZONE  N= %d , E= %d , ",NUM_POINT,NUM_CELL);
	if(CELL_POINT[0][0]==4)
		fprintf(fp_write, "ZONETYPE=FEQUADRILATERAL\n");
	else if(CELL_POINT[0][0]==3)
		fprintf(fp_write, "ZONETYPE=FETRIANGLE\n");
	else
		{
			printf("NON ZONETYPE!");
			exit(1);
		}	
	fprintf(fp_write, "DATAPACKING=BLOCK\n");
	fprintf(fp_write, "VARLOCATION=([1-2]=NODAL, [3]=CELLCENTERED)\n");
	for(k = 0; k < NUM_POINT; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", X[k]);
		}
	for(k = 0; k < NUM_POINT; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", Y[k]);
		}
	for(k = 0; k < NUM_CELL; ++k)
		{
			fprintf(fp_write, "\t%.16lf\n", CC[k]);
		}
	for(k = 0; k < NUM_CELL; ++k)
		{
			for(i = 1; i <= CELL_POINT[k][0]; ++i)
				{
					fprintf(fp_write, "\t%d", CELL_POINT[k][i]+1);
				}
			fprintf(fp_write, "\n");
		}
	fclose(fp_write);
}
