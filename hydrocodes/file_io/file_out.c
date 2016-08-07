#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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

extern double config[200];




/* This function write the solution into an output file.
 * It is quite simple so there will be no more comments.
 */
void _1D_file_write_TEC(int m, int N, double * RHO[], double * U[], double * P[], double * X, double * cpu_time, double * config, char * example, char * problem)
{
  FILE * fp_write;
  char file_data[100] = "";


//===================Write solution File=========================
  strcat(file_data, "../../../data_out/one-dim/");
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

 char rho_data[100] = "";
  strcpy(rho_data, file_data);
  strcat(rho_data, "/\0");
  strcat(rho_data, "RHO.tec\0");
 char u_data[100] = "";
  strcpy(u_data, file_data);
  strcat(u_data, "/\0");
  strcat(u_data, "U.tec\0");
 char p_data[100] = "";
  strcpy(p_data, file_data);
  strcat(p_data, "/\0");
  strcat(p_data, "P.tec\0");



  int j;


  if((fp_write = fopen(rho_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  fprintf(fp_write, "VARIABLES = \"X\",  \"RHO\"\n");
  fprintf(fp_write, "ZONE  I= %d , f=point\n",m);

    for(j = 0; j < m; ++j)
      fprintf(fp_write, "%.18lf\t %.18lf\n", 0.5 * (X[j] + X[j+1]), RHO[N][j]);
  fclose(fp_write);


  if((fp_write = fopen(u_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  fprintf(fp_write, "VARIABLES = \"X\",  \"U\"\n");
  fprintf(fp_write, "ZONE  I= %d , f=point\n",m);

    for(j = 0; j < m; ++j)
      fprintf(fp_write, "%.18lf\t %.18lf\n", 0.5 * (X[j] + X[j+1]), U[N][j]);
  fclose(fp_write);


  if((fp_write = fopen(p_data, "w")) == 0)
  {
    printf("Cannot open solution output file!\n");
    exit(1);
  }
  fprintf(fp_write, "VARIABLES = \"X\",  \"P\"\n");
  fprintf(fp_write, "ZONE  I= %d , f=point\n",m);

    for(j = 0; j < m; ++j)
      fprintf(fp_write, "%.18lf\t %.18lf\n", 0.5 * (X[j] + X[j+1]), P[N][j]);
  fclose(fp_write);

}
