/*!\file file_in.c
 * \author Du Zhifang, Lei Xin
 * \brief This file is a collection of functions used to read files
 and check whether they are qualified. 
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <stdbool.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "../include/var_struc.h"
#include "../include/file_io.h"


#ifndef N_CONF
#define N_CONF 400 
#endif



static int file_pre_read(FILE * fp, char * add, _Bool rc)
{
	double n_x = config[13], n_y = config[14], n_z = config[15];

	
	int num = 0;
	
	/* "flag" helps us to count.
	 *
	 * - flase: when we read a number-using character (0, 1, 2, ..., e, E, minus sign and dot).
	 * -  true: when we read a non-number-using character.
	 */
	_Bool flag = false;

	int r_count, c_count;

	int ch;

	while ((ch = getc(fp)) != EOF)
		{		   
			if (rc && ch == '\n' && !isinf(n_x) && !isinf(n_y))
				{
					if (num%(int)n_x == 0)
						{
							if(num/(int)n_x == r_count || num/(int)n_x == (r_count+1))
								r_count = num/(int)config[13];
							else
								{
									fprintf(stderr, "Row test is failed in the file '%s'!\n", add);
									exit(2);
								}
						}
					else
						{
							fprintf(stderr,"Row test is failed in the file '%s'!\n", add);
							exit(2);
						}
				}
			
			if( isspace(ch)&&flag )
				{
					++num;
					flag = 0;
				}			
			else if((ch == 46||ch == 45||ch == 69||ch == 101||isdigit(ch)) && (!flag))
				flag = 1;
			else if(isspace(ch)&&(!flag))
				continue;
			else if((ch == 46||ch == 45||ch == 69||ch == 101||isdigit(ch)) && (flag) )
				continue;
			else
				{
					fprintf(stderr, "Input contains illegal character(ASCII=%d, num=%d)!\n", ch, num);
					exit(2);
				}
		}
	
	return num;
}


void initialize(char * add, double * F)
{
	FILE * fp;
	int num;  
	int file_read_state;


	// Open the initial data file.
	if((fp = fopen(add, "r")) == NULL)
		{
			fprintf(stderr, "Can't open initial data file: %s!\n", add);
			exit(1);
		}
	if(isinf(config[3]))
		config[3] = num;
	else if((num = file_pre_read(fp, add, 2)) != (int)config[3])
		{
			printf("Unequal! number of grid in '%s' is %d, but in config is %d.", add, num, (int)config[3]);
			exit(3);
		}
	fseek(fp, 0L, SEEK_SET); 
  

	F = (double *)malloc(num * sizeof(double));
	file_read_state = file_read(fp, F, num);
	fclose(fp);
	if(file_read_state)
		{
			if(file_read_state == num)
				printf("Error on file reading!\n");
			else
				printf("The %d-th entry in the file '%s' is not a number.\n", file_read_state, add);
			exit(2);
		}

}


/* This function read the configuration data file,
 * and store the configuration data in the array
 * "config".
 * config[0] is the polytropic index
 * config[1] is the total time
 * config[2] is the spatial grid size
 * config[3] is the largest value can be seen as zero
 * config[4] is the maximum number of time steps
 */
void configurate(char * add)
{
	double config_pre[N_CONF*2];		
	FILE * fp_data;

	//open the configuration data file.
	if((fp_data = fopen(add, "r")) == 0)
		{
			fprintf(stderr,"Cannot open configuration file!\n");
			exit(1);
		}

	int num, state;

	//read the configuration data file.
	if((num = file_pre_read(fp_data, add, 0)) == 0||num > N_CONF*2)
		{
			printf("Configuration file error.\n");
			exit(2);
		}
	else if (num%2)
		{
			printf("Number of the configuration file is not even.");
			exit(2);
		}

	fseek(fp_data, 0L, SEEK_SET);
	state = file_read(fp_data, config_pre, num);
	fclose(fp_data);


	int n;  
	if(state)
		{
			if(state == num)
				printf("Error on configuration data file reading!\n");
			else
				printf("The %d-th entry in the configuration file is not a number.\n", state);
			exit(2);
		}
	else
		{
			for(int i=0; i<num/2; i+=2)  
				{
					if((n = (int)config_pre[i])>200 || n<0)
						{
							printf("The configuration file index bounds.");
							exit(2);
						}
					else if(!isinf(config[n]))
						{
							printf("The configuration file index repeats.");
							exit(2);
						}
					else
						{
							config[n] = config_pre[i+1];
						}
				}
		}

	config_check;
}


int config_check(void)
{
	if(isinf(config[1])||isinf(config[5]))
		{
			printf("The total time and the maximum number of time steps must be set.");
			exit(2);
		}
	
	// The default largest value can be seen as zero is set as 1e-9.
	if(isinf(config[4]))
		config[4] = pow(0.1,9);
	else if(config[4] < 0)
		{
			printf("eps(%lf) should be positive.\n", config[4]);
			exit(2);
		}

	// The default largest value can be seen as zero is set as 1e-9.	
	if(isinf(config[0]))
		config[0] = 1.4;	
	else if(config[0] < (1.0 + config[4]))
		{
			printf("The constant of the perfect gas(%lf) should be larger than 1.0.\n", config[0]);
			exit(2);
		}

	printf("t_all         = %g\n", config[1]);
	printf("max_tim_step  = %d\n", (int)config[5]);
	printf("eps           = %g\n", config[4]);
	printf("gamma         = %g\n", config[0]);


	if((!isinf(config[13]))&&(!isinf(config[14])))
		{
			printf("n_x, n_y      =  %d, %d\n", (int)config[13], (int)config[14]);
		}
	
	return 0;
}
