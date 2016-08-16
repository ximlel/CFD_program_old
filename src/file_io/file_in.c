/*!\file file_in.c
 * \author Du Zhifang, Lei Xin
 * \brief This file is a collection of functions used to read files
 and check whether they are qualified. 
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>

#include "../include/file_io.h"


extern double config[200];


/*!\brief This function counts how many numbers are there in the initial data file. 
 * \author Du Zhifang, Lei Xin
 * \param[in] fp The pointer of the file to read in.
 * \param[in] test_lc Whether there is test for the range of data in the initial data file:
 * - 0 or 1 : no test,
 * - 2 : 2-D test (CR separates row),
 * - 3 : 3-D test (CR separates row, blank line separates column in x-y plane).
 * \return The number of the numbers in the initial data file.  
 * \retval -1 The given number of column is not coincided with that in the data file.
 * \retval 
 */
int file_pre_read(FILE * fp, int test_rc)
{
	int num = 0;
	
	/* \brief It helps us to count.
	 *
	 * We read characters one by one from the data file.
	 * The value of "flag" is:
	 * - 1: when we read a number-using character (0, 1, 2, ..., e, E, minus sign and the dot).
	 * - 0: when we read a non-number-using character.
	 */
	int flag = 0;

	int row_count;
	int colun_count;

	char ch;

	while((ch = getc(fp)) != EOF)
		{		   
			if(test_rc&&(ch == '\n')&&(!isinf(config[13]))&&(!isinf(config[14])))
				{
					if(num%(int)config[13]==0)
						{
							if(num/(int)config[13]==row_count||num/(int)config[13]==row_count+1)
								row_count = num/(int)config[13];
							else
								{
									printf("Row test is failed.\n");
									exit(1);
								}
						}
					else
						{
							printf("Row test is failed.\n");
							exit(1);
						}
				}
			
			if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (flag))
				{
					++num;
					flag = 0;
				}			
			else if( ((ch == 46)||(ch == 45)||(ch == 69)||(ch == 101)||((ch > 47) && (ch < 58))) && (!flag) )
				flag = 1;
			else if(((ch == ' ') || (ch == '\t') || (ch == '\n')) && (!flag))
				continue;
			else if( ((ch == 46)||(ch == 45)||(ch == 69)||(ch == 101)||((ch > 47) && (ch < 58))) && (flag) )
				continue;
			else
				{
					printf("Input contains illigal character(ASCII=%d, num=%d)!\n", (int)ch, num);
					exit(1);
				}
		}
	
	return num;
}


/* This function reads the initial data file. The function 
 * initialize return a pointer pointing to the position of
 * a block of memory consisting (m+1) variables* of type
 * double. The value of first of these variables is m.
 * The following m variables are the initial value.
 */
void initialize(char * add, double * F)
{
	FILE * fp;
	int num;  
	char ch;
	int file_read_state;


	//open the initial data file
	if((fp = fopen(add, "r")) == 0)
		{
			printf("Cannot open initial data file: '%s'!\n",add);
			exit(1);
		}
	if(isinf(config[3]))
		config[3] = num;
	else if((num = file_pre_read(fp)) != (int)config[3])
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
	double config_pre[400];		
	FILE * fp_data;

	//open the configuration data file.
	//printf("%s\n", add);
	if((fp_data = fopen(add, "r")) == 0)
		{
			printf("Cannot open configuration file!\n");
			exit(1);
		}

	int n_conf, state;

	//read the configuration data file.
	if((n_conf = file_pre_read(fp_data)) == 0||n_conf > 400)
		{
			printf("Configuration file error.\n");
			exit(2);
		}
	else if (n_conf%2)
		{
			printf("Number of the configuration file is not even.");
			exit(2);
		}

	fseek(fp_data, 0L, SEEK_SET);
	state = file_read(fp_data, config_pre, n_conf);
	fclose(fp_data);


	int n;  
	if(state)
		{
			if(state == n_conf)
				printf("Error on configuration data file reading!\n");
			else
				printf("The %d-th entry in the configuration file is not a number.\n", state);
			exit(2);
		}
	else
		{
			for(int i=0; i<n_conf/2; i+=2)  
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
