/*!\file file_in.c
 * \brief This file is a collection of functions used to read files
 and check whether they are qualified. 
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <dirent.h>

#include "../include/var_struc.h"
#include "../include/file_io.h"

#define EPS 0.0000000001


/*!\brief This function counts how many numbers are there in the initial data file. 
 * \param[in] fp The pointer of the file to read in.
 * \return The number of the numbers in the initial data file.
 */

#define STR_FLU_INI(sfv)							\
	do {											\
		strcpy(add, add_mkdir);						\
		strcat(add, "/" #sfv ".txt");				\
		flu_var_init(add, FV.sfv);					\
	} while(0)

void flu_conf_load(char *example)
{
	const int DIM = (int)config[0];
	
	// Find the directory of input data.
	char add_mkdir[FILENAME_MAX];
	switch(DIM)
		{
		case 1 :
			strcpy(add_mkdir, "../data_in/one-dim/");
			break;
		case 2 :
			strcpy(add_mkdir, "../data_in/two-dim/");
			break;
		default :
			fprintf(stderr, "Strange computational dimension!\n");
			break;
		}	
	DIR * dir_test = NULL;
	strcat(add_mkdir, example);

	dir_test = opendir(add_mkdir);
	if(dir_test == NULL)
		{
			fprintf(stderr, "Input directory is not exist!\n");
			exit(1);
		}
	closedir(dir_test);

	// We read the initial data file. 
	char add[FILENAME_MAX];
	strcpy(add, add_mkdir);
	strcat(add, "/config.txt");
	configurate(add);

	printf("%s is configurated, dimension = %d .\n", example, DIM);


			strcpy(add, add_mkdir);					
		strcat(add, "/RHO.txt");				
		flu_var_init(add, FV.RHO);

		//STR_FLU_INI(RHO);
	STR_FLU_INI(U);
	STR_FLU_INI(P);	
	if(DIM > 1)
		{		
			STR_FLU_INI(V);
		}	
	if(!isinf(config[2]))
		switch((int)config[2])
			{
			case 2 :			 							
				STR_FLU_INI(PHI);
				break;				
			}	
	printf("%s is initialized, grid number = %d .\n", example, (int)config[3]);
}

static int flu_var_count(FILE * fp, const char * add)
{	
	int num = 0; // Data number.
	
	// "flag" helps us to count.	
	// - flase: when we read a number-using character (0, 1, 2, ..., e, E, minus sign and dot).
	// -  true: when we read a non-number-using character.
	int flag = 0;

	int ch;


	while ((ch = getc(fp)) != EOF)
		{		
			// Count the data number.
			if (ch == 45 || ch == 46 || ch == 69 || ch == 101 || isdigit(ch))
				flag = 1;
			else if (!isspace(ch))
				{
					fprintf(stderr, "Input contains illegal character(ASCII=%d) in the file '%s'!\n", ch, add);
					exit(2);
				}
			else if (flag)
				{
					num++;
					flag = 0;
				}				
		}
	
	rewind(fp);
	
	return num;
}


#define RANGE_WRONG														\
	do {																\
		fprintf(stderr,"Data range is irregular in the file '%s'!\n", add); \
		exit(2);														\
	} while(0)

/*!\brief This function counts how many numbers are there in the initial data file. 
 * \param[in] fp The pointer of the file to read in.
 */
static void flu_var_read(FILE * fp, double * F, char * add)
{
	const int num_all = (int)config[3];
	const double n_x = config[13], n_y = config[14], n_z = config[15];

	static int D = 0;
	
	if (D == 0 && !isinf(n_x) && !isinf(n_y))
		{
			if (n_x < 1.0 || n_y < 1.0 || n_z < 1.0)
				{
					fprintf(stderr, "For the structural mesh, data number in some dimension < 1!\n");
					exit(2);
				}
			// Test whether the total number of data is matched.
			else if (num_all != (int)n_x * (int)n_y * (isinf(n_z) ? 1 : (int)n_z))
				{
					fprintf(stderr, "Data number is't matched to the structural mesh!\n");
					exit(2);
				}
			D = isinf(n_z) ? 2 : 3;
			if (D != (int)config[0])
				{
					fprintf(stderr, "Dimension is't matched to the structural mesh!\n");
					exit(2);
				}
		}
	

	int r_count = 0, c_count = 0;
	div_t r_div, c_div;

	
	int num = 0;

	int idx = 0;
	char flo_num[50]; // String to store floating numbers.
	char *endptr;
	
	int flag = 0;

	int ch;


	while ((ch = getc(fp)) != EOF)
		{			
			// Read the data number.
			if (ch == 45 || ch == 46 || ch == 69 || ch == 101 || isdigit(ch))
				{							
					flo_num[idx++] = (char)ch;
					flag = 1;
				}
			else if (!isspace(ch))
				{
					fprintf(stderr, "Input contains illegal character(ASCII=%d) in the file '%s'!\n", ch, add);
					exit(2);
				}
			else if (flag)
				{
					flo_num[idx] = '\0';
					idx = 0;
					if (num >= num_all)
						{
							num++;
							break;
						}											
					F[num++] = strtod(flo_num, &endptr);
					if (strlen(endptr) != 0)
						{
							fprintf(stderr,"Reading Sth. that isn't a floating number in the file '%s'!\n", add);
							exit(2);
						}
					flag = 0;
				}

			// Test whether the data range is regular and matched to the structual mesh.
			if (ch == '\n' && D > 0)
				{
					r_div = div(num, (int)n_x);
					c_div = div(num, (int)n_x * (int)n_y);
					if (r_div.rem != 0)
						RANGE_WRONG;
					else if (r_div.quot == (r_count+1))								
						r_count = r_div.quot;
					else if (r_div.quot != r_count)
						RANGE_WRONG;
					else if (D == 2)
						; 
					else if (c_div.rem == 0 && (c_div.quot == c_count || c_div.quot == (c_count+1)))
						c_count = c_div.quot;
					else						
						RANGE_WRONG;
				}		
		}
	
	// Test whether the total number of data is matched.
	if (num != num_all)
		{
			fprintf(stderr, "Data number isn't equal to the given total number in the file '%s'!\n", add);
			exit(2);
		}	
}


void flu_var_init(char * add, double * F)
{
	FILE * fp;

	// Open the initial data file.
	if ((fp = fopen(add, "r")) == NULL)
		{
			perror(add);
			exit(1);
		}
	
	if (isinf(config[3]))		
		config[3] = (double)flu_var_count(fp, add);
	else if (config[3] < 1.0)
		CONF_ERR(3);
	
	F = malloc((int)config[3] * sizeof(double));

	flu_var_read(fp, F, add);
	for (int i = 0; i<30; i++)
		printf("%d,%lf\n",i, F[i]);
	printf("ASAA:%s\n", add);
	fclose(fp);
}


static void config_read(FILE * fp)
{	
	char one_line[200]; // String to store one line.
	char *endptr;

	int i; // Index of config[*].
	
	double temp;

	while (fgets(one_line, sizeof(one_line), fp) != NULL)
		{
			// A line that doesn't begin with digits is a comment.
			i =strtol(one_line, &endptr, 10);			
			// If the value of config[i] doesn't exit, it is 0 by default.
			if (0 < i && i < N_CONF)
				{
					temp = strtod(endptr, NULL);
					if(fabs(config[i] - temp) > EPS)
						CONF_INI(i,temp);
					config[i] = temp;
				}
		}
}


void configurate(char * add)
{
	FILE * fp;

	//open the configuration data file.
	if((fp = fopen(add, "r")) == NULL)
		{
			perror(add);
			exit(1);
		}
	
	config_read(fp);
	fclose(fp);
}


void config_check(void)
{
	if(isinf(config[1]) && isinf(config[5]))
		{
			printf("The total time or the maximum number of time steps must be setted!\n");
			exit(2);
		}

	if(isinf(config[4]))
		config[4] = EPS;
	else if(config[4] < 0.0 || config[4] > 0.1)
		{
			printf("eps(%lf) should in (0, 0.1) !\n", config[4]);
			exit(2);
		}
	
	if(isinf(config[6]))
		config[6] = 1.4;	
	else if(config[6] < (1.0 + config[4]))
		{
			printf("The constant of the perfect gas(%lf) should be larger than 1.0 !\n", config[6]);
			exit(2);
		}	
}
