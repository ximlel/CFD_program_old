
//extern double * U0;

double str2num(const char * number);


/* This function examine whether a string
 * represents a real number.
 * Transform the string represents a
 * negtive number into a string represents
 * a positive one and return its' sign.
 * It returns 0 if the string do not
 * represents a real number.
 */
int format_string(char * str);


/* This function reads the initial data file
 * to generate the initial data.
 * It returns 0 if successfully read the file,
 * while returns the index of the wrong entry.
 */
int file_read(FILE * fp, double * U, int num);



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
static int file_pre_read(FILE * fp, char * file_add, _Bool rc);


/* This function reads the initial data file. The function 
 * initialize return a pointer pointing to the position of
 * a block of memory consisting (m+1) variables* of type
 * double. The value of first of these variables is m.
 * The following m variables are the initial value.
 */
void initialize(char * add, double * F);

/* This function read the configuration data file,
 * and store the configuration data in the array "config".
 * config[?] stand for what? (see  Configuration_instructions.pdf)
 */
void configurate(char * add);

int config_check();


void file_write_VTK(int NUM_POINT, double * X, double * Y, int NUM_CELL, int * CELL_POINT[], double * RHO, double * U, double * V, double * P, double * cpu_time, double * config, char * example, char * problem);

void file_write_TEC(int NUM_POINT, double * X, double * Y, int NUM_CELL, int * CELL_POINT[], double * RHO, double * U, double * V, double * P, double * cpu_time, double * config, char * example, char * problem);

void file_two_species_write_TEC(int NUM_POINT, double * X, double * Y, int NUM_CELL, int * CELL_POINT[], double * RHO, double * U, double * V, double * P, double * CC, double * cpu_time, double * config, char * example, char * problem);
