void flu_conf_load(char *example);

/* This function reads the initial data file. The function 
 * initialize return a pointer pointing to the position of
 * a block of memory consisting (m+1) variables* of type
 * double. The value of first of these variables is m.
 * The following m variables are the initial value.
 */
void flu_var_init(char * add, double * F);

/* This function read the configuration data file,
 * and store the configuration data in the array "config".
 * config[?] stand for what? (see  Configuration_instructions.pdf)
 */
void configurate(char * add);

void config_check(void);


//void file_write_VTK(int NUM_POINT, double * X, double * Y, int NUM_CELL, int * CELL_POINT[], double * RHO, double * U, double * V, double * P, double * cpu_time, double * config, char * example, char * problem);

//void file_write_TEC(int NUM_POINT, double * X, double * Y, int NUM_CELL, int * CELL_POINT[], double * RHO, double * U, double * V, double * P, double * cpu_time, double * config, char * example, char * problem);

//void file_two_species_write_TEC(int NUM_POINT, double * X, double * Y, int NUM_CELL, int * CELL_POINT[], double * RHO, double * U, double * V, double * P, double * CC, double * cpu_time, double * config, char * example, char * problem);
