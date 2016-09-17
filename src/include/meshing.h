void Sod_mesh(struct mesh_var * mv);

void Shear_mesh(struct mesh_var * mv);

void free_mesh(struct mesh_var * mv);

void RMI_mesh(struct mesh_var * mv);

void cylinder_mesh(struct mesh_var * mv);

void odd_even_mesh(struct mesh_var * mv);

void odd_even_periodic_mesh(struct mesh_var * mv);

void odd_even_inflow_mesh(struct mesh_var * mv);

void rand_disturb_inflow_mesh(struct mesh_var * mv);

void free_1D_mesh(struct mesh_var * mv);


int msh_read(FILE * fp, struct mesh_var * mv);


struct mesh_var mesh_load(const char *example, const char *mesh_name);


void period_cell_modi(struct mesh_var * mv);

void period_ghost(struct cell_var * cv, struct mesh_var mv, double t);



/*
int Sod_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * gamma, double * config, int m, int n);
int Free_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * gamma, double * config, int m, int n);
int odd_even_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * gamma, double * config, int m, int n);
int odd_even_EW_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * gamma, double * config, int m, int n);
int odd_even_EW_upstream_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * gamma, double * config, int m, int n);
int Shear_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * gamma, double * config, int m, int n);
int Cylinder_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * gamma, double * config, int m, int n);
int RMI_mesh
(int * CELL_POINT[], double * X, double * Y, int * BOUNDARY_POINT[],
 double * gamma, double * config, int m, int n);
*/
