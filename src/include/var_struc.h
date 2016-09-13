#ifndef VARSTRUC_H
#define VARSTRUC_H


#define EPS 0.0000000001


#define N_CONF 400

extern double config[];

#define CONF_INI(i,j) printf("%3d-th configuration = %g .\n", i, j)

#define CONF_ERR(n)														\
	do {																\
		fprintf(stderr, "Error in the %d-th value of the configuration!\n", n); \
		exit(2);														\
	} while (0)


struct flu_var {
	double *RHO, *U, *V, *W, *P, *PHI, *gamma;
};


struct mesh_var {
	int num_pt, num_ghost, *cell_type, **cell_pt, num_border[10], *border_pt, *border_cond;
	double *X, *Y, *Z;
	void (*bc)(struct flu_var * FV, int i, int **cell_pt, double t);
};

struct cell_var {
	int **cell_cell;
	double **n_x, **n_y, **n_z;
	double **F_rho, **F_e, **F_phi, **F_u, **F_v, **F_w;
	double *U_rho, *U_e, *U_phi, *U_u, *U_v, *U_w;
	double *X_c, *Y_c, *Z_c;
	double *vol;
	double *gradx_rho, *grady_rho, *gradz_rho;
	double *gradx_phi, *grady_phi, *gradz_phi;
	double *gradx_e, *grady_e, *gradz_e;
	double *gradx_u, *grady_u, *gradz_u;
	double *gradx_v, *grady_v, *gradz_v;
	double *gradx_w, *grady_w, *gradz_w;
};

	
#endif
