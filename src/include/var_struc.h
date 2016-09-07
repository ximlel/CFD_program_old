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
	void (*bc)(struct flu_var * FV, double t);
};

#endif
