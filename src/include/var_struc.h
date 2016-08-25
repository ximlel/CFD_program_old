#ifndef VARSTRUC_H
#define VARSTRUC_H

extern double config[];

#define N_CONF 400

#define CONF_ERR(n)													\
	do {																\
		fprintf(stderr, "Error in the %d-th value of the configuration!\n", n); \
		exit(2);														\
	} while (0)

#define CONF_INI(i,j) printf("%3d-th configuration = %g .\n", i, j)


struct flu_var {
	double *RHO, *U, *V, *W, *P, *PHI, *gamma;
};

extern struct flu_var FV;


struct mesh_var {
	int **cell_point, *border[2];
	double *X, *Y, *Z;
};

#endif
