double Riemann_solver_exact(double * U_star, double * P_star, double gamma, double u_L, double u_R, double p_L, double p_R, double c_L, double c_R, int * CRW, double tol, int N);

void ROE_solver(double *F, double gamma, double P_L, double RHO_L, double U_L, double P_R, double RHO_R, double U_R, double *lambda_max, double delta);
