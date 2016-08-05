void ROE_solver(double *F, double gamma, double P_L, double RHO_L, double U_L, double V_L,double n_x, double n_y, double P_R, double RHO_R, double U_R,double V_R,double *lambda_max, double delta);

void HLL_solver(double *F, double gamma, double P_L, double RHO_L, double U_L, double V_L,double n_x, double n_y, double P_R, double RHO_R, double U_R,double V_R,double *lambda_max);

double Riemann_solver_exact(double * U_star, double * P_star, double gamma, double u_L, double u_R, double p_L, double p_R, double c_L, double c_R, int * CRW, double tol, int N);

void Roe_Goundov_solver(double *F, double gamma, double P_L, double RHO_L, double U_L, double P_R, double RHO_R, double U_R, double *lambda_max, double delta);

void Roe_HLL_solver(double *V_mk, double *F, double gamma, double P_L, double RHO_L, double U_L, double V_L, double P_R, double RHO_R, double U_R, double V_R, double *lambda_max, double delta);

void GLACE(double (* M_pc)[2], double gamma, double P, double RHO, double l_pc, double n_pc_x, double n_pc_y);

void EUCCLHYD(double (* M_pc)[2], double gamma, double P, double RHO, double posi_l_pc, double posi_n_pc_x, double posi_n_pc_y, double nega_l_pc, double nega_n_pc_x, double nega_n_pc_y);

void EUCCLHYD_dt(double (* M_pc)[2], double gamma, double P, double RHO, double posi_l_pc, double posi_n_pc_x, double posi_n_pc_y, double nega_l_pc, double nega_n_pc_x, double nega_n_pc_y, double dt_posi_l_n_pc_x, double dt_posi_l_n_pc_y, double dt_posi_n_pc_x, double dt_posi_n_pc_y, double dt_nega_l_n_pc_x, double dt_nega_l_n_pc_y, double dt_nega_n_pc_x, double dt_nega_n_pc_y);

void GLACE_dt(double (* M_pc)[2], double gamma, double P, double RHO, double l_pc, double n_pc_x, double n_pc_y, double dt_l_n_pc_x, double dt_l_n_pc_y, double dt_n_pc_x, double dt_n_pc_y);

void EUCCLHYD_iteration(double (* M_pc)[2], double gamma, double P, double RHO, double posi_l_pc, double posi_n_pc_x, double posi_n_pc_y, double nega_l_pc, double nega_n_pc_x, double nega_n_pc_y, double U, double V, double U_p_x, double U_p_y);

void EUCCLHYD_dt_iteration(double (* M_pc)[2], double gamma, double P, double RHO, double posi_l_pc, double posi_n_pc_x, double posi_n_pc_y, double nega_l_pc, double nega_n_pc_x, double nega_n_pc_y, double dt_posi_l_n_pc_x, double dt_posi_l_n_pc_y, double dt_posi_n_pc_x, double dt_posi_n_pc_y, double dt_nega_l_n_pc_x, double dt_nega_l_n_pc_y, double dt_nega_n_pc_x, double dt_nega_n_pc_y);
