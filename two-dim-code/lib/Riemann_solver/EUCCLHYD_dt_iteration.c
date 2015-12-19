#include <math.h>
#include <stdio.h>



void EUCCLHYD_dt_iteration(double (* M_pc)[2], double gamma, double P, double RHO, double posi_l_pc, double posi_n_pc_x, double posi_n_pc_y, double nega_l_pc, double nega_n_pc_x, double nega_n_pc_y, double dt_posi_l_n_pc_x, double dt_posi_l_n_pc_y, double dt_posi_n_pc_x, double dt_posi_n_pc_y, double dt_nega_l_n_pc_x, double dt_nega_l_n_pc_y, double dt_nega_n_pc_x, double dt_nega_n_pc_y)
{
	double z;
	z = RHO * sqrt(gamma * P / RHO);
	M_pc[0][0] = z * (dt_posi_l_n_pc_x * posi_n_pc_x + posi_l_pc * posi_n_pc_x * dt_posi_n_pc_x + dt_nega_l_n_pc_x * nega_n_pc_x + nega_l_pc * nega_n_pc_x * dt_nega_n_pc_x);
	M_pc[0][1] = z * (dt_posi_l_n_pc_x * posi_n_pc_y + posi_l_pc * posi_n_pc_x * dt_posi_n_pc_y + dt_nega_l_n_pc_x * nega_n_pc_y + nega_l_pc * nega_n_pc_x * dt_nega_n_pc_y);
	M_pc[1][0] = z * (dt_posi_l_n_pc_y * posi_n_pc_x + posi_l_pc * posi_n_pc_y * dt_posi_n_pc_x + dt_nega_l_n_pc_y * nega_n_pc_x + nega_l_pc * nega_n_pc_y * dt_nega_n_pc_x);
	M_pc[1][1] = z * (dt_posi_l_n_pc_y * posi_n_pc_y + posi_l_pc * posi_n_pc_y * dt_posi_n_pc_y + dt_nega_l_n_pc_y * nega_n_pc_y + nega_l_pc * nega_n_pc_y * dt_nega_n_pc_y);
}
