#include <math.h>
#include <stdio.h>



void EUCCLHYD_iteration(double (* M_pc)[2], double gamma, double P, double RHO, double posi_l_pc, double posi_n_pc_x, double posi_n_pc_y, double nega_l_pc, double nega_n_pc_x, double nega_n_pc_y, double U, double V, double U_p_x, double U_p_y)
{
	double z_posi, z_nega;
	z_posi = RHO * (sqrt(gamma * P / RHO) + (gamma+1)/2 * fabs((U_p_x-U)*posi_n_pc_x+(U_p_y-V)*posi_n_pc_y));
	z_nega = RHO * (sqrt(gamma * P / RHO) + (gamma+1)/2 * fabs((U_p_x-U)*nega_n_pc_x+(U_p_y-V)*nega_n_pc_y));
	M_pc[0][0] = z_posi * posi_l_pc * posi_n_pc_x * posi_n_pc_x + z_nega * nega_l_pc * nega_n_pc_x * nega_n_pc_x;
	M_pc[0][1] = z_posi * posi_l_pc * posi_n_pc_x * posi_n_pc_y + z_nega * nega_l_pc * nega_n_pc_x * nega_n_pc_y;
	M_pc[1][0] = z_posi * posi_l_pc * posi_n_pc_y * posi_n_pc_x + z_nega * nega_l_pc * nega_n_pc_y * nega_n_pc_x;
	M_pc[1][1] = z_posi * posi_l_pc * posi_n_pc_y * posi_n_pc_y + z_nega * nega_l_pc * nega_n_pc_y * nega_n_pc_y;
}
