#include <math.h>
#include <stdio.h>



void EUCCLHYD(double (* M_pc)[2], double gamma, double P, double RHO, double posi_l_pc, double posi_n_pc_x, double posi_n_pc_y, double nega_l_pc, double nega_n_pc_x, double nega_n_pc_y)
{
	double z;
	z = RHO * sqrt(gamma * P / RHO);
	M_pc[0][0] = z * (posi_l_pc * posi_n_pc_x * posi_n_pc_x + nega_l_pc * nega_n_pc_x * nega_n_pc_x);
	M_pc[0][1] = z * (posi_l_pc * posi_n_pc_x * posi_n_pc_y + nega_l_pc * nega_n_pc_x * nega_n_pc_y);
	M_pc[1][0] = z * (posi_l_pc * posi_n_pc_y * posi_n_pc_x + nega_l_pc * nega_n_pc_y * nega_n_pc_x);
	M_pc[1][1] = z * (posi_l_pc * posi_n_pc_y * posi_n_pc_y + nega_l_pc * nega_n_pc_y * nega_n_pc_y);
}
