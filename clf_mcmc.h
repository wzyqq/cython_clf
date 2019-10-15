#ifndef __CLF_H__
#define __CLF_H__

#define column_per_file 9
#define boxsize 500
#define lf_file_number 6
#define xi_file_number 4
#define ndemension 11
#define chains 1
#define walks 1
#define mag_max 17
#define pi     3.141592653589793
#define distance_number 300001
#define wp_number 2
#define node_x5  1000
#define sqrt_2pi   2.506628275
#define r_num  431
#define mask_ratio 0.231

extern char Lf_mean[200];
extern char Inverse_matrix[200];
extern char Wp_mean[200];
extern char Inverse_matrix_xi[200];
extern double  Halo_file_redshift[lf_file_number];
extern double  Xi_redshift[xi_file_number];
extern double  Xi_redshift_name[xi_file_number];
extern double  Start_mag;
extern double  Start_mag_xi[xi_file_number];
extern double  End_mag;
extern int Mag_min;
extern int Mag_min_xi[xi_file_number];
extern int Mag_min_xi_acc[xi_file_number+1];
extern int Mag_min_xi_all;
extern int Start_point;
extern struct subhalo *Halodata[lf_file_number];
extern double Volume;
extern double Volume_inv;
extern double *Distribution;
extern double *Gauss_distribution;
extern double *Mass_bin;
extern double *Redshift_bin;
extern double *Redshift_bin_add;
extern double *Distance_bin_add;
extern int Mass_size, Z_size;
extern int Z_size_add;
extern int Z_size_cut[lf_file_number+1];
extern int Gauss_number;
extern double *N_h,*N_s;
extern double *Xicc, *Xics, *Xisc, *Xiss;
extern double *Mass_left, *Mass_right, *Mass_mid;
extern int Xi_nbin;
extern int Rp_bin;
extern int Rpi_bin;
extern int Nwp;
extern char **File_corr;
extern int Start_point_xi;
extern int Perpendicular_min;
extern double *Perpendicular;

void getlf(double *param, double *phi_analysis);
void getwp(double *param, double *w_analysis);
void getMCMC(double *init, int *argc_address, char ***argv_address);

struct subhalo
{
	int subhalo_id;
	double subhalo_mass, pos_x, pos_y ,pos_z, vel_x, vel_y, vel_z, acc_mass, redshift;
};
struct z_distance
{
	double redshift, distance;
};

struct z_distance Distance_data[distance_number];


#endif
