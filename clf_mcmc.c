#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_expint.h>
#include <mpi.h>
#include "clf_mcmc.h"

void getlf(double *param, double *phi_analysis)
{
	int i, j, k, num;
	double magnitude_abs_star0=0, m10=0, alpha0=0, beta0=0, gamma1=0, gamma2=0, gamma3=0, gamma4=0, gamma5=0, sigma_m=0, error=0;
	int temp_z=0;
	double magnitude_abs_star[lf_file_number]={0}, m1[lf_file_number]={0}, alpha[lf_file_number]={0}, beta[lf_file_number]={0};
	double posx=0;
	int  mag_loc=0;
	double *mag_abs;
	double *mag_abs_gauss;	
	double mag_app=0;


//	printf("good\n");
	magnitude_abs_star0 = *(param+0);
	m10 				= pow( 10.0,*(param+1) );	
	alpha0 			  	= *(param+2);	
	beta0 				= *(param+3);
	gamma1 				= *(param+4);
	gamma2 	  			= *(param+5);
	gamma3 				= *(param+6);
	gamma4 	  			= *(param+7);
	gamma5 				= *(param+8);
	sigma_m             = *(param+9);
	error 	  			= *(param+10);

	mag_abs  = 		(double *)malloc(Mass_size*Z_size*sizeof(double));
	mag_abs_gauss = (double *)malloc(Gauss_number*sizeof(double));


	for(i=0; i<Mass_size; i++)
	{
		for(j=0; j<Z_size; j++)
		{
			mag_abs[i*Z_size+j]=0;
		}
	}

	for(i=0; i<lf_file_number; i++)
	{
		magnitude_abs_star[i]	= magnitude_abs_star0 + gamma1*Halo_file_redshift[i];
//		printf("%lf\n", gamma1*Halo_file_redshift[i]);
//		printf("%lf\n", magnitude_abs_star[i]);
	
		m1[i] 					= pow( 10.0, log10(m10) + gamma2*Halo_file_redshift[i] );
//		printf("%lf\n", log10(m10) + gamma2*Halo_file_redshift[i]);
//		printf("%lf\n", m1[i]);

		alpha[i]				= alpha0 + gamma3*Halo_file_redshift[i];
//		printf("%lf\n", gamma3*Halo_file_redshift[i]);
//		printf("%lf\n", alpha[i]);

		beta[i]					= pow( 10.0, 
		GSL_MIN( log10(beta0) + gamma4*Halo_file_redshift[i]+gamma5*Halo_file_redshift[i]*Halo_file_redshift[i], 2.0) );
//		printf("%lf\n", log10(beta0) + gamma4*Halo_file_redshift[i]+gamma5*Halo_file_redshift[i]*Halo_file_redshift[i]);
//		printf("%lf\n", beta[i]);

		
		for(j=0; j<Mass_size; j++)
		{
			for(k=Z_size_cut[i]; k<Z_size_cut[i+1]; k++)	
			{
//				printf("%le\t%lf\n", Mass_bin[j], Redshift_bin[k]);
				mag_abs[j*Z_size+k] = magnitude_abs_star[i] - (alpha[i]+beta[i])*
				log10(Mass_bin[j]/m1[i]) + beta[i]*log10(1+(Mass_bin[j]/m1[i]));	
//				printf("k=%d\n", k);
//				printf("mag_abs[%d,%d]=%lf\n", j,k, mag_abs[j*Z_size+k]);			
				temp_z  = (int)(( Redshift_bin[k] )*100000 + 0.5);
//				printf("temp_z1=%d\n", temp_z1);
				if(temp_z>=138832)
				{
					break;
				}
				
				posx    = Distance_data[ temp_z ].distance;
//				printf("posx=%lf\n", posx);		
				for(num=0; num<Gauss_number; num++)
				{
					mag_abs_gauss[num] = mag_abs[j*Z_size+k]+8.0/Gauss_number*sigma_m*(num+0.5-Gauss_number/2);
//					printf("mag_abs=%lf\n", mag_abs);
//					printf("mag_abs_gauss[%d]=%lf\n", num, mag_abs_gauss[num]);

					mag_app = mag_abs_gauss[num] + 5*log10( posx*(1+Redshift_bin[k]) ) + 25;
			
//					printf("mag_abs[%d,%d]=%lf\n", j,k, mag_abs[j*Z_size+k]);			
//					if(num==0)
//					{
//						printf("mag_app=%lf\n", mag_app);	
//					}
					if( (mag_app>Start_mag)&&(mag_app<End_mag) )
					{
/*						printf("num=%d\n", num);
						printf("temp_z=%d\n", temp_z);							
						printf("posx=%lf\n", posx);
						printf("mag_abs[%d,%d]=%lf\n", j,k, mag_abs[j*Z_size+k]);			
						printf("mag_app=%lf\n", mag_app);	
*/						mag_loc = (int)(mag_app-Start_mag);
/*						if( (num==0) )
						{
							printf("num=%d\n", num);
							printf("log10(Mass_bin[j])=%lf\tRedshift_bin[k]=%lf\n", log10(Mass_bin[j]), Redshift_bin[k]);
							printf("mag_abs[%d,%d]=%lf\n", j, k, mag_abs[j*Z_size+k]);
							printf("mag_app=%lf\n", mag_app);	
//							printf("pos_loc=%d\n", pos_loc);		
//							printf("Gauss_distribution[%d]=%lf\n", num, Gauss_distribution[num]);																				
						}*/
//						if(num==0)
//						{
//							printf("\n");
//							printf("mag_loc=%d\n", mag_loc);	
//							printf("Gauss_distribution[%d]=%lf\n", num, Gauss_distribution[num]);
//							printf("mask_ratio%lf\n", mask_ratio);
//							printf("Distribution[%d,%d]=%lf\n", j, k, Distribution[j*Z_size+k]);
//							printf("\n");
//							printf("%d\t%lf\n", k,Distance_bin_add[k]);
//						}

						phi_analysis[ mag_loc ] = phi_analysis[ mag_loc ] + 
						Gauss_distribution[num]*Volume_inv*Distribution[j*Z_size+k]; 
//						printf("phi_analysis[%d]=%le\n", Mag_min_acc[pos_loc1]+mag_loc1, phi_analysis[ Mag_min_acc[pos_loc1]+mag_loc1 ]);
					}									
				}

			}
		}
	}

/*	for(i=0; i<Mag_min_all; i++)
	{
		printf("%lf\n", phi_analysis[i]*Volume);
		printf("%le\n", phi_analysis[i]);
	}		
*/
	free(mag_abs);
	free(mag_abs_gauss);	
}

void getwp(double *param, double *w_analysis)
{
	int i, j, k, num, p, q, s, m;
	double magnitude_abs_star0=0, m10=0, alpha0=0, beta0=0, gamma1=0, gamma2=0, gamma3=0, gamma4=0, gamma5=0, sigma_m=0, error=0;
	double magnitude_abs_star[xi_file_number]={0}, m1[xi_file_number]={0}, alpha[xi_file_number]={0}, beta[xi_file_number]={0};
	double mag_abs;
	double *mag_abs_gauss;
	int *temp_mag;
	double *temp_mag_f;	
	double *p_cen;
	double *xi_cc, *xi_cs, *xi_sc, *xi_ss;
	double *xi_gg;
	double *w1, *w2, *w;
	double high_limit;			
	double *n_halo;
	double xxi5[node_x5]={0};	
	double xwi5[node_x5]={0};	
	int lwp = 0;
	double *temp_error_distance_high, *temp_error_distance_low, *error_distance_high, *error_distance_low, *error_distance;
	int *loc_error_distance_high, *loc_error_distance_low;
	gsl_integration_glfixed_table *table_x5;

	magnitude_abs_star0 = *(param+0);
	m10 				= pow( 10.0,*(param+1) );	
	alpha0 			  	= *(param+2);	
	beta0 				= *(param+3);
	gamma1 				= *(param+4);
	gamma2 	  			= *(param+5);
	gamma3 				= *(param+6);
	gamma4 	  			= *(param+7);
	gamma5 				= *(param+8);
	sigma_m             = *(param+9);
	error 	  			= *(param+10);

	p_cen = (double *)malloc(Mag_min_xi_all*Xi_nbin*sizeof(double));
	xi_cc = (double *)malloc(Mag_min_xi_all*Rp_bin*Rpi_bin*sizeof(double));
	xi_cs = (double *)malloc(Mag_min_xi_all*Rp_bin*Rpi_bin*sizeof(double));
	xi_sc = (double *)malloc(Mag_min_xi_all*Rp_bin*Rpi_bin*sizeof(double));
	xi_ss = (double *)malloc(Mag_min_xi_all*Rp_bin*Rpi_bin*sizeof(double));
	xi_gg = (double *)malloc(Mag_min_xi_all*Rp_bin*Rpi_bin*sizeof(double));
	n_halo= (double *)malloc(Mag_min_xi_all*sizeof(double));
	w     = (double *)malloc(Mag_min_xi_all*wp_number*Rpi_bin*sizeof(double));
	w1    = (double *)malloc(Mag_min_xi_all*wp_number*Rpi_bin*sizeof(double));
	w2    = (double *)malloc(Mag_min_xi_all*wp_number*Rpi_bin*sizeof(double));
	mag_abs_gauss = (double *)malloc(Gauss_number*sizeof(double));
	temp_mag = (int *)malloc(Gauss_number*sizeof(int));
	temp_mag_f = (double *)malloc(Gauss_number*sizeof(double));


	temp_error_distance_high = (double *)malloc(Mag_min_xi_all*sizeof(double));
	temp_error_distance_low  = (double *)malloc(Mag_min_xi_all*sizeof(double));
	error_distance_high      = (double *)malloc(Mag_min_xi_all*sizeof(double));
	error_distance_low       = (double *)malloc(Mag_min_xi_all*sizeof(double));
	error_distance           = (double *)malloc(Mag_min_xi_all*sizeof(double));
	loc_error_distance_high  = (int *)malloc(Mag_min_xi_all*sizeof(int));
	loc_error_distance_low   = (int *)malloc(Mag_min_xi_all*sizeof(int));



	for(i=0; i<Mag_min_xi_all*Rp_bin*Rpi_bin; i++)
	{
		xi_cc[i] = xi_cs[i] = xi_sc[i] = xi_ss[i] = xi_gg[i] = 0;
	}
	for(i=0; i<Mag_min_xi_all*wp_number*Rpi_bin; i++)
	{
		w[i] = 0;
	}
	for(i=0; i<Mag_min_xi_all*Xi_nbin; i++)
	{
		p_cen[i] = 0;
	}
	for(i=0; i<Mag_min_xi_all; i++)
	{
		n_halo[i] = 0;
	}

	for(i=0; i<xi_file_number; i++)
	{
		for(j=Mag_min_xi_acc[i]; j<Mag_min_xi_acc[i+1]; j++)
		{
			temp_error_distance_high[j] = Xi_redshift[i]+error*(1+Xi_redshift[i])*0.5;
			temp_error_distance_low[j]  = Xi_redshift[i]-error*(1+Xi_redshift[i])*0.5;
			loc_error_distance_high[j]  = (int)(( temp_error_distance_high[j] )*100000 + 0.5);
			loc_error_distance_low[j]   = (int)(( temp_error_distance_low[j]  )*100000 + 0.5);
			error_distance_high[j]      = Distance_data[ loc_error_distance_high[j] ].distance;
			error_distance_low[j]       = Distance_data[ loc_error_distance_low[j]  ].distance;
			error_distance[j]           = error_distance_high[j] - error_distance_low[j];
//			printf("j=%d\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\n", j, temp_error_distance_high[j], temp_error_distance_low[j], loc_error_distance_high[j], 
//				loc_error_distance_low[j], error_distance_high[j], error_distance_low[j], error_distance[j]);
		}

		magnitude_abs_star[i]	= magnitude_abs_star0 + gamma1*Xi_redshift[i];
//		printf("%lf\n", gamma1*Halo_file_redshift[i]);
//		printf("%lf\n", magnitude_abs_star[i]);
	
		m1[i] 					= pow( 10.0, log10(m10) + gamma2*Xi_redshift[i] );
//		printf("%lf\n", log10(m10) + gamma2*Halo_file_redshift[i]);
//		printf("%lf\n", m1[i]);

		alpha[i]				= alpha0 + gamma3*Xi_redshift[i];
//		printf("%lf\n", gamma3*Halo_file_redshift[i]);
//		printf("%lf\n", alpha[i]);

		beta[i]					= pow( 10.0, 
		GSL_MIN( log10(beta0) + gamma4*Xi_redshift[i]+gamma5*Xi_redshift[i]*Xi_redshift[i], 2.0) );
//		printf("%lf\n", log10(beta0) + gamma4*Halo_file_redshift[i]+gamma5*Halo_file_redshift[i]*Halo_file_redshift[i]);
//		printf("%lf\n", beta[i]);

//		printf("Mag_min_xi[%d]=%d\n", i, Mag_min_xi[i]);
		for(j=0; j<Xi_nbin; j++)
		{
			mag_abs = magnitude_abs_star[i] - (alpha[i]+beta[i])*
			log10(Mass_mid[j]/m1[i]) + beta[i]*log10(1+(Mass_mid[j]/m1[i]));
			
			for(num=0; num<Gauss_number; num++)
			{
				mag_abs_gauss[num] = mag_abs+8.0/Gauss_number*sigma_m*(num+0.5-Gauss_number/2);
      			temp_mag_f[num] = (mag_abs_gauss[num] - Start_mag_xi[i]);
//				printf("mag_abs=%lf\n", mag_abs);
//				printf("mag_abs_gauss[%d]=%lf\n", num, mag_abs_gauss[num]);
//				printf("i=%d,j=%d,temp_mag_f[%d]=%lf\n", i, j, num, temp_mag_f[num]);
				if( (temp_mag_f[num]>=0) && (temp_mag_f[num]<Mag_min_xi[i]) )
				{
					temp_mag[num] = (int)temp_mag_f[num];
//					printf("\n");
//					printf("temp_mag[%d]=%d\n", num, temp_mag[num]);
//					printf("gauss=%lf\n", Gauss_distribution[num]);
					p_cen[ (Mag_min_xi_acc[i]+temp_mag[num])*Xi_nbin+j ] = p_cen[ (Mag_min_xi_acc[i]+temp_mag[num])*Xi_nbin+j ] + Gauss_distribution[num];
//					printf("p_cen[%d]=%lf\n", (Mag_min_xi_acc[i]+temp_mag[num])*Xi_nbin+j, p_cen[ (Mag_min_xi_acc[i]+temp_mag[num])*Xi_nbin+j ]);	
				}					
			}
		}
	}


	for(k=0; k<Mag_min_xi_all; k++)
	{
		for(j=0; j<Xi_nbin; j++)				
		{
//			p_cen[k*Xi_nbin+j] = 1;
//			n_halo[k] = n_halo[k] + (N_h[k*Xi_nbin+j])*p_cen[k*Xi_nbin+j];			
			n_halo[k] = n_halo[k] + (N_h[k*Xi_nbin+j]+N_s[k*Xi_nbin+j])*p_cen[k*Xi_nbin+j];
//			printf("p_cen[%d] = %lf\n", k*Xi_nbin+j, p_cen[k*Xi_nbin+j]);
		}
//		n_halo[k] = 0.009122;
//		printf("n_halo[%d]=%le\n", k, n_halo[k]);
	}
	
		
	for(i=0; i<Mag_min_xi_all; i++)
	{
		if(n_halo[i]>pow(10.0,-9.0))		
		{
			for(j=0; j<Rp_bin; j++)
			{
				for(k=0; k<Rpi_bin; k++)
				{
					lwp = 0;

					for(p=0; p<Xi_nbin; p++)
					{
						for(q=p; q<Xi_nbin; q++)
						{
//							if(p>98)
//							{
								xi_cc[i*Rp_bin*Rpi_bin+j*Rpi_bin+k] = xi_cc[i*Rp_bin*Rpi_bin+j*Rpi_bin+k] + 
								p_cen[i*Xi_nbin+p]*p_cen[i*Xi_nbin+q]*
								Xicc[i*Rp_bin*Rpi_bin*Nwp+j*Rpi_bin*Nwp+k*Nwp+lwp]
								*N_h[i*Xi_nbin+p]*N_h[i*Xi_nbin+q]/n_halo[i]/n_halo[i];

								xi_cs[i*Rp_bin*Rpi_bin+j*Rpi_bin+k] = xi_cs[i*Rp_bin*Rpi_bin+j*Rpi_bin+k] + 
								p_cen[i*Xi_nbin+p]*p_cen[i*Xi_nbin+q]*
								Xics[i*Rp_bin*Rpi_bin*Nwp+j*Rpi_bin*Nwp+k*Nwp+lwp]
								*N_h[i*Xi_nbin+p]*N_s[i*Xi_nbin+q]/n_halo[i]/n_halo[i];
						
								xi_sc[i*Rp_bin*Rpi_bin+j*Rpi_bin+k] = xi_sc[i*Rp_bin*Rpi_bin+j*Rpi_bin+k] + 
								p_cen[i*Xi_nbin+p]*p_cen[i*Xi_nbin+q]*
								Xisc[i*Rp_bin*Rpi_bin*Nwp+j*Rpi_bin*Nwp+k*Nwp+lwp]
								*N_s[i*Xi_nbin+p]*N_h[i*Xi_nbin+q]/n_halo[i]/n_halo[i];
						
								xi_ss[i*Rp_bin*Rpi_bin+j*Rpi_bin+k] = xi_ss[i*Rp_bin*Rpi_bin+j*Rpi_bin+k] + 
								p_cen[i*Xi_nbin+p]*p_cen[i*Xi_nbin+q]*
								Xiss[i*Rp_bin*Rpi_bin*Nwp+j*Rpi_bin*Nwp+k*Nwp+lwp]
								*N_s[i*Xi_nbin+p]*N_s[i*Xi_nbin+q]/n_halo[i]/n_halo[i];
//							}
								lwp = lwp + 1;

						}
					}

					xi_gg[i*Rp_bin*Rpi_bin+j*Rpi_bin+k] = xi_cc[i*Rp_bin*Rpi_bin+j*Rpi_bin+k] + 
					xi_cs[i*Rp_bin*Rpi_bin+j*Rpi_bin+k] + xi_sc[i*Rp_bin*Rpi_bin+j*Rpi_bin+k] + 
					xi_ss[i*Rp_bin*Rpi_bin+j*Rpi_bin+k];
//					printf("lwp=%d\n", lwp);
//					printf("xi_cc[%d,%d,%d]=%le\n", i,j,k, xi_cc[i*Rp_bin*Rpi_bin+j*Rpi_bin+k]);
//					printf("xi_gg[%d,%d,%d]=%le\n", i,j,k, xi_gg[i*Rp_bin*Rpi_bin+j*Rpi_bin+k]);
				}

			}	
	
		}
		else
		{
			for(j=0; j<Rp_bin; j++)
			{
				for(k=0; k<Rpi_bin; k++)
				{
					xi_gg[i*Rp_bin*Rpi_bin+j*Rpi_bin+k] = 0;
					xi_cc[i*Rp_bin*Rpi_bin+j*Rpi_bin+k] = 0; 
					xi_cs[i*Rp_bin*Rpi_bin+j*Rpi_bin+k] = 0;
					xi_sc[i*Rp_bin*Rpi_bin+j*Rpi_bin+k] = 0; 
					xi_ss[i*Rp_bin*Rpi_bin+j*Rpi_bin+k] = 0;					
				}
			}
		}
	}

/*	
	double *wpcc, *wpgg;
	wpcc = (double *)malloc(Mag_min_xi_all*Rp_bin*sizeof(double));
	wpgg = (double *)malloc(Mag_min_xi_all*Rp_bin*sizeof(double));

	for(i=0; i<Mag_min_xi_all; i++)
	{
		wpcc[i] = 0;
	}

	for(i=0; i<Mag_min_xi_all; i++)
	{
		for(j=0; j<Rp_bin; j++)
		{
			for(k=0; k<Rpi_bin; k++)
			{

				wpcc[i*Rp_bin+j] = wpcc[i*Rp_bin+j] + 4*xi_cc[i*Rp_bin*Rpi_bin+j*Rpi_bin+k];
				wpgg[i*Rp_bin+j] = wpgg[i*Rp_bin+j] + 4*xi_gg[i*Rp_bin*Rpi_bin+j*Rpi_bin+k];			
			}

			printf("wpcc[%d,%d]=%lf\n", i, j, wpcc[i*Rp_bin+j]);
//			printf("wpgg[%d,%d]=%lf\n", i, j, wpgg[i*Rp_bin+j]);

		}	
	}
*/

	table_x5 = gsl_integration_glfixed_table_alloc(node_x5);


   	FILE *fcorr_x;
    double *corr_x;
  
    corr_x = (double *)malloc(Rpi_bin*sizeof(double));

	fcorr_x = fopen("cubic_x", "r");

	for(i=0; i<Rpi_bin; i++)
	{
		fscanf(fcorr_x, "%lf", &corr_x[i]);
	}

    fclose(fcorr_x);


	for(s=0; s<Mag_min_xi_all; s++)
    {
/*
    FILE *fcorr_x;
    FILE *fcorr_2d_new;
    double corr_x[52]={0};
    double corr_y[52]={0};
    double corr_2d_new[52][28]={0};
	fcorr_x = fopen("cubic_x", "r");
	fcorr_2d_new = fopen("cubic_y_test", "r");

	for(i=0; i<52; i++)
	{
		fscanf(fcorr_x, "%lf", &corr_x[i]);
	}


	for(i=0; i<52; i++)
	{
		for(j=0; j<28; j++)
		{
			fscanf(fcorr_2d_new, "%lf", &corr_2d_new[i][j]);
		}
	}

	Rp_bin = 28;
	Rpi_bin=52;
	error_distance[s]=200;
*/
    	double *corr_y;
    	corr_y = (double *)malloc(Rpi_bin*sizeof(double));

		for(m=1; m<=wp_number; m++)
		{
			num = m-1;
			high_limit=50*m;

			for(i=0; i<node_x5; i++)
			{
				gsl_integration_glfixed_point(0,high_limit,    i,&xxi5[i],&xwi5[i],table_x5);			
//				printf("%lf\t%lf\n", xxi5[i], xwi5[i]);
			}
			for(j=0; j<Rp_bin; j++)	
			{ 
				w1[s*wp_number*Rp_bin+num*Rp_bin+j]=w2[s*wp_number*Rp_bin+num*Rp_bin+j]=0;
			}

			for(j=0; j<Rp_bin; j++)			
			{
//			for(k=0; k<52; k++)
//			{
//				corr_y[k] = corr_2d_new[k][j];
//			}

				for(k=0; k<Rpi_bin; k++)
				{
					corr_y[k] = xi_gg[s*Rp_bin*Rpi_bin+j*Rpi_bin+k];
//					printf("%lf\n", corr_y[k]);
				}

				gsl_interp_accel *acc9 = gsl_interp_accel_alloc ();
				gsl_spline *spline9 = gsl_spline_alloc (gsl_interp_cspline, Rpi_bin);
				gsl_spline_init (spline9, corr_x, corr_y, Rpi_bin);

				for(i=0; i<node_x5; i++)		
				{ 	
					w1[s*wp_number*Rp_bin+num*Rp_bin+j] = w1[s*wp_number*Rp_bin+num*Rp_bin+j]+gsl_spline_eval(spline9, xxi5[i], acc9)*gsl_sf_erf( (high_limit+xxi5[i])/sqrt(2)/error_distance[s] )*xwi5[i];
//					w2[num][j] = w2[num][j]+corr_2d[i][j]*2;
				}
				for(i=0; i<node_x5; i++)
				{
					w2[s*wp_number*Rp_bin+num*Rp_bin+j] = w2[s*wp_number*Rp_bin+num*Rp_bin+j]+gsl_spline_eval(spline9, xxi5[i], acc9)*gsl_sf_erf( (high_limit-xxi5[i])/sqrt(2)/error_distance[s] )*xwi5[i];
//					w4[num][j] = w4[num][j]+corr_2d[i][j]*2;
				}	
			
				w[s*wp_number*Rp_bin+num*Rp_bin+j] = w1[s*wp_number*Rp_bin+num*Rp_bin+j] + w2[s*wp_number*Rp_bin+num*Rp_bin+j]; 		
//				printf("w[%d][%d][%d]=%lf\t%lf\t%lf\n", s, num, j, w1[s*wp_number*Rp_bin+num*Rp_bin+j], w2[s*wp_number*Rp_bin+num*Rp_bin+j], w[s*wp_number*Rp_bin+num*Rp_bin+j]);
		

				gsl_spline_free (spline9);
				gsl_interp_accel_free (acc9);	
			}
	
		}

		free(corr_y);
	}
	for(k=0; k<Mag_min_xi_all; k++)
	{
		for(i=0; i<wp_number*Perpendicular_min; i++)
		{
			w_analysis[k*wp_number*Perpendicular_min+i] = w[k*wp_number*Rp_bin+(i/Perpendicular_min)*Rp_bin+(i-(i/Perpendicular_min)*Perpendicular_min+Start_point_xi-1)%Rp_bin];
//			printf("%d\n", k*wp_number*Rp_bin);
//			printf("%d\n", (i/Perpendicular_min)*Rp_bin);
//			printf("%d\n", (i+Start_point_xi-1)%Rp_bin );
//			printf("index=%d\n", k*wp_number*Rp_bin+(i/Perpendicular_min)*Rp_bin+(i-(i/Perpendicular_min)*Perpendicular_min+Start_point_xi-1)%Rp_bin );
//			printf("hope%le\n", w_analysis[k*wp_number*Perpendicular_min+i]);
		}
	}
	
	free(p_cen);
	free(xi_cc);
	free(xi_cs);
	free(xi_sc);
	free(xi_ss);
	free(xi_gg);
	free(n_halo);
    free(w);
    free(w1);
    free(w2);
    free(mag_abs_gauss);
    free(temp_mag);
    free(temp_mag_f);    
    free(temp_error_distance_high);
    free(temp_error_distance_low);
    free(error_distance_high);
    free(error_distance_low);
    free(error_distance);
    free(loc_error_distance_high);
    free(loc_error_distance_low);
}
