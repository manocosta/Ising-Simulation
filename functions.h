#ifndef N_SUM_S_VECINOS
#define N_SUM_S_VECINOS 5
#endif

#ifndef N_DSK
#define N_DSK 2
#endif

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

void get_arguments(int argc, char** argv, 
				int *seed, int *iRows, float *init_T, float *stop_T, float *d_T,
				float *J, float *B, float *prob, bool *init_HOT,
				int *k_descor, int *k_termal, int *N_meds,
				char* file_path, char* filename, bool print);
				
double calc_E_inicial(int *red, int iRows, float J, float B);
double myrandom();
double myrand();
double maximo(double a, double b);
double minimo(double a, double b); 
double modulo(double x);
double poblar_red(double* red, int size_red, double prob, double seed);
double poblar_ising(int* red, int size_red, double prob, double seed, bool init_HOT);
double print_matrix(int* red, int num_col, int num_lines);
double print_d_matrix(double* red, int num_col, int num_lines);
int calc_dH_table(double* d_Hs, double J, double B);
int calc_w_table(double* w, double J, double B, double T);
double get_w(int ii, int jj, double* w_tab, int* red, int iRows);
double get_dH(int ii, int jj, double* dH_tab, int* red, int iRows);
int Metropolis_step(int *red, int iRows, double *w_tab, double *dH_tab,
				 double *Energia, double *Magnetizacion);
int update_vec(double *vec, int len, double new_val);
int termalizacion(int *red, int iRows, double *w_tab, double *dH_tab,
				double *Energia, double *Magnetizacion, int k_termal);
double corr_Ising(int* red0,int* red, int iRed, 
			int sum_spins, int sum_spins0, double std_x0);
double calc_std2_Ising(int* xs, int n);
int descorr_Metropolis(int *red, int iRows, double *w_tab, double *dH_tab,
				double *Energia, double *Magnetizacion, int N_descor);


#include "Imp_sampling_funct.h"




// double etiquetar_frag(double* red, int size_red, double* clase);
// double descartar_frag_perc(double* red, int size_red, double* top, double* bottom, double* perc, double* n_perc);
// void size_of_frag(double *y,double *x,int n,double a,double b,int m);
// void histograma(double *y,double *x,int n,double a,double b,int m);
#endif

