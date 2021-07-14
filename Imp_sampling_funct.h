#ifndef IMP_SAM_FUNCT
#define IMP_SAM_FUNCT
double trial(double x0,double delta);
int correlation(double *c,double *x,int n);
double corr_redes(int* red0,int* red);
double normal_dist(double* xs, int n, double delta, double x0);
double standar_dev(double* xs, int n);
double calc_std2(double* xs, int n);
double mean_val(double* xs, int n);
#endif

#ifndef DELTA
#define DELTA 0.75
#endif

#ifndef N
#define N 10000
#endif


#include <math.h>
#ifndef MPI
#define M_PI 3.14159265358979323846
#endif
