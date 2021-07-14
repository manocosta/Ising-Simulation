#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>

#define D_PROB 0.001
#define N_DESC 5000
#define N_PROM 5000
#define N N_DESC*N_PROM

#define N_REPS 500000
#define N_SUM_S_VECINOS 5
#define N_DSK 2
#define J_MAX 0.60001
#define D_J 0.001
#define B_MAX 2.00001
#define FILE_NAME "Ising(J).txt"


// gcc -Wall -O3 -o File.e File.c -lm
#include "Default_params.h"
#include "functions.h"

int main(int argc, char **argv)
{
	char file_path[200]=FILE_DIR, filename[100]=FILE_NAME, fullpath[300];
	// char mass_file[255];
	int i;
	int accept=0;
	double Magnetizacion; // Es la Magnetizacion total
	double Energia; // magn = Magnetizacion/(L*L) y energ = Energia/(L*L)
	double T, M_prom=0, accept_prom=0, E_prom=0, m_prom=0, e_prom=0;
	// double *d_Hs, w;
	// double *corr, corr_prom; //definimos el puntero como float para que pueda haber números decimales en la matriz.
	
	/* Seteo los parametros */
	bool init_HOT=true;
	float init_T=INIT_T, stop_T=STOP_T, d_T=D_T, prob=PROB, J=DEFAULT_J, B=DEFAULT_B;
	// double target_std=TARGET_STD;
	int seed=SEED, iRows=L, k_descor=K_DESCOR, k_termal=K_TERMALIZACION, N_meds=N_MEDS;
	get_arguments(argc, argv, &seed, &iRows, &init_T, &stop_T, &d_T, &J, &B,
			&prob, &init_HOT, &k_descor, &k_termal, &N_meds,
			file_path, filename, true);
	
	char bool_value[2][12] = {"false", "true"};
	int init_HOT_temp = 1;
    if (init_HOT==false) init_HOT_temp=0;
	int Temp_down=1;
    if (init_HOT==false) Temp_down=-1;
	int red[iRows*iRows];//, red0[iRows*iRows];
	T = init_T;

	// corr = (double*) malloc(iRows*iRows * sizeof(double));

	strcpy(fullpath, file_path);
	strcat(fullpath, filename);
	FILE *fp = fopen(fullpath,"w");
	fprintf(fp, "SEED = %d, L = %d, INIT_T = %f, STOP_T = %f, D_T = %f\n", seed, iRows, init_T, stop_T, d_T);
	fprintf(fp, "J = %f, B = %f, PROB = %f, INIT_HOT = %s\n FILE_NAME = %s\n", J, B, prob, bool_value[init_HOT_temp], filename);

	fprintf(fp, "T = %lf, J = %f, B = %f\n", T, J, B);  //,         M_prom/(L*L),     M^2_prom,        std_M         corr_prom,     p_accept_prom     L\n"
	
	fprintf(fp, "T\taccept\tEnergia\tenerg\tMagnetizacion\tmagn\n");
	
	// seed = time(NULL);

	poblar_ising(red, iRows, prob, seed, init_HOT);
	// printf("Matriz:\n(seed = %d)\n", (SEED));
	// print_matrix(red,iRows,iRows);
	
	Energia = calc_E_inicial(red, iRows, J, B);
	// energ = Energia/(iRows*iRows);
	
	Magnetizacion=0;
	for(i=0; i<iRows*iRows; i++){
	/* Ahora vamos a Calcular el valor inicial de la magnetizacion */
		Magnetizacion += *(red+i);
	}
	// magn = Magnetizacion/(iRows*iRows);
	
	double dH_tab[N_DSK*N_SUM_S_VECINOS], w_tab[N_DSK*N_SUM_S_VECINOS];;
	calc_dH_table(dH_tab, J, B);
	printf("dH:\n");
	print_d_matrix(dH_tab,N_DSK,N_SUM_S_VECINOS);
	
	calc_w_table(w_tab, J, B, T);
	printf("w:\n");
	print_d_matrix(w_tab,N_DSK,N_SUM_S_VECINOS);


	B = 0;
	/* Loop sobre Temperaturas */
	for (T=init_T; Temp_down*T>Temp_down*stop_T; T+=-Temp_down*d_T)
	{
		calc_dH_table(dH_tab, J, B);
		calc_w_table(w_tab, J, B, T);
		/* Termalizacion */
		termalizacion(red, iRows, w_tab, dH_tab, &Energia, &Magnetizacion, k_termal);
		// magn = Magnetizacion/(iRows*iRows);
		// energ = Energia/(iRows*iRows);
		// J = 0;
		/* Loop sobre mediciones */
		M_prom=0, accept=0, accept_prom=0, E_prom=0, m_prom=0, e_prom=0;
		for (i=0; i<N_meds; i++)
		{
			// accept = 0;
		/* Loop sobre pasos de Metropolis */
			accept += descorr_Metropolis(red, iRows, w_tab, dH_tab,
				&Energia, &Magnetizacion, k_descor);
			
			accept_prom = accept;
			M_prom += abs(Magnetizacion);
			E_prom += Energia;
		}
		M_prom = M_prom/(double)(N_meds);
		m_prom = M_prom/(double)(iRows*iRows);
		E_prom = E_prom/(double)N_meds;
		e_prom = E_prom/(double)(iRows*iRows);
		accept_prom = accept_prom/(double)(N_meds*k_descor*iRows*iRows);

		fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
				T, accept_prom, E_prom, e_prom, M_prom, m_prom);
		printf("T = %lf", T);
		
	} /* Final del Loop sobre Temperaturas */
	fclose(fp);
	// printf("Matriz:\n(seed = %d, p_accept = %.3lf, M = %.3lf)\n", (seed),((double)accept/N), M);
	// print_matrix(red,iRows,iRows);
	// printf("Press ENTER to finsh\n");
	// getchar();
	return 0;
}

#include "functions.c"