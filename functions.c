
#include "Imp_sampling_funct.c"

void get_arguments(int argc, char** argv, 
				int *seed, int *iRows, float *init_T, float *stop_T, float *d_T,
				float *J, float *B, float *prob, bool *init_HOT,
				int *k_descor, int *k_termal, int *N_meds,
				char* file_path, char* filename, bool print)
{
/*
	Llamado al programa:
			./programa.exe SEED 260572 L 32 N_REPS 50000 PROB 0.5
	 argc es el numero de elementos en argv
	 argv son los argumentos con los que se llama el programa
*/
	int i, init_HOT_temp=1;
	char argument[50];
	char bool_value[2][12] = {"false", "true"};
	
	for (i=1; i<argc-1; i++){
		sscanf(argv[i], "%s", argument);
		if (strcmp(argument, "SEED")==0) sscanf(argv[i+1], "%d", seed);

		if (strcmp(argument,"L")==0) sscanf(argv[i+1], "%d", iRows);

		// if (strcmp(argument,"N_REPS")==0) sscanf(argv[i+1], "%d", n_reps);

		if (strcmp(argument,"PROB")==0) sscanf(argv[i+1], "%f", prob);

		if (strcmp(argument,"INIT_T")==0) sscanf(argv[i+1], "%f", init_T);

		if (strcmp(argument,"STOP_T")==0) sscanf(argv[i+1], "%f", stop_T);

		if (strcmp(argument,"J")==0) sscanf(argv[i+1], "%f", J);

		if (strcmp(argument,"B")==0) sscanf(argv[i+1], "%f", B);
	
		if (strcmp(argument,"FILE")==0) sscanf(argv[i+1], "%s", filename);

		if (strcmp(argument,"PATH")==0) sscanf(argv[i+1], "%s", file_path);

		if (strcmp(argument,"INIT_HOT")==0) {
			sscanf(argv[i+1], "%d", &init_HOT_temp);
			if (init_HOT_temp==0){
				*init_HOT = false;  //If init_HOT==0 entonces Isisng se inicia con todos los espines iguales (baja temp)
			}else{ 
				init_HOT_temp = 1;
				*init_HOT = true;
			}
		}
		if (strcmp(argument,"K_DESCOR")==0) sscanf(argv[i+1], "%d", k_descor);

		if (strcmp(argument,"K_TERMAL")==0) sscanf(argv[i+1], "%d", k_termal);

		if (strcmp(argument,"N_MEDS")==0) sscanf(argv[i+1], "%d", N_meds);
	}
	if (print==true)
	{
		printf("SEED = %d, L = %d, INIT_T = %f, STOP_T = %f, D_T = %f\n", *seed, *iRows, *init_T, *stop_T, *d_T);
		printf("J = %f, B = %f, PROB = %f, INIT_HOT = %s\n FILE_NAME = %s\n", *J, *B, *prob, bool_value[init_HOT_temp], filename);
	}
}


double calc_E_inicial(int *red, int iRows, float J, float B)
{
	int s, s_down, s_right;
	double E_init = 0;

	for (int i=0; i<iRows; i++)
	{
		for (int j=0; j<iRows; j++)
		{
			s = *(red+iRows*i+j);
			s_down = *(red + iRows*((i+1+iRows)%iRows) + j);
			s_right = *(red + iRows*i + (j+1+iRows)%iRows);
			E_init += -J*(s_down + s_right)*s - B*s; 
			// printf("i:%d, j:%d, s:%d, s_r:%d, i_d:%d, j_d:%d, s_d:%d, E=%f\n", i, j, s, s_right, ((i+1+iRows)%iRows), j, s_down, E_init);
		}
	}
	return E_init;
}

double myrandom()
{
	double r;

	r = (double)rand()/(double)RAND_MAX;

	return r;
}
double myrand()
{
	double r;

	r = (double)rand()/(double)RAND_MAX;

	return r;
}
double maximo(double a, double b)
{
	double max;

	if (a>b) max=a;
	else max=b;

	return max;
}
double minimo(double a, double b)
{
	double min;

	if (a<b) min=a;
	else min=b;

	return min;
}
double modulo(double x)
{
	double abs;

	if (x<0) abs=-x;
	else abs = x;

	return abs;
}

double poblar_red(double* red, int iRows, double prob, double seed)
{
	
	int i, n;
	n = iRows;
	srand(seed); //setea que el seed para rand() sea nuestro SEED
	for (i=0;i<n*n;i++) if (prob<myrandom()) *(red+i)=0; else *(red+i)=1;

	return 0;
}

double poblar_ising(int* red, int iRows, double prob, double seed, bool init_HOT)
{
	int i, n;
	n = iRows;
	srand(seed); //setea que el seed para rand() sea nuestro SEED
	if (init_HOT==true){
		for (i=0;i<n*n;i++) if (prob<myrandom()) *(red+i)=-1; else *(red+i)=1;
	}else{
		for (i=0;i<n*n;i++) *(red+i)=1;
	}
	return 0;
}
double print_d_matrix(double* red, int num_col, int num_lines)
{
	int i, j;

	for(i = 0; i<num_lines; i++)
	{
  		for(j = 0; j<num_col; j++)
  		{
    		printf("%E   ",*(red+ (i)+j*num_lines)); // \t es tabular. Para que me queden espaciados los elementos de cada fila en pantalla.
  		}
  		printf("\n");// Paso al renglón siguiente de la pantalla para la próxima fila.
	}
	return 0;
}
double print_matrix(int* red, int num_col, int num_lines)
{
	int i, j;

	for(i = 0; i<num_lines; i++)
	{
  		for(j = 0; j<num_col; j++)
  		{
    		printf("%i\t",*(red+ (i)+j*num_lines)); // \t es tabular. Para que me queden espaciados los elementos de cada fila en pantalla.
  		}
  		printf("\n");// Paso al renglón siguiente de la pantalla para la próxima fila.
	}
	return 0;
}

int calc_dH_table(double* d_Hs, double J, double B)
{
	int i, j, dsk, sum_s;
	for (i=0;i<N_DSK;i++)
	{
		for (j=0;j<N_SUM_S_VECINOS;j++)
		{
			dsk = (4*i - 2);
			sum_s = 2*(j-2);
			*(d_Hs+(i*N_SUM_S_VECINOS + j)) = dsk*(-sum_s*J-B);
		}			
	}
	return 0;
}

int calc_w_table(double* w, double J, double B, double T)
{
	int i, j, dsk, sum_s;
	double d_Hs;
	for (i=0;i<N_DSK;i++)
	{
		for (j=0;j<N_SUM_S_VECINOS;j++)
		{
			dsk = (4*i-2);
			sum_s = 2*(j-2); 
			d_Hs = dsk*(-sum_s*J - B);
			*(w + (i*N_SUM_S_VECINOS + j)) = exp(-d_Hs/T);
		}			
	}
	return 0;
}

double get_w(int ii, int jj, double* w_tab, int* red, int iRows)
{
	int s, s_up, s_left, s_right, s_down, h;
	s = *(red+iRows*ii+jj);
	s_up = *(red+iRows*(((ii-1)+iRows)%iRows)+jj);
	s_left = *(red+iRows*ii+(jj-1+iRows)%iRows);
	s_right = *(red+iRows*ii+(jj+1+iRows)%iRows);
	s_down = *(red+iRows*(((ii+1)+iRows)%iRows)+jj);
	h = (((2-2*s)/4)*N_SUM_S_VECINOS+(s_up+s_left+s_right+s_down)/2+2);
	// printf("i:%d, j:%d, s:%d, s_r:%d, i_d:%d, j_d:%d, s_d:%d, E=%f\n", ii, jj, s, s_right, ((ii+1+iRows)%iRows), jj, s_down, E_init);
	return *(w_tab+h);
}

double get_dH(int ii, int jj, double* dH_tab, int* red, int iRows)
{
	int s, s_up, s_left, s_right, s_down, h;
	s = *(red+iRows*ii+jj);
	s_up = *(red+iRows*(((ii-1)+iRows)%iRows)+jj);
	s_left = *(red+iRows*ii+(jj-1+iRows)%iRows);
	s_right = *(red+iRows*ii+(jj+1+iRows)%iRows);
	s_down = *(red+iRows*(((ii+1)+iRows)%iRows)+jj);
	h = (((2-2*s)/4)*N_SUM_S_VECINOS+(s_up+s_left+s_right+s_down)/2+2);
	// printf("i:%d, j:%d, s:%d, s_r:%d, i_d:%d, j_d:%d, s_d:%d, dE=%f\n", ii, jj, s, s_right, ((ii+1+iRows)%iRows), jj, s_down, *(dH_tab + h));

	return *(dH_tab + h);
}

int Metropolis_step(int *red, int iRows, double *w_tab, double *dH_tab,
				 double *Energia, double *Magnetizacion)
{
	int ii = rand()%iRows; // Numero random entre 0 y iRows
	int jj = rand()%iRows; // Numero random entre 0 y iRows
	int s = *(red+iRows*ii+jj); //El spin en que nos paramos
	double w = get_w(ii, jj, w_tab, red, iRows);
	double p = myrand();
	if (p<w)
	{
		*Energia += get_dH(ii, jj, dH_tab, red, iRows);
		*(red+iRows*ii+jj) = -s;
		*Magnetizacion += -2*s;
		// printf("Matriz:\n(seed = %d, i = %i)\n", (seed), i);
		// print_matrix(red,iRows,iRows);
		return 1; // El spin fue modificado
	}
	return 0; // No se modifico el spin
}

int update_vec(double *vec, int len, double new_val)
{
	int j;
	for (j=0; j<len-1; j++)
	{
		*(vec+j+1) = *(vec + j);
	}
	*vec = new_val;
	return 0;
}

int descorr_Metropolis(int *red, int iRows, double *w_tab, double *dH_tab,
				double *Energia, double *Magnetizacion, int k_descor)
{
	int i, accept=0, len_comp = k_descor*iRows*iRows;
	for (i=0; i<len_comp; i++)
	{
		accept += Metropolis_step(red, iRows, w_tab, dH_tab,
			Energia, Magnetizacion);
	}
	return accept;
}

int termalizacion(int *red, int iRows, double *w_tab, double *dH_tab,
				double *Energia, double *Magnetizacion, int k_termal)
{
	int i, accept=0, len_comp=k_termal*iRows*iRows;
	// double vec_M[len_comp], vec_E[len_comp];//, vec_M2[len_comp], vec_E2[len_comp];
	// double std2_M, std2_E, std2;
	for (i=0; accept<len_comp; i++)
	{
		accept += Metropolis_step(red, iRows, w_tab, dH_tab, Energia, Magnetizacion);
		// update_vec(vec_M, len_comp, *Magnetizacion);
		// update_vec(vec_M2, len_comp, *Magnetizacion*(*Magnetizacion));
		// update_vec(vec_E, len_comp, *Energia);
		// update_vec(vec_E2, len_comp, *Energia*(*Energia));
	}
	// std2_M = calc_std2(vec_M, len_comp);
	// std2_E = calc_std2(vec_E, len_comp);
	// std2 = 0.5*(std2_M+std2_E);
	// while (std2>target_std*target_std)
	// {
	// 	for (j=0; j<100; j++) Metropolis_step(red, iRows, w_tab, dH_tab, Energia, Magnetizacion);
	// 	accept += Metropolis_step(red, iRows, w_tab, dH_tab, Energia, Magnetizacion);
	// 	update_vec(vec_M, len_comp, *Magnetizacion);
	// 	update_vec(vec_E, len_comp, *Energia);
	// 	std2_M = calc_std2(vec_M, len_comp);
	// 	std2_E = calc_std2(vec_E, len_comp);
	// 	std2 = 0.5*(std2_M+std2_E);
	// 	i++;
	// }
	return accept;
}

double corr_Ising(int* red0,int* red, int iRed, 
			int sum_spins, int sum_spins0, double std_x0)
{
	int i;
	double xi_medio = sum_spins/(double)(iRed*iRed);
	double x0_medio = sum_spins0/(double)(iRed*iRed);
	double xi,x0,sum_up, std_xi=0.0, corr=0.0;
	// for (k=0;k<n;k++)
	// {
	sum_up = 0.0;
	xi = 0;
	x0 = 0;
	for (i=0;i<iRed*iRed;i++)
	{
		xi = *(red0+i);
		x0 = *(red+i);
		sum_up += (xi-xi_medio)*(x0-x0_medio);
		std_xi += (xi-xi_medio)*(xi-xi_medio);
	}
	std_xi = sqrt(std_xi);
	corr = (sum_up)/(std_xi*std_x0);
	return corr;
}

double calc_std2_Ising(int* xs, int n)
{
	int i;
	double std2=0, mean_x=0;
	for (i=0; i<n; i++)
	{
		mean_x += *(xs+i);
	}
	mean_x = mean_x/n;
	for (i=0; i<n; i++)
	{
		std2 += (*(xs+i)-mean_x)*(*(xs+i)-mean_x);
	}
	// mean_x2 = mean_x2/n;
	// std2 = (mean_x2-mean_x*mean_x);
	return std2;
}

