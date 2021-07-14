
double trial(double x0,double delta)
{
	double p,x;
	p = myrand();
	x = 2.0*delta*(p-0.5)+x0;
	return x;
}

int correlation(double *c,double *x,int n)
{
	int i,k;
	double xi,xk,s0,s1,s2,sk;
	for (k=0;k<n;k++)
	{
		s0 = 0.0;
		s1 = 0.0;
		s2 = 0.0;
		sk = 0.0;
		for (i=0;i<N-n;i++)
		{
			xi = *(x+i);
			xk = *(x+i+k);
			s1 += xi/(double)(N-n);
			s0 += xi*xk/(double)(N-n);
			s2 += xi*xi/(double)(N-n);
			sk += xk/(double)(N-n);
		}
	// *(c+k) = (s0-s1*s1)/(s2-s1*s1);
		*(c+k) = (s0-s1*sk)/(s2-s1*s1);
	}
	return 1;
}

double corr_redes(int* red0,int* red)
{
	int i;
	double xi,xk,s0,s1,s2,sk, corr;
	// for (k=0;k<n;k++)
	// {
	s0 = 0.0;
	s1 = 0.0;
	s2 = 0.0;
	sk = 0.0;
	for (i=0;i<L*L;i++)
	{
		xi = *(red0+i);
		xk = *(red+i);
		s1 += xi/(double)(L*L);
		s0 += xi*xk/(double)(L*L);
		s2 += xi*xi/(double)(L*L);
		sk += xk/(double)(L*L);
	}
	// *(c+k) = (s0-s1*s1)/(s2-s1*s1);
	corr = (s0-s1*sk)/(s2-s1*s1);
	return corr;
}

double normal_dist(double* xs, int n, double delta, double x0)
{
	double p, x, w, accept=0;
	int i;	
	accept=0;
	*(xs)=x0;
	for(i=0;i<n-1;i++)
	{
		p = myrand();
		x = trial(x0,delta);
		w = exp(-0.5*(x*x-x0*x0));
		if (p<w) 
		{
			x0 = x;
			accept++;
		}
		*(xs+i+1)=x0;
	}
	accept = accept/N;
	return accept;
}
double standar_dev(double* xs, int n)
{
	int i;
	double std=0, mean_x=0, mean_x2=0;
	for (i=0; i<n; i++)
	{
		mean_x += *(xs+i);
		mean_x2 += *(xs+i)*(*(xs+i));
	}
	mean_x = mean_x/n;
	mean_x2 = mean_x2/n;
	std = sqrt(mean_x2-mean_x*mean_x);
	return std;
}


double calc_std2(double* xs, int n)
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

double mean_val(double* xs, int n)
{
	int i;
	double mean_x=0;
	for (i=0; i<n; i++)
	{
		mean_x += *(xs+i);
	}
	mean_x = mean_x/n;
	return mean_x;
}