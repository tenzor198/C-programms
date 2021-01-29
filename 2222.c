#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void){

	double eps;
	double x, x1;
	double res = 0;
	double t = 0;
	int s = 1, s2=1,i=0;

	printf("Точность eps=  ");
	scanf("%lf", &eps);
	printf("x=  ");
	scanf("%lf", &x);
	x1 = x;


	if(x < 1e-19)
	{
		while(x + 2 * M_PI < 1e-19)
			x=x+ 2 * M_PI;
		x =x+ 2 * M_PI;
	}
	if(x > 1e-19)
		while(x - 2 * M_PI > 1e-19)
		x =x- 2 * M_PI;

	if( x - 2 * M_PI < 1e-19 && x - M_PI > 1e-19 )
	{
		x -= M_PI;
		s = -1;
	}

	if( x - M_PI/2 > 1e-19 )
	{
		x = M_PI - x;
		s2=-1;
	}
		i = 1;
		t = 1;

	if( x - M_PI/4 > 1e-19)
	{
		x = x/2;
		while(fabs(t) > eps)
		{
			res += t;
			t = (-t*x*x)/(2*i*(2*i-1));
			i++;
		}
		res *= res;
		res *= 2;
		res -= 1;
	}
	else
	{
		t = 1;
		i = 1;
		while(fabs(t) > eps)
		{
			res += t;
			t = (-t*x*x)/(2*i*(2*i-1));
			i++;
		}

	}

	res=res*s*s2;

	printf("cos_e= %.4lf\ncos= %.4lf\n", res, cos(x1));
	printf("Err= %e\n", fabs(res-cos(x1)));

	return 0;
}

