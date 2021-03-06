#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double dx(double t, double y) {

	return y;
}
double dy(double t, double x) {

        return -x;
}


int main() {
	double y=1.,x=0.,x_=0.,y_=1.,t=0.,h=0.01,err=0.,e=0.;
	double 	p[8] = {0.0,5179.0/57600, 0.0, 7571.0/16695, 393.0/640,-92097.0/339200, 187.0/2100, 1.0/40},
		p1[8]={0.,35.0/384.0, 0., 500.0/1113, 125.0/192, -2187.0/6784.0, 11.0/84.0, 0};
	double k[8],q[8],k1=0.,q1=0.,N=10000 * M_PI +1.e-13;;
	static const double a2=0.2,a3=0.3,a4=0.8,a5=8.0/9.0,
	b21=1.0/5,
	b31=3.0/40.0,b32=9.0/40.0,
	b41=44.0/45.0,b42=-56.0/15.0,b43=32.0/9.0,
	b51=19372.0/6561.0,b52=-25360.0/2187.0,b53=64448.0/6561.0,b54=-212.0/729.0,
	b61=9017.0/3168.0,b62=-355.0/33.0,b63=46732.0/5247.0,b64=49.0/176.0,b65=-5103.0/18656.0,
	b71=35.0/384.0, b72=0.0, b73=500.0/1113.0,b74=125.0/192.0,b75=-2187.0/6784.0,b76=11.0/84.0;
	FILE *fw;
	if(!(fw=fopen("res.txt","w")))
	{
		printf("Can't open\n");
		return -1;
	}
	for(;t<10000 * M_PI -h+1.e-13;)
	{
		k[0]=0;
		k[1] = h*dx(t, y);
		k[2] = h*dx(t + a2*h, y + b21*k[1]);
		k[3] = h*dx(t + a3*h, y + b31*k[1] + b32*k[2]);
		k[4] = h*dx(t + a4*h, y + b41*k[1] + b42*k[2] + b43*k[3]);
		k[5] = h*dx(t + a5*h, y + b51*k[1] + b52*k[2] + b53*k[3] + b54*k[4] );
		k[6] = h*dx(t + h, y + b61*k[1] + b62*k[2] + b63*k[3] + b64*k[4] + b65*k[5]);
		k[7] = h*dx(t + h, y + b71*k[1] + b72*k[2] + b73*k[3] + b74*k[4] + b75*k[5] + b76*k[6]);

		x+=p[1]*k[1] + p[2]*k[2] + p[3]*k[3] + p[4]*k[4] + p[5]*k[5] + p[6]*k[6] + p[7]*k[7];
//		x_+=p1[1]*k[1] + p1[2]*k[2] + p1[3]*k[3] + p1[4]*k[4] + p1[5]*k[5] + p1[6]*k[6] + p1[7]*k[7];


                q[0]=0;
                q[1] = h*dy(t,x);
		q[2] = h*dy(t + a2*h, x + b21*q[1]);
		q[3] = h*dy(t + a3*h, x + b31*q[1] + b32*q[2]);
		q[4] = h*dy(t + a4*h, x + b41*q[1] + b42*q[2] + b43*q[3]);
		q[5] = h*dy(t + a5*h, x + b51*q[1] + b52*q[2] + b53*q[3] + b54*q[4] );
		q[6] = h*dy(t + h, x + b61*q[1] + b62*q[2] + b63*q[3] + b64*q[4] + b65*q[5]);
		q[7] = h*dy(t + h, x + b71*q[1] + b72*q[2] + b73*q[3] + b74*q[4] + b75*q[5] + b76*q[6]);
/*
		for(i=1;i<8;i++){
			q1+=q[i]*p[i];
		}
		y_=y_+q1;
		y_+=p1[1]*q[1] + p1[2]*q[2] + p1[3]*q[3] + p1[4]*q[4] + p1[5]*q[5] + p1[6]*q[6] + p1[7]*q[7];
*/
		y+=p[1]*q[1] + p[2]*q[2] + p[3]*q[3] + p[4]*q[4] + p[5]*q[5] + p[6]*q[6] + p[7]*q[7];
		t+=h;



//		printf("%lf %lf\n", x, y);
		fprintf(fw, "%lf %lf\n",x, y);
	}

	printf("Погрешность= %lf\n", fabs(x*x+y*y-1));

	fclose(fw);
	return 0;
}
