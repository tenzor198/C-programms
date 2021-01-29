#include <stdio.h>
#include <math.h>
#include <stdlib.h>
double Norm(double *x, double *y, int n);
void SelectMainRow(double* A, int cur, int n, int* c);
int Jor(double* P,double* x,double* vek1,int n,int* c);
void matxmat(double *A, double *B, int n);

int main ()
{
	double P[9]={7.24e6,0.0,0.0,
		     0.0,1.2e3,0.0,
		     0.0,0.0,9.2e-8 };
	double b[3]={1.0,
		     0.0,
                     0.0 };
	double D[9]={1.0/3,2.0/3,-2.0/3,
		    -2.0/3,2.0/3,1.0/3,
		     2.0/3,1.0/3,2.0/3 };
	double D1[9]={1.0/3, -2.0/3, 2.0/3,
		      2.0/3,  2.0/3, 1.0/3,
		     -2.0/3, 1.0/3, 2.0/3};
	double alp =1.0/7.24e6;
	double alpopt = 2.0/(7.24e6 + 9.2e-8);
	int n=3,*ind,N,i,j,k,q;
	double *ans, *x, *xopt, *x_cot,*tp, *xopt_cot;

	FILE *file1;
	FILE *file2;

    	ans = (double*)malloc(n*sizeof(double));
       	ind =(int*)malloc(n*sizeof(int));
	tp = (double*)malloc(n*sizeof(double));
	x = (double*)malloc(n*sizeof(double));
	xopt = (double*)malloc(n*sizeof(double));
	x_cot = (double*)malloc(n*sizeof(double));
	xopt_cot = (double*)malloc(n*sizeof(double));
	for(i=0; i<n; i++)
	{
		ind[i]=0;
		x[i]=0;
		xopt[i]=0;
		tp[i]=0;
	}
	printf("N: ");
	scanf("%d",&N);
	file1 = fopen("answ1.txt","w");
	file2 = fopen("answ2.txt","w");


//	Jor(P,b,ans,n,ind);
/*
	for(i=0;i<3;i++)
		printf("%lf\n",ans[i]); 

*/
//	ans[0]=	1.0/(7.24e6);
/*	ans[0]=	0.0;

	ans[1]=ans[2]=0.0;

	b[0]=0.0;
	b[1]=b[2]=0.0;
*/
	for(i=0; i<n; i++)
	{
		for(j=0; j<n; j++)
			tp[i] += D[i*n+j] *b[j];
	}
	for(i = 0; i<n; i++) 
	{
		b[i] = tp[i];
		tp[i]=0;
	}

	matxmat(D, P, n);
//	matxmat(D, D1, n);
	matxmat(D, D1, n);

	Jor(D,b,ans,n,ind);


	for(q=1; q<N+1; q++)
	{
		for(i=0; i<n; i++)
		{
			x_cot[i] = x[i];
			xopt_cot[i] = xopt[i];
		}


		for(i=0; i<n; i++) 
		{
			for(j=0; j<n; j++)
				tp[i] += D[i*n+j] *x_cot[j];
		}
		for(i = 0; i<n; i++)
		{
			x_cot[i] = tp[i];
			tp[i]=0;
		}

		for(i=0; i<n; i++) 
		{
			for(j=0; j<n; j++)
				tp[i] += D[i*n+j] *xopt_cot[j];
		}
		for(i = 0; i<n; i++) 
		{
			xopt_cot[i] = tp[i];
			tp[i]=0;
		}



		for(i=0; i<n; i++) 
		{
			x_cot[i]-=b[i];
			xopt_cot[i]-=b[i];
		}
		for(i=0; i<n; i++)
		{
			x_cot[i] *= alp;
			xopt_cot[i] *= alpopt;
		}
		for(i=0; i<n; i++)
		{
			x[i]-=x_cot[i];
			xopt[i]-=xopt_cot[i];
		}
		fprintf(file1, "%d %e %e\n", q, Norm(ans,xopt,n), Norm(ans,x,n));
//		fprintf(file1, "%d %e %e\n", q, Norm(ans,x,n),Norm(ans,xopt,n));
	}

	for(i=0;i<n;i++) {
		x[i]=0;
		x_cot[i]=0;
	}


	for(alp=0; alp<0.000002; alp+=0.0000001){

		for(j=1; j<N+1; j++){
			for(i=0; i<n; i++)
				x_cot[i] = x[i];

			matxmat(D, x_cot, n);


			for(i=0; i<n; i++)
				x_cot[i] -= b[i];

			for(i=0; i<n; i++)

				x_cot[i] *= alp;
			for(i=0; i<n; i++)
				x[i] -=x_cot[i];
		}

		fprintf(file2, "%e %e\n", alp, Norm(ans,x,n));
		for(q=0;q<n;q++)
			x[q]=0;
	}





	fclose(file1);
	free(ans);
	free(tp);
	free(ind);
	free(x);
	free(x_cot);
	free(xopt);
	free(xopt_cot);
	return 1;
}
void SelectMainRow(double* A, int cur, int n, int* c){
	int i, maxrow, temp;
	maxrow = cur;
	for(i=cur; i<n; i++)
		if(fabs(A[cur*n+c[i]]) > fabs(A[cur*n + c[maxrow]])){
			maxrow = i;
		}
	temp=c[cur];
	c[cur]=c[maxrow];
	c[maxrow] = temp;
}


int Jor(double* A,double* x,double* vek1,int n,int* c){
	int i,j,k;
	double tmp;
	for(i=0; i<n; i++)
		c[i] = i;
	for(i=0; i<n; i++){
		SelectMainRow(A, i, n, c);
		if( fabs( A[i*n+c[i]] )<1e-18)
			return -1;
		tmp = A[i*n+c[i]];
		x[i] /= tmp;
		for(j=0; j<n; j++)
			A[i*n+c[j]]/=tmp;
		for(j=0; j<n; j++){
			if(j!=i){
				tmp = A[j*n+c[i]];
				x[j] -= x[i]*tmp;
				for(k=i; k<n; k++)
					A[j*n+c[k]] -= A[i*n+c[k]]*tmp;
			}
		}
	}
	for(i=0; i<n; i++)
		vek1[c[i]]=x[i];
	for(i=0; i<n; i++) {
		x[i]=vek1[i];
	}
	return 0;
}

void matxmat(double *A, double *B, int n){
	int i,j,k;
	double *C = (double*)malloc((n*n)*sizeof(double));


	for(i = 0; i<n; i++)
		for(j=0; j<n; j++)
		{
			C[i*n + j] = 0;
			for(k=0; k<n; k++)
				C[ i*n+j ]+=A[i*n+k]*B[k*n+j];
		}
	for(i = 0; i<n; i++)
		for(j=0; j<n; j++)
			A[i*n + j]=C[ i*n + j];
}

double Norm(double * x, double *y, int n){
	int i=0;
	double tmp=0,tmp1=0;
       	for(i=0;i<n;i++){
    		tmp=fabs(x[i]-y[i]);
            if(tmp>tmp1) tmp1=tmp;
       	}
	return tmp1;
}

