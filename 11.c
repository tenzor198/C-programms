#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int Jor(double *mat,double *b,double *x,int n, int* c);
void SelectMainRow(double* A, int cur, int n, int* c);
void matxmat(double *mat1,double *mat2,double *mat,int n);
int main()
{

	int n=3,i,j,*ind,i1,i2,k,j2;
	double b[3]={1.0,
		     0.0,
                     0.0};

	double P[9]={7.24e6,0.0,0.0,
		     0.0,1.2e3,0.0,
		     0.0,0.0,9.2e-8};

	double D[9]={1.0/3,2.0/3,-2.0/3,
		    -2.0/3,2.0/3,1.0/3,
		     2.0/3,1.0/3,2.0/3};

	double D1[9]={1.0/3, -2.0/3, 2.0/3,
		      2.0/3,  2.0/3, 1.0/3,
		     -2.0/3, 1.0/3, 2.0/3};

        double *y,*x,*x_alp,cur,cur1,*A_alp, alp,*c,*A,*tmp,*tmp2,l,r;
	FILE *file;
	scanf("%lf %lf",&l, &r);
	file = fopen("answer.txt","w");

       	ind =(int*)malloc(n*sizeof(int));
       	c =(double*)malloc(n*sizeof(double));
       	y =(double*)malloc(n*sizeof(double));
	x_alp =(double*)malloc(n*sizeof(double));
	tmp =(double*)malloc((n*n)*sizeof(double));
        tmp2 =(double*)malloc((n)*sizeof(double));
	A =(double*)malloc((n*n)*sizeof(double));
       	x =(double*)malloc(n*sizeof(double));

        for (i = 0; i < n; i++)
                for (j = 0; j < n; j++) 
                	tmp[i*n+j] = A[i*n + j] = 0;

	matxmat(D, P, tmp, n);

	matxmat(tmp, D1, A, n);
///matD*tmp
        for(i=0; i<n; i++)
		c[i]=0;
        for(i=0; i<n; i++)
                for(j=0; j<n; j++)
                        tmp2[i]=c[i]+=D[i*n + j] * b[j];

	Jor(P,b,y,n,ind);
//*                        tmp2[i]=c[i]+=D[i*n + j] * b[j];

	Jor(P,b,y,n,ind);
*/
        for(i=0; i<n; i++)
                for(j=0; j<n; j++) 
                        x[i]+=D[i*n + j] * y[j];

        	for(i1=0;i1<n;i1++) 
			tmp2[i1] = c[i1];

	for(alp = 0.001+l; alp <= 0.001+r; alp += 0.001)
	{
       		A_alp =(double*)malloc((n*n)*sizeof(double));
//            		cur=fabs(x[j2]-x_alp[j2]);
//			tmp2[i1] = c[i1];
        	for(i1=0;i1<n;i1++) 
			tmp2[i1] = c[i1];

        	for(i1=0;i1<n;i1++) 
		{
                	for(i2=0;i2<n;i2++) 
               			A_alp[i1*n + i2] = A[i1*n + i2];
		}
//                	A_alp[k*n + k] +=alp;
        	for(k=0;k<n;k++)
                	A_alp[k*n + k] +=alp;

		Jor(A_alp,tmp2,x_alp,n,ind);
		cur=cur1=0;
//////          cur=0;
        	for(j2=0; j2<n; j2++) 
		{
            		cur=fabs(x[j2]-x_alp[j2]);
            		if(cur>cur1)
				cur1=cur;
        	}
        	fprintf(file,"%f %e\n",alp, cur1);

		printf("%f %e\n",alp, cur1);
		free(A_alp);
	}

	free(A);
	free(x_alp);
	free(ind);
	free(x);
	free(tmp);
	free(c);
	free(y);
	free(tmp2);
	fclose(file);
        return 0;
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
	for(i=0; i<n; i++)
		x[i]=vek1[i];
	return 0;
}

void matxmat(double *mat1,double *mat2,double *mat,int n)
{
	int i,j,l;
	for ( i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			for (l = 0; l < n; l++) 
				mat[i*n + j] += mat1[i*n + l] * mat2[l*n + j];
}


