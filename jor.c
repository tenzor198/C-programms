#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#define EPS 1e-100

void GenOfMat(double *A, int n, int formula);
void SelectMainRow(double* A, int cur, int n, int* r);
int Jordan(double* A, double* x, int n, int* c, double* vek1);
double NormaNev(double *A, double *b, double *x, int n);
void PrintMat(double *A, int n,int count);
struct timespec diff(struct timespec start, struct timespec end);
struct timespec diff(struct timespec start, struct timespec end)
{
    struct timespec temp;
    if ((end.tv_nsec-start.tv_nsec)<0) {
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return temp;
}
int main(int argc, char * argv[]){
	int i=0,
	    j=0;
	int n,formula,count,l,k;
	double* x;
	double* b;
	double* A;
	double* B;
	int * c;
	double* vek1;
	FILE *f;
	struct timespec time_start, time_end;
	if (argc == 1 || argc==3 ||argc >4) {
		printf("Неверное количество аргументов\n");
		return -1;
	}
	switch(argc) {

		case 2: {
			f=fopen(argv[1], "r");
 			if (f == NULL){
    				printf("Error: 2 \n");
				return -2;
			}
			if(fscanf(f, "%d", &n)==EOF) {
				printf("Error: 3 \n");
				return -3;
			}
			if(n <= 0) {
				printf("Error: 4 \n");
				return -4;
			}
			if(fscanf(f, "%d", &count)==EOF){
				printf("Error: 5 \n");
				return -5;
			}
			if(count <= 0 || count>n) {
	 			printf("Error: 6 \n");
				return -6;
			}
			A = (double*)malloc( (n*n) * sizeof(double));
			B = (double*)malloc( (n*n) * sizeof(double));
			b = (double*)malloc(n * sizeof(double));
			for(i=0; i<n; i++) {
				for(j=0; j<n; j++) {
					l=fscanf(f,"%lf", &A[i*n+j]);
					if(l==EOF || l==0) {
						printf("Error: 7 \n");
						return -7;
					}
				}
				k=fscanf(f,"%lf", &b[i]);
				if(k==EOF || k==0) {
						printf("Error: 8 \n");
						return -8;
				}
			}
			x = (double*)malloc(n * sizeof(double));
			for(i=0;i<n;i++){
					x[i]=b[i];
			}
			fclose(f);
			f=fopen(argv[1], "r");
			fscanf(f, "%d", &n);
			fscanf(f, "%d", &count);
			for(i=0; i<n; i++) {
				for(j=0; j<n; j++)
					fscanf(f,"%lf", &B[i*n+j]);
				fscanf(f,"%lf", &b[i]);
			}
			fclose(f);
			break;
		}

		case 4: {
			sscanf(argv[1], "%d", & n);
			sscanf(argv[2], "%d", & formula);
			sscanf(argv[3], "%d", & count);
			A = (double*)malloc( (n*n) * sizeof(double));
			B = (double*)malloc( (n*n) * sizeof(double));
			GenOfMat(A,n,formula);
			GenOfMat(B,n,formula);
			if(n <= 0){
		 		printf("Error 1 \n");
				return -1;
			}
			if(formula != 1 && formula!=2 ){
				printf("Error 2 \n");
				return -2;
			}
			if(count <= 0 || count>n){
		 		printf("Error 3 \n");
				return -3;
			}
			b = (double*)malloc(n * sizeof(double));
			x = (double*)malloc(n * sizeof(double));
			for(i=0;i<n;i++){
				if(i==1)
					b[i]=x[i]=1;
				else
					b[i]=x[i]=0;
			}
			break;
		}

		default: {
			printf("Неверное количество аргументов\n");
			return -1;
		}

	}
	printf("\n");
	PrintMat(A,n,count);
	printf("\n");
	printf("b= ");
	c = (int*)malloc(n * sizeof(int));
	vek1 = (double*)malloc(n * sizeof(double));
	for(i=0;i<count;i++)
		printf("%.2lf ",b[i]);
	printf("\n");

	clock_gettime(CLOCK_MONOTONIC, &time_start);

	if(Jordan(A, x, n, c, vek1)== -1)
		printf("Матрица вырождена\n");
	else{
		printf("x= ");
		for(j=0;j<count;j++)
			printf("%.2lf ",x[j]);
		printf("\n");
	}
	clock_gettime(CLOCK_MONOTONIC, &time_end);
	printf("Norma = %e ", NormaNev(B,b,x,n));
	printf("\n");

	free(vek1);
	free(c);
	free(A);
	free(b);
	free(x);
	free(B);

	time_end =  diff(time_start, time_end);
	printf("Time  = %lf \n\n", time_end.tv_sec+(time_end.tv_nsec)/1000000000.0);
	return 0;
}

int Jordan(double* A, double* x, int n, int* c, double* vek1){
	int i,j,k;
	double tmp;
	for(i=0; i<n; i++)
		c[i] = i;
	for(i=0; i<n; i++){
		SelectMainRow(A, i, n, c);
		if( fabs( A[i*n+c[i]] )<EPS)
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

void GenOfMat(double *A, int n, int formula){
	int i,j;
	if (formula==1){
		for(i=0;i<n;i++){
			for(j=0;j<n;j++)
				A[i*n+j]=fabs(i-j);
		}
	}
	if(formula==2){
		for(i=0;i<n;i++){
			for(j=0;j<n;j++)
				A[i*n+j]=(double)(1./(1+i+j));
		}
	}
}

double NormaNev(double* A, double* b, double* x, int n){
	int i,j;
	double temp=0, ans=0;
	for(i=0; i<n; i++) {
		for(j=0; j<n; j++) {
			temp+=A[i*n+j]*x[j];
		}
		b[i]-=temp;
		temp = 0;
	}
	for(i=0; i<n; i++)
		ans+=fabs(b[i]);
	return ans;
}

void PrintMat(double *A, int n,int count){
	int i,j;
	for(i=0;i<count;i++){
		for(j=0;j<count;j++)
			printf("%.2lf ",A[i*n+j]);
		printf("\n");
	}
}

