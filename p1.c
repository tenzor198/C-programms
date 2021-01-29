#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>

#define EPS 1e-17

void GenOfMat(double *A, int n, int formula);
void SelectMainRow(double* A, int cur, int n, int* r);
int Jordan(int n, double *a, double *b, double *x,int *c ,int my_rank, int total_threads);
double NormaNev(double *A, double *b, double *x, int n);
void PrintMat(double *A, int n,int count);
void synchronize(int total_threads);
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

typedef struct
{	int n;
	int count;
	int formula;
	double *A;
	double *b;
	double *b1;
	double *x;
	int *c;
	int my_rank;
	int total_threads;
} ARGS;


pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

void *Solution(void *p_arg)
{
	struct timespec time_start1, time_end1;


	ARGS *arg = (ARGS*)p_arg;

	clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time_start1);


	Jordan(arg->n, arg->A, arg->b, arg->x,arg->c, arg->my_rank, arg->total_threads);

	clock_gettime(CLOCK_THREAD_CPUTIME_ID, &time_end1);

	pthread_mutex_lock(&mutex);
	pthread_mutex_unlock(&mutex);


	time_end1 =  diff(time_start1, time_end1);
	printf("Time thread  = %lf "  " %d\n\n", time_end1.tv_sec+(time_end1.tv_nsec)/1000000000.0  ,  arg->my_rank+1);

	return NULL;
}

int main(int argc, char * argv[]){
	int i=0,
	    j=0;
	int n,formula,count,l,k;
	double* x;
	double* b;
	double* b1;
	double* A;
	double* B;
	int * c;
	FILE *f;

	int total_threads;
	pthread_t *threads;
	ARGS *args;
	struct timespec time_start, time_end;

	switch(argc) {
		case 2: {
			f=fopen(argv[1], "r");
 			if (f == NULL){
    				printf("Error: 2 \n");
				return -2;
			}
			if(fscanf(f, "%d", &total_threads)==EOF) {
				printf("Error: 3a \n");
				return -3;
			}
 			if(total_threads <= 0) {
       				printf("Error: 3b \n");
				return -3;
			}

			if(fscanf(f, "%d", &n)==EOF) {
				printf("Error: 4 \n");
				return -4;
			}
			if(n <= 0) {
				printf("Error: 5 \n");
				return -5;
			}
			if(fscanf(f, "%d", &count)==EOF){
				printf("Error: 6 \n");
				return -6;
			}
			if(count <= 0 || count>n) {
	 			printf("Error: 7 \n");
				return -7;
			}
			A = (double*)malloc( (n*n) * sizeof(double));
			B = (double*)malloc( (n*n) * sizeof(double));
			b = (double*)malloc(n * sizeof(double));
			b1 = (double*)malloc(n * sizeof(double));
			for(i=0; i<n; i++) {
				for(j=0; j<n; j++) {
					l=fscanf(f,"%lf", &A[i*n+j]);
					if(l==EOF || l==0) {
						printf("Error: 7a \n");
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
			fscanf(f, "%d", &total_threads);
			fscanf(f, "%d", &n);
			fscanf(f, "%d", &count);
			for(i=0; i<n; i++) {
				for(j=0; j<n; j++)
					fscanf(f,"%lf", &B[i*n+j]);
				fscanf(f,"%lf", &b1[i]);
			}
			fclose(f);
			break;
		}
		case 5: {
			sscanf(argv[1], "%d", & total_threads);
			sscanf(argv[2], "%d", & n);
			sscanf(argv[3], "%d", & count);
			sscanf(argv[4], "%d", & formula);

			A = (double*)malloc( (n*n) * sizeof(double));
			B = (double*)malloc( (n*n) * sizeof(double));

		 	for(i=0;i<n*n;i++)
				A[i]=B[i]=0;

			GenOfMat(A,n,formula);
			GenOfMat(B,n,formula);

			if(total_threads <= 0) {
				printf("Error: 3b \n");
				return -3;
			}

			if(n <= 0){
		 		printf("Error 9 \n");
				return -9;
			}
			if(count <= 0 || count>n){
		 		printf("Error 10 \n");
				return -10;
			}
			if(formula != 1 && formula!=2 ){
				printf("Error 2 \n");
				return -2;
			}

			b = (double*)malloc(n * sizeof(double));
			b1 = (double*)malloc(n * sizeof(double));
			x = (double*)malloc(n * sizeof(double));

			for(i=0;i<n;i++) {
				if(i==1)
					b[i]=x[i]=b1[i]=1;
				else
					b[i]=x[i]=b1[i]=0;
			}
			break;
		}
		default: {
			printf("Неверное количество аргументов\n");
			return -1;

		}
	}


	c = (int*)malloc(n * sizeof(int));
	threads = (pthread_t*)malloc(total_threads * sizeof(pthread_t));
	args = (ARGS*)malloc(total_threads * sizeof(ARGS));

	if (!(A && b && b1 && x && c && threads && args))
	{
		printf("Not enough memory!\n");

		if (A) free(A);
		if (b) free(b);
		if (b1) free(b1);
		if (x) free(x);
		if (c) free(c);
		if (threads) free(threads);
		if (args) free(args);

		return -6;
	}

	printf("\n");
	PrintMat(A,n,count);

	printf("b= ");
	for (i = 0; i < count-1; ++i)
		printf("%.2lf ",  b[i]);
	printf("%.2lf \n\n" ,  b[n-1]);




	for (i = 0; i < total_threads; i++)
	{
		args[i].n = n;
		args[i].count = count;
		args[i].formula = formula;
		args[i].A = A;
		args[i].b = b;
		args[i].b1 = b1;
		args[i].x = x;
		args[i].c = c;
		args[i].my_rank = i;
		args[i].total_threads = total_threads;
	}



	clock_gettime(CLOCK_MONOTONIC, &time_start);

	for (i = 0; i < total_threads; i++)
		if (pthread_create(threads + i, 0, Solution, args + i))
		{
			printf("Cannot create thread %d!\n", i);

			if (A) free(A);
			if (b) free(b);
			if (b1) free(b1);
			if (x) free(x);
			if (c) free(c);
			if (threads) free(threads);
			if (args) free(args);

			return -7;
		}

	for (i = 0; i < total_threads; i++)
		if (pthread_join(threads[i], 0))
		{
			printf("Cannot wait thread %d!\n", i);

			if (A) free(A);
			if (b) free(b);
			if (b1) free(b1);
			if (x) free(x);
			if (c) free(c);
			if (threads) free(threads);
			if (args) free(args);

			return -8;
		}
	clock_gettime(CLOCK_MONOTONIC, &time_end);

	for(j=0;j<n;j++)
			printf("%.2lf ",x[j]);
	printf("\n\n");


	printf("Norma = %e \n\n", NormaNev(B,b1,x,n));


	time_end =  diff(time_start, time_end);
	printf("Time  = %lf \n\n", time_end.tv_sec+(time_end.tv_nsec)/1000000000.0);

	free(c);
	free(A);
	free(b);
	free(b1);
	free(x);
	free(B);
	free(threads);
	free(args);

	return 0;
}

void synchronize(int total_threads)
{
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	static pthread_cond_t condvar_in = PTHREAD_COND_INITIALIZER;
	static pthread_cond_t condvar_out = PTHREAD_COND_INITIALIZER;
	static int threads_in = 0;
	static int threads_out = 0;


	pthread_mutex_lock(&mutex);

	threads_in++;
	if (threads_in >= total_threads)
	{
		threads_out = 0;
		pthread_cond_broadcast(&condvar_in);
	} else
		while (threads_in < total_threads)
			pthread_cond_wait(&condvar_in,&mutex);

	threads_out++;
	if (threads_out >= total_threads)
	{
		threads_in = 0;
		pthread_cond_broadcast(&condvar_out);
	} else
		while (threads_out < total_threads)
			pthread_cond_wait(&condvar_out,&mutex);

	pthread_mutex_unlock(&mutex);

}

int Jordan(int n, double *a, double *b, double *x, int *index, int my_rank, int total_threads)
{
	int i, j, k;
	int first_row;
	int last_row;
	double tmp;

	for (i = 0; i < n; i++)
		index[i] = i;

	for (i = 0; i < n; i++)
	{
		if (my_rank == 0)
		{
			k = i;
			for (j = i+1 ; j < n; j++)
				if (fabs(a[i * n + k]) < fabs(a[i * n + j]))
					k = j;

			j = index[i];
			index[i] = index[k];
			index[k] = j;


			for (j = 0; j < n; j++)
			{
				tmp = a[j * n + i];
				a[j * n + i] = a[j * n + k];
				a[j * n + k] = tmp;
			}

			tmp = a[i * n + i];

			tmp = 1.0/tmp;
			for (j = i; j < n; j++)
				a[i * n + j] *= tmp;
			b[i] *= tmp;
		}
		synchronize(total_threads);

		first_row = i * my_rank;
		first_row = first_row/total_threads;
		last_row = i * (my_rank + 1);
		last_row = last_row/total_threads;

		for (j = first_row; j < last_row; j++)
		{
			tmp = a[j * n + i];
			for (k = i; k < n; k++)
				a[j * n + k] -= tmp * a[i * n + k];
			b[j] -= tmp * b[i];
		}

		first_row = (n - i-1 ) * my_rank;
		first_row = first_row/total_threads + i + 1;
		last_row = (n - i-1 ) * (my_rank + 1);
		last_row = last_row/total_threads + i + 1;

		for (j = first_row; j < last_row; j++)
		{
			tmp = a[j * n + i];
			for (k = i; k < n; k++)
				a[j * n + k] -= tmp * a[i * n + k];
			b[j] -= tmp * b[i];
		}
		synchronize(total_threads);
	}

	if (my_rank == 0)
		for (i = n-1; i >=0; i--)
			x[index[i]] = b[i];

	return 0;
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
	if(count!=n)
		printf("%.2lf \n\n",A[(n-1)*n+n-1]);
}

