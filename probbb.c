#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <string.h>
#include <mpi.h>

#define EPS 1e-16
void generation(double *A, int n, int formula)
{
	if (formula == 1)
	{
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				A[n*i + j] = fabs(i - j);
	}
	if (formula == 2)
	{
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				A[n*i + j] = 1 / (double)(i + j + 1);
	}
	if (formula == 3)
	{
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				A[i*n + j] = 1;
		//memset(A, 1, n*n*sizeof(double));
		for (int i = 0; i < n; i++)
			A[(i + 1)*(n - 1)] = -1;
	}
}

int main(int argc, char **argv)
{
        int rank, size;
  	MPI_Status status;
        int n, N, *c;
    	double *A, *x, *y, *X;
	n = atoi(argv[1]);
    	MPI_Init(&argc, &argv);
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    	MPI_Comm_size(MPI_COMM_WORLD, &size);
  	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
   	if (rank < n % size)
        	N = n/size + 1;
    	else
        	N = n/size;
//	printf("%d\n", nCols);

	MPI_Barrier(MPI_COMM_WORLD);

    	A = (double*)malloc(n*N*sizeof(double));
    	y = (double*)malloc(n*sizeof(double));
    	c = (int*)malloc(n*sizeof(int));
    	x = (double*)malloc(n*sizeof(double));	//res
    	X = (double*)malloc(n*sizeof(double));	//res
        MPI_Barrier(MPI_COMM_WORLD);

	if(rank == 0)
	{
		int k = 0;
		double *AA;
		AA = (double*)malloc(n*n*sizeof(double));
		generation(AA, n, 1);
		for(int i=0; i<n; i++)
		{
			for(int j=0; j<n; j++)
				printf("%.2lf ", AA[i*n + j]);
			printf("\n");
		}
			printf("\n");
		int  temp = 0;
		int n1 = 0;
		temp+= N;
		for(int i = 0; i < size; i++)
		{
			
			if(i%size != 0)
			{
	 			MPI_Recv(&n1 , 1, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//				printf("nCols = %d\n", n1);
				for(int j=0; j<n; j++)
				{
//					printf("%.1lf\n", AA[temp + j*n]);
					MPI_Send(AA + temp + j*n, n1,  MPI_DOUBLE, i%size, 667, MPI_COMM_WORLD);
				}
				temp += n1;
			}
			else
			{
				for(int i = 0; i < n; i++)
					for(int j=0; j < N; j++)
						A[i*N + j] = AA[i*n + j];
			}
		}
	}
	else
	{
		MPI_Send(&N, 1,  MPI_DOUBLE, 0, 667, MPI_COMM_WORLD);
		for(int i=0; i<n; i++)
		{
 			MPI_Recv(A + i*N  , N, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//			printf("poluchil: %.1lf\n", A[i*nCols]);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0)
	{
		for(int i=0; i<n; i++)
		{
			for(int j=0; j<N; j++)
				printf("%.2lf ", A[i*N+j]);
			printf("\n");
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
			printf("\n");
			printf("\n");

	if(rank == 1)
	{
		for(int i=0; i<n; i++)
		{
			for(int j=0; j<N; j++)
				printf("%.2lf ", A[i*N+j]);
			printf("\n");
		}
	}

	for(int i=0; i<n; i++)
		c[i] = i;
	int i,j,k,maxc,numb,p;
	double tmp, maxv, orB;
	for(p=0; p < size; p++) {

		if(p == rank)
		{

			for(k=0; k < N;  k++)
			{

				maxv = A[(N*p + k-1)*n +k];
				maxc = N*p +k-1;

					printf("%.2lf =maxv   ", maxv);
					printf("%d =maxc   ", maxc);
				for( i = N*p + k + 1; i < n; i++)

					if( fabs(A[i*n + k]) > fabs(maxv) )
					{
						maxv = A[i*n + k];
						maxc = i;
//					printf("%.2lf =maxv   ", maxv);
//					printf("%d =maxc   ", maxc);

					}

				for( i = 0; i < N; i++)
				{
					tmp = A[ (N*p + k)*n + i];
					A[ (N*p+k)*n + i] = A[ maxc*n + i];
					A[ maxc*n +i] = tmp;
				}

				if( maxc != N*p + k)
				{
					tmp = c[N*p + k];
					c[N*p + k] = c[maxc];
					c[maxc] = tmp;
					tmp = x[N * p + k];
					x[N * p + k] = x[maxc];
					x[maxc] = tmp;
				}

				for(j = N*p + k; j < n-1; j++) {
					y[j] = A[(j+1)*n + k]/A[(p*N+k)*n+k];
			//		printf("%.2lf dxfcgvhbjnkm   ", A[(p*N+k)*n+k]);

				}
		               	for(int proc = p + 1;proc < size; proc++)
                		{
                   			MPI_Send(y,n,MPI_DOUBLE,proc,1,MPI_COMM_WORLD);
                    			MPI_Send(&k,1,MPI_DOUBLE,proc,2,MPI_COMM_WORLD);
                    			MPI_Send(&maxc,1,MPI_INT,proc,22,MPI_COMM_WORLD);
                    			MPI_Send(c,n,MPI_INT,proc,33,MPI_COMM_WORLD);
                		}
////////////////////////


				for(int i = N*p + k + 1; i < n; i++)
	        		{
            				for(int j = k; j < N; j++)
					{
                				A[i*n+j] = A[i*n + j] - A[(N*p + k)*n + j] * y[i - 1];
            					if(fabs(A[i*n + j]) < 1e-10)
							A[i*n + j] = 0;
                        		}
            				x[i] = x[i] - x[N*p + k] * y[i - 1];
        			}

        			if(p != size - 1)
				{
            				for(int proc = p + 1;proc < size; proc++)

                			MPI_Send(x,n,MPI_DOUBLE,proc,3,MPI_COMM_WORLD);
        			}

    			}
    		}
/*	if(rank == 0)
	{
		for(int i=0; i<n; i++)
		{
			for(int j=0; j<N; j++)
				printf("%.2lf ", A[i*N+j]);
			printf("\n");
		}
	}
*/
    	if(rank > p)
   	{
		for(int k = 0; k < N; k++)
		{
			MPI_Recv(&numb,1,MPI_DOUBLE,p,2,MPI_COMM_WORLD,&status);
        	        MPI_Recv(y,n,MPI_DOUBLE,p,1,MPI_COMM_WORLD,&status);
               		MPI_Recv(x,n,MPI_DOUBLE,p,3,MPI_COMM_WORLD,&status);
                	MPI_Recv(c,n,MPI_INT,p,33,MPI_COMM_WORLD,&status);
                	MPI_Recv(&maxc,1,MPI_INT,p,22,MPI_COMM_WORLD,&status);
                	for(int i = 0; i < N; i++)
			{
                    		tmp = A[maxc*n + i];
                    		A[maxc*n + i] = A[(N*p + numb)*n + i];
                    		A[(N*p + numb)*n + i] = tmp;
			}
                	for(int i = numb + 1; i < n; i++)
                        {
                		for(int j = numb; j < N; j++)
                                	A[i*n + j] = A[i*n + j] - A[k*n + j] * y[i - 1];//к потому что по диагонали
                            //B[i] = B[i] - B[k]*L[i - 1];
                        }
            	}

    	}
	}


//-------так как последний свободный вектор, верный
	MPI_Bcast(x,n,MPI_DOUBLE,size - 1,MPI_COMM_WORLD);
	MPI_Bcast(c,n,MPI_INT,size - 1, MPI_COMM_WORLD);
	orB = x[0];
//-----обратный ход
	for(int p = size - 1; p >=0; p--)
	{
    		if(p == rank)
        	{

            		X[N*p + N - 1] = x[N*p + N - 1]/ A[(N*p + N-1)*n + (N-1)];

			for(int i = N - 2; i >= 0; i--){
				for(int j = i + 1; j < N; j++)
					x[N*p + i] = x[N*p + i] - A[(N*p + i)*n + j] * X[N*p+j];
				X[N*p + i] = x[N*p + i]/A[(N*p + i)*n + i];
                }
                for(int i = 0; i < N*p; i++)
                    for(int j = 0; j < N;j++)
                        x[i] = x[i] - A[i*n + j] * X[N*p + j];
                if(p>0)
		{
                	MPI_Send(x,n,MPI_DOUBLE,p-1,4,MPI_COMM_WORLD);
                	MPI_Send(X,n,MPI_DOUBLE,p-1,5,MPI_COMM_WORLD);
		}

        }
    	if(rank == p - 1)
    	{
            MPI_Recv(x,n,MPI_DOUBLE,p,4,MPI_COMM_WORLD,&status);
            MPI_Recv(X,n,MPI_DOUBLE,p,5,MPI_COMM_WORLD,&status);

    	}
	}


	double ch = 0;
	MPI_Bcast(X,n,MPI_DOUBLE,0,MPI_COMM_WORLD);

//printVector(count,X);
	for(int p = 0; p < size; p++)
	{
    		if(rank == p) 
		{
        		for(int i = 0; i < N; i++)
            			ch = ch +  A[0*n + i] * X[N*p + i];
        		if(p != size - 1)
				MPI_Send(&ch,1,MPI_DOUBLE,p + 1,223,MPI_COMM_WORLD);
        		if(p == size - 1)
				printf("pogreshnost = %2.16f\n",orB - ch);
    		}
    		if(rank == p + 1)
        		MPI_Recv(&ch,1,MPI_DOUBLE,p,223,MPI_COMM_WORLD,&status);
    	}













    	MPI_Finalize();
	return 0;
}
