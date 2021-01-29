#include <stdio.h>
#include <math.h>
double Obus(int n);
int main(void){
	FILE *file;
	int n, N,i;
	if ((file = fopen("result.txt","w")) == NULL)
	        printf("ERROR\n");
	scanf("%d %d", &n, &N);
	for(i=1; i<N; i++){
		printf("%d %.2lf\n", i+n, Obus(i+n));
		fprintf(file,"%d %lf\n", i+n, Obus(i+n));
	}
	return 0;
}
double Obus(int N){
	double max,min,h,cur=0;
	int i;
	h = 1.0/N;
	max=min = (2.0/h*h)*(1-cos((M_PI)/N));
	for(i=1; i<N; i++){
		cur = (2.0/h*h)*(1-cos((M_PI*i)/N));
		if(fabs(max)<fabs(cur))
			max= cur;
		if(fabs(min)>fabs(cur))
			min= cur;
	}
	cur=fabs(max/min);

	return cur;
}
