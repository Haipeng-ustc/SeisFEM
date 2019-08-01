#include <stdio.h>
#include <stdlib.h>

#define N 10000000
int main()
{
	int *a = (int *)malloc(N * sizeof(int));
	int *b = (int *)malloc(N * sizeof(int));
	int *c = (int *)malloc(N * sizeof(int));
	int i;

	for(i=0;i<N;i++)
	{
		a[i] = i;
		b[i] = i;
		c[i] = a[i] + b[i];
		printf("%-5d, %-5d, %-5d\n",a[i],b[i],c[i]); 
	}


	return 0;
}



