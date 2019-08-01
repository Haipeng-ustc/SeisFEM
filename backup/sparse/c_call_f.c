#include <stdio.h>
#include <stdlib.h>
#define nnz 9
#define n_row 4

//void myadd_( int *, int (*)[], int (*)[], int (*)[]);

int main()
{
        // int ia[N] = {0,0,1,1,2,2,3,3,3};
	// int ja[N] = {0,1,1,2,4,2,3,1,1};
	// double a[N] = {1.0,7.0,2.0,8.0,5.0,3.0,9.0,6.0,4.0};
	
	int a[4] = {1,2,3,4};
	int b[4] = {1,2,3,4};
	int c[4] = {0,0,0,0};
	int i;
	int n = 4;	
	for(i=0;i<n;i++) printf("%d, %d, %d\n",a[i],b[i],c[i]);

	myadd_(&n, &a, &b, &c);



	printf("after add\n");
	for(i=0;i<n;i++) printf("%d, %d, %d\n",a[i],b[i],c[i]);

}



