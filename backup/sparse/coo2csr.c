#include <stdio.h>
#include <stdlib.h>
#define nnz 9
#define n_row 4
void main()

{

	// the order!!!!!!!!!!!!!
        int Ai[nnz] = {0,0,1,1,2,2,2,3,3};
	int Aj[nnz] = {0,1,1,2,0,2,3,1,3};
	double A[nnz] = {1.0,7.0,2.0,8.0,5.0,3.0,9.0,6.0,4.0};
	
	int Bp[n_row + 1];
	int Bj[nnz];
	double Bx[nnz];
	int i;
	int cumsum, temp, last, dest, row;

	// the size of Bp is row + 1	
	
 	for(i = 0; i < n_row + 1; i++)
	{
		Bp[i] = 0;	
	}
	
	for(i = 0; i < nnz; i++)
	{
		Bp[ Ai[i] ]++;	
	}

	
	for(i = 0, cumsum = 0; i < n_row; i++)
	{
		temp = Bp[i];
		Bp[i] = cumsum;
		cumsum +=temp;
	} 
	Bp[ n_row + 1] = nnz;

	for(i = 0; i < nnz; i++)
	{
		row = Ai[i];
		dest = Bp[row];
		Bj[dest] = Aj[i];
		Bx[dest] = A[i];
		Bp[row] ++;
	}
	

	for(i = 0, last = 0; i < n_row + 1; i++)
	{
		temp =  Bp[i];
		Bp[i] = last;
		last = temp;
	}
	
	for(i = 0; i < n_row + 1; i++)
	{
		printf("%5d", Bp[i]);	
	}
	printf("\n");

	for(i = 0;i < nnz; i++)
	{
		printf("%5d", Bj[i]);	
	}
	printf("\n");

	for(i = 0;i < nnz; i++)
	{
		printf("%5f	", Bx[i]);	
	}
	printf("\n");


}
