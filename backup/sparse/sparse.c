#include <stdio.h>
#include <stdlib.h>
#define nnz 9
#define n_row 4

//void csr_matvec(int Ap_size, int Ax_size, int Ap[], int Aj[], double Ax[], int x_size[], double x[], double y[])

// construct this funtion at last


void coo2csr(int nnz, int Ai[], int Aj, double Ax[], int n_row, int Bp[], int Bj[], double Bx[])
{
	
	// Converts from COO (Ai, Aj, Ax) into CSR (Bp, Bj, Bx)
        // Row and column indices are *not* assumed to be ordered.
        // Duplicate entries are carried over to the CSR representation.
	
	
	int i, cumsum, temp, row, dest;
	
	// the size of Bp is n_row + 1
	
	for(i = 0; i < n_row + 1; i++)
	{
		Bp[i] = 0;	// set zero
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
	for(i = 0; i < nnz; i++)
	{
		row = Ai[i];
		dest = Bp[row];
		Bj[dest] = Aj[i];
		Bx[dest] = Ax[i];
		Bp[row] ++;
	}
	for( i = n_row ; i > 0 ; i-- )
	{
		Bp[i] = Bp[ i - 1 ];
	} 
	Bp[0] = 0;
		

}


void main()
{
	
        int Ai[nnz] = {0,0,1,1,2,2,2,3,3};
	int Aj[nnz] = {0,1,1,2,0,2,3,1,3};
	double Ax[nnz] = {1.0,7.0,2.0,8.0,5.0,3.0,9.0,6.0,4.0};
	
	int Bp[n_row + 1];
	int Bj[nnz];
	double Bx[nnz];

	int i, cumsum, temp, row, dest;
	
	// the size of Bp is n_row + 1
	
	for(i = 0; i < n_row + 1; i++)
	{
		Bp[i] = 0;	// set zero
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
	for(i = 0; i < nnz; i++)
	{
		row = Ai[i];
		dest = Bp[row];
		Bj[dest] = Aj[i];
		Bx[dest] = Ax[i];
		Bp[row] ++;
	}
	for( i = n_row ; i > 0 ; i-- )
	{
		Bp[i] = Bp[ i - 1 ];
	} 
	Bp[0] = 0;
	

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
