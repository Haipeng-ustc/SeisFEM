#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#define N 18

//void get_csr_size_(int *, int *, int (*)[], int (*)[], double (*)[], int *);

int main()
{
	int Ai[N] = {0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 2, 2, 2, 1, 1, 1};
	int Aj[N] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 1, 3, 2, 1, 3, 2, 1};
	double Ax[N] = {1.0, -0.5, -0.5, -0.5, 0.5, 0.0, -0.5, 0.0, 0.5, 1.0, -0.5, -0.5, -0.5, 0.5, 0.0, -0.5, 0.0, 0.5};
	int Bp[N] = {0};
	int Bj[N] = {0};
	double Bx[N] = {0.0};

	int nnz, node_num, Bp_size, csr_size;
	int i;
	nnz = N;
	node_num = N;
	csr_size = 0;
	Bp_size = 4 + 1;

	__coo2csr_lib_MOD_coo2csr_canonical(&nnz, &Bp_size, &csr_size, &Ai, &Aj, &Ax, &Bp, &Bj, &Bx);

	printf("\ncrs size is %d\n", csr_size);
	printf("\ncoo format \n");
	for (i = 0; i < N; i++)
	{
		printf("%d %d %f \n", Ai[i], Aj[i], Ax[i]);
	}

	printf("\ncsr format \n");
	for (i = 0; i < Bp_size; i++)
	{
		printf("%d \n", Bp[i]);
	}
	for (i = 0; i < csr_size; i++)
	{
		printf("%d %f \n", Bj[i], Bx[i]);
	}
	printf("Normal End of Excuation");
	//free ( Bj );
	//free ( Bx );

	return 0;
}
