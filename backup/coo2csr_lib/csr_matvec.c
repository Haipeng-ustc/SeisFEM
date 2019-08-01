void csr_matvec(int Bp_size, int *Bp, int *Bj, double *Bx, double *x, double *y)

{

	int i, k;
	double t;
	
 	for(i = 0; i < Bp_size - 1; i++) // Bp_size - 1 is node_num, which is the size of x and y
	{
		t = 0.0;
		for( k = Bp[i]; k <= Bp[i+1]-1; k++)
		{
			t = t + Bx[k] * x[ Bj[k] ];
		}
		y[i] = t;
	}


}
