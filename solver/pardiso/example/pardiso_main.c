#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pardiso_unsym.c"
int main()
{
    /* Matrix data. */
    FILE *fp_mass_p, *fp_mass_j, *fp_mass_x, *fp_x;
    int i;
    int n = 185033;
    int nnz = 1291781;
    int *ia;
    int *ja;
    double *a;
    double *rhs;
    double *x;
    ia = (int *)malloc((n + 1) * sizeof(int));
    ja = (int *)malloc(nnz * sizeof(int));
    a = (double *)malloc(nnz * sizeof(double));
    x = (double *)malloc(n * sizeof(double));
    rhs = (double *)malloc(n * sizeof(double));

    /* Set right hand side to i. */
    for (i = 0; i < n; i++)
    {
        rhs[i] = 1.0 * i;
    }

    fp_mass_p = fopen("mass_p.dat", "r");
    fp_mass_j = fopen("mass_j.dat", "r");
    fp_mass_x = fopen("mass_x.dat", "r");
    //   csr to csc
    for (i = 0; i <= n; i++)
    {
        fscanf(fp_mass_p, "%d\n", &ia[i]);
    }
    for (i = 0; i < nnz; i++)
    {
        fscanf(fp_mass_j, "%d\n", &ja[i]);
        fscanf(fp_mass_x, "%lf\n", &a[i]);
    }
    fclose(fp_mass_p);
    fclose(fp_mass_j);
    fclose(fp_mass_x);
    pardiso_unsym(nnz, n, ia, ja, a, rhs, x);
    fp_x = fopen("result.dat", "w");
    for (i = 0; i < n; i++)
    {
        fprintf(fp_x, "%lf\n",x[i]);
    }
    fclose(fp_x);
}
