#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#define N 1000000
int main()
{

    int i;
    double *a, *b, *c;
    a = (double *)malloc(N * sizeof(double));
    b = (double *)malloc(N * sizeof(double));
    c = (double *)malloc(N * sizeof(double));

    double start_time, run_time;

    for (i = 0; i < N; i++)
    {
        a[i] = 1.0 * i;
        b[i] = 2.0 * i;
    }

    start_time = omp_get_wtime();
    #pragma omp parallel for private(i)
    for (i = 0; i < N; i++)
    {
        c[i] = a[i] + b[i] * a[i] / (a[i] + b[i] + 1.0);
    }
    run_time = omp_get_wtime() - start_time;
    printf("%f\n", run_time);

    free(a);
    free(b);
    free(c);
    return 0;
}
