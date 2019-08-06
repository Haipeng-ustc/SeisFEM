#include <stdio.h>
#include "mgmres.c"

int main()
{
    int i;
    int n = 4;
    int nz_num = 9;
    int ia[5] = {0,2,4,7,9};
    int ja[9] = {0,1,1,2,0,2,3,1,3};
    double a[9] ={1.0,7.0,2.0,8.0,5.0,3.0,9.0,6.0,4.0};
    double x[4] = {0};
    double rhs[4] = {1.0,2.0,3.0,4.0};
    int itr_max, mr;
    int test;
    double tol_abs = 1.0e-08;
    double tol_rel = 1.0e-08;
    int rhs_zero_flag = 0;
    itr_max = 50;
    mr = 40;
    pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x, rhs, itr_max, mr, tol_abs, tol_rel );

    for (i=0; i<n; i++) {
        printf("%.10lf\n",x[i]);
    }
    
    
}
