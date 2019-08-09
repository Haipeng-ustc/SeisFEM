/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/*! @file superlu.c
 * \brief a small 5x5 example
 * 
 * <pre>
 * * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 * </pre>
 */
#include "slu_ddefs.h"

int main()
{
    /*
 * Purpose
 * =======
 * 
 * This is the small 5x5 example used in the Sections 2 and 3 of the 
 * Users' Guide to illustrate how to call a SuperLU routine, and the
 * matrix data structures used by SuperLU.
 *
 */
    SuperMatrix A, L, U, B;
    double *a, *rhs;
    double s, u, p, e, r, l;
    int *asub, *xa;
    int *perm_r; /* row permutations from partial pivoting */
    int *perm_c; /* column permutation vector */
    int nrhs, info, i, m, n, nnz, permc_spec;
    FILE *fp_mass_p, *fp_mass_j, *fp_mass_x;
    superlu_options_t options;
    SuperLUStat_t stat;

    /* Initialize matrix A. */
    m = n = 441;
    nnz = 2921;
    if (!(a = doubleMalloc(nnz)))
        ABORT("Malloc fails for a[].");
    if (!(asub = intMalloc(nnz)))
        ABORT("Malloc fails for asub[].");
    if (!(xa = intMalloc(n + 1)))
        ABORT("Malloc fails for xa[].");
    fp_mass_p = fopen("mass_p.dat", "r");
    fp_mass_j = fopen("mass_j.dat", "r");
    fp_mass_x = fopen("mass_x.dat", "r");
  //   csr to csc
     for (i = 0; i <= m; i++)
    {
        fscanf(fp_mass_p, "%d\n", &xa[i]);
    }
    for (i = 0; i < nnz; i++)
    {
        fscanf(fp_mass_j, "%d\n", &asub[i]);
        fscanf(fp_mass_j, "%lf\n", &a[i]);
    }
    fclose(fp_mass_p);
    fclose(fp_mass_j);
    fclose(fp_mass_x);

    /*  a[0] = 1.0; a[1] = 7.0; a[2] = 2.0; a[3] = 8.0; a[4] = 5.0; a[5] = 3.0;
    a[6] = 9.0; a[7] = 6.0; a[8] = 4.0;
    asub[0] = 0; asub[1] = 1; asub[2] = 1; asub[3] = 2;
    asub[4] = 0; asub[5] = 2; asub[6] = 3; asub[7] = 1;
    asub[8] = 3;
    xa[0] = 0; xa[1] = 2; xa[2] = 4; xa[3] = 7; xa[4] = 9;
    */

    /* Create matrix A in the format expected by SuperLU. */
    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

    /* Create right-hand side matrix B. */
    nrhs = 1;
    if (!(rhs = doubleMalloc(m * nrhs)))
        ABORT("Malloc fails for rhs[].");
    for (i = 0; i < m; ++i)
        rhs[i] =  1.0;
    dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);

    if (!(perm_r = intMalloc(m)))
        ABORT("Malloc fails for perm_r[].");
    if (!(perm_c = intMalloc(n)))
        ABORT("Malloc fails for perm_c[].");

    /* Set the default input options. */
    set_default_options(&options);
    options.ColPerm = NATURAL; //COLAMD to speed up

    /* Initialize the statistics variables. */
    StatInit(&stat);

    /* Solve the linear system. */
    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

    //dPrint_CompCol_Matrix("A", &A);
    //dPrint_CompCol_Matrix("U", &U);
    //dPrint_SuperNode_Matrix("L", &L);
    //dPrint_SuperNode_Matrix("L", &L);
    //print_int_vec("\nperm_r", m, perm_r);
    dPrint_Dense_Matrix("B", &B);

    /* De-allocate storage */
    SUPERLU_FREE(rhs);
    SUPERLU_FREE(perm_r);
    SUPERLU_FREE(perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree(&stat);

    return 0;
}
