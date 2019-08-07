void elastic_mpml(char *type, int node_num, int element_num, int element_order, int *element_node, double **node_xy, int nnz, int csr_p_size, int step, double dt,
                  double f0, double t0, double edge_size, double xmin, double xmax, double ymin, double ymax, int pml_nx, int pml_ny, int source_node, char *solver)
{

    /*    stiffness matrix List:
 
        I      Stiffness_TYPE
        -  --------------------------------------
        0  dphi_dx * dphi_dx + dphi_dy * dphi_dy
        1  dphi_dx * dphi_dx
        2  dphi_dy * dphi_dy
        3  dphi_dx * dphi_dy
        4  dphi_dy * dphi_dx
        5  phi     * dphi_dx
        6  phi     * dphi_dy
 */
    /***************************************
     coo(i,j,x) and csr(p,j,x) arrays
     ****************************************/
    int mass_csr_size = 0;                                                                                                      // csr size for mass : j and x
    int stif1_csr_size = 0, stif2_csr_size = 0, stif3_csr_size = 0, stif4_csr_size = 0, stif5_csr_size = 0, stif6_csr_size = 0; // csr size for stif1-6: j and x
    int *mass_coo_i = NULL;
    int *mass_coo_j = NULL;
    int *mass_csr_p = NULL;
    int *mass_csr_j = NULL;
    int *mass_csr_j_temp = NULL;
    int *stif_coo_i = NULL;
    int *stif_coo_j = NULL;
    int *stif_csr_j_temp = NULL;
    int *stif1_csr_p = NULL;
    int *stif2_csr_p = NULL;
    int *stif3_csr_p = NULL;
    int *stif4_csr_p = NULL;
    int *stif5_csr_p = NULL;
    int *stif6_csr_p = NULL;
    int *stif1_csr_j = NULL;
    int *stif2_csr_j = NULL;
    int *stif3_csr_j = NULL;
    int *stif4_csr_j = NULL;
    int *stif5_csr_j = NULL;
    int *stif6_csr_j = NULL;
    double *mass_lump = NULL;
    double *mass_coo_x = NULL;
    double *mass_csr_x = NULL;
    double *mass_csr_x_temp = NULL;
    double *stif_coo_x = NULL;
    double *stif_csr_x_temp = NULL;
    double *stif1_csr_x = NULL;
    double *stif2_csr_x = NULL;
    double *stif3_csr_x = NULL;
    double *stif4_csr_x = NULL;
    double *stif5_csr_x = NULL;
    double *stif6_csr_x = NULL;

    /***************************************
     model density and velocity parameters
     ****************************************/
    double *rho = NULL;
    double *vp = NULL;
    double *vs = NULL;
    double **c = NULL;
    double Angle_force = 135.0;
    double pi = 3.1415926535898;
    /***************************************
             time evolution parameters
     ****************************************/
    int i, it;
    double time;
    double vp_max;
    double point_source;
    double *Energy_u = NULL, *Energy_w = NULL;
    double *U_now = NULL, *W_now = NULL;
    double *U1_now = NULL, *U2_now = NULL, *U3_now = NULL;
    double *W1_now = NULL, *W2_now = NULL, *W3_now = NULL;
    double *U1t_now = NULL, *U2t_now = NULL, *U3t_now = NULL;
    double *W1t_now = NULL, *W2t_now = NULL, *W3t_now = NULL;
    double *U1tt_now = NULL, *U2tt_now = NULL, *U3tt_now = NULL;
    double *W1tt_now = NULL, *W2tt_now = NULL, *W3tt_now = NULL;
    double *U1tt_new = NULL, *U2tt_new = NULL, *U3tt_new = NULL;
    double *W1tt_new = NULL, *W2tt_new = NULL, *W3tt_new = NULL;
    double *Lx1_now = NULL, *Lx2_now = NULL, *Lx3_now = NULL, *Lx4_now = NULL;
    double *Ly1_now = NULL, *Ly2_now = NULL, *Ly3_now = NULL, *Ly4_now = NULL;
    double *rhs_u1 = NULL, *rhs_u2 = NULL, *rhs_u3 = NULL, *rhs_u4 = NULL, *rhs_u5 = NULL, *rhs_u6 = NULL, *rhs_u7 = NULL;
    double *rhs_w1 = NULL, *rhs_w2 = NULL, *rhs_w3 = NULL, *rhs_w4 = NULL, *rhs_w5 = NULL, *rhs_w6 = NULL, *rhs_w7 = NULL;
    double *equ1_1 = NULL, *equ1_2 = NULL, *equ1_3 = NULL, *equ1_4 = NULL, *equ1_5 = NULL;
    double *equ2_1 = NULL, *equ2_2 = NULL, *equ2_3 = NULL, *equ2_4 = NULL;
    double *equ3_1 = NULL, *equ3_2 = NULL, *equ3_3 = NULL, *equ3_4 = NULL, *equ3_5 = NULL;
    double *equ4_1 = NULL, *equ4_2 = NULL;
    double *equ5_1 = NULL, *equ5_2 = NULL;
    double *equ6_1 = NULL, *equ6_2 = NULL;
    double *equ7_1 = NULL, *equ7_2 = NULL;
    double delta = 1.5, alpha = 1.0;
    double tol_abs = 1.0e-08;   // a relative tolerance comparing the current residual to the initial residual.
    double tol_rel = 1.0e-08;   // an absolute tolerance applied to the current residual.
    int itr_max = 100, mr = 50; // the maximum number of outer and inner iterations to take.

    /***************************************
     abosorbing bc, mpml parameters
     ****************************************/
    int use_mpml_xmin = 1;
    int use_mpml_xmax = 1;
    int use_mpml_ymin = 1;
    int use_mpml_ymax = 1; // use_pml_ymax = 0, free surface on ymax.
    double *mpml_dx = NULL;
    double *mpml_dy = NULL;
    double *mpml_dxx = NULL;
    double *mpml_dyy = NULL;
    double *mpml_dxx_pyx = NULL;
    double *mpml_dyy_pxy = NULL;

    /***************************************
                  file pointers
     ****************************************/
    FILE *fp_u, *fp_w, *fp_Energy_u, *fp_Energy_w;

    /***************************************
          allocate the dynamic arrays
     ****************************************/
    mass_lump = (double *)malloc(node_num * sizeof(double));
    mass_coo_i = (int *)malloc(nnz * sizeof(int));
    mass_coo_j = (int *)malloc(nnz * sizeof(int));
    stif_coo_i = (int *)malloc(nnz * sizeof(int));
    stif_coo_j = (int *)malloc(nnz * sizeof(int));
    mass_csr_p = (int *)malloc(csr_p_size * sizeof(int));
    stif1_csr_p = (int *)malloc(csr_p_size * sizeof(int));
    stif2_csr_p = (int *)malloc(csr_p_size * sizeof(int));
    stif3_csr_p = (int *)malloc(csr_p_size * sizeof(int));
    stif4_csr_p = (int *)malloc(csr_p_size * sizeof(int));
    stif5_csr_p = (int *)malloc(csr_p_size * sizeof(int));
    stif6_csr_p = (int *)malloc(csr_p_size * sizeof(int));
    stif_csr_j_temp = (int *)malloc(nnz * sizeof(int));
    mass_csr_j_temp = (int *)malloc(nnz * sizeof(int));
    mass_coo_x = (double *)malloc(nnz * sizeof(double));
    stif_coo_x = (double *)malloc(nnz * sizeof(double));
    mass_csr_x_temp = (double *)malloc(nnz * sizeof(double));
    stif_csr_x_temp = (double *)malloc(nnz * sizeof(double));

    rho = (double *)malloc(node_num * sizeof(double));
    vp = (double *)malloc(node_num * sizeof(double));
    vs = (double *)malloc(node_num * sizeof(double));
    c = (double **)malloc(4 * sizeof(double));
    for (i = 0; i < 4; i++)
        c[i] = malloc(node_num * sizeof(double));
    mpml_dx = (double *)malloc(node_num * sizeof(double));
    mpml_dy = (double *)malloc(node_num * sizeof(double));
    mpml_dxx = (double *)malloc(node_num * sizeof(double));
    mpml_dyy = (double *)malloc(node_num * sizeof(double));
    mpml_dxx_pyx = (double *)malloc(node_num * sizeof(double));
    mpml_dyy_pxy = (double *)malloc(node_num * sizeof(double));
    // time evolution
    Energy_u = (double *)malloc(step * sizeof(double));
    Energy_w = (double *)malloc(step * sizeof(double));
    U_now = (double *)malloc(node_num * sizeof(double));
    W_now = (double *)malloc(node_num * sizeof(double));
    U1_now = (double *)malloc(node_num * sizeof(double));
    U2_now = (double *)malloc(node_num * sizeof(double));
    U3_now = (double *)malloc(node_num * sizeof(double));
    W1_now = (double *)malloc(node_num * sizeof(double));
    W2_now = (double *)malloc(node_num * sizeof(double));
    W3_now = (double *)malloc(node_num * sizeof(double));
    U1t_now = (double *)malloc(node_num * sizeof(double));
    U2t_now = (double *)malloc(node_num * sizeof(double));
    U3t_now = (double *)malloc(node_num * sizeof(double));
    W1t_now = (double *)malloc(node_num * sizeof(double));
    W2t_now = (double *)malloc(node_num * sizeof(double));
    W3t_now = (double *)malloc(node_num * sizeof(double));
    U1tt_now = (double *)malloc(node_num * sizeof(double));
    U2tt_now = (double *)malloc(node_num * sizeof(double));
    U3tt_now = (double *)malloc(node_num * sizeof(double));
    W1tt_now = (double *)malloc(node_num * sizeof(double));
    W2tt_now = (double *)malloc(node_num * sizeof(double));
    W3tt_now = (double *)malloc(node_num * sizeof(double));
    U1tt_new = (double *)malloc(node_num * sizeof(double));
    U2tt_new = (double *)malloc(node_num * sizeof(double));
    U3tt_new = (double *)malloc(node_num * sizeof(double));
    W1tt_new = (double *)malloc(node_num * sizeof(double));
    W2tt_new = (double *)malloc(node_num * sizeof(double));
    W3tt_new = (double *)malloc(node_num * sizeof(double));
    Lx1_now = (double *)malloc(node_num * sizeof(double));
    Lx2_now = (double *)malloc(node_num * sizeof(double));
    Lx3_now = (double *)malloc(node_num * sizeof(double));
    Lx4_now = (double *)malloc(node_num * sizeof(double));
    Ly1_now = (double *)malloc(node_num * sizeof(double));
    Ly2_now = (double *)malloc(node_num * sizeof(double));
    Ly3_now = (double *)malloc(node_num * sizeof(double));
    Ly4_now = (double *)malloc(node_num * sizeof(double));
    equ1_1 = (double *)malloc(node_num * sizeof(double));
    equ1_2 = (double *)malloc(node_num * sizeof(double));
    equ1_3 = (double *)malloc(node_num * sizeof(double));
    equ1_4 = (double *)malloc(node_num * sizeof(double));
    equ1_5 = (double *)malloc(node_num * sizeof(double));
    equ2_1 = (double *)malloc(node_num * sizeof(double));
    equ2_2 = (double *)malloc(node_num * sizeof(double));
    equ2_3 = (double *)malloc(node_num * sizeof(double));
    equ2_4 = (double *)malloc(node_num * sizeof(double));
    equ3_1 = (double *)malloc(node_num * sizeof(double));
    equ3_2 = (double *)malloc(node_num * sizeof(double));
    equ3_3 = (double *)malloc(node_num * sizeof(double));
    equ3_4 = (double *)malloc(node_num * sizeof(double));
    equ3_5 = (double *)malloc(node_num * sizeof(double));
    equ4_1 = (double *)malloc(node_num * sizeof(double));
    equ4_2 = (double *)malloc(node_num * sizeof(double));
    equ5_1 = (double *)malloc(node_num * sizeof(double));
    equ5_2 = (double *)malloc(node_num * sizeof(double));
    equ6_1 = (double *)malloc(node_num * sizeof(double));
    equ6_2 = (double *)malloc(node_num * sizeof(double));
    equ7_1 = (double *)malloc(node_num * sizeof(double));
    equ7_2 = (double *)malloc(node_num * sizeof(double));
    rhs_u1 = (double *)malloc(node_num * sizeof(double));
    rhs_u2 = (double *)malloc(node_num * sizeof(double));
    rhs_u3 = (double *)malloc(node_num * sizeof(double));
    rhs_u4 = (double *)malloc(node_num * sizeof(double));
    rhs_u5 = (double *)malloc(node_num * sizeof(double));
    rhs_u6 = (double *)malloc(node_num * sizeof(double));
    rhs_u7 = (double *)malloc(node_num * sizeof(double));
    rhs_w1 = (double *)malloc(node_num * sizeof(double));
    rhs_w2 = (double *)malloc(node_num * sizeof(double));
    rhs_w3 = (double *)malloc(node_num * sizeof(double));
    rhs_w4 = (double *)malloc(node_num * sizeof(double));
    rhs_w5 = (double *)malloc(node_num * sizeof(double));
    rhs_w6 = (double *)malloc(node_num * sizeof(double));
    rhs_w7 = (double *)malloc(node_num * sizeof(double));

    elastic_model(node_num, element_num, element_order, element_node, node_xy, rho, vp, vs, c);
    vp_max = vp[0];
    for (i = 0; i < node_num; i++)
    {
        if (vp[i] > vp_max)
            vp_max = vp[i];
    }

    abc_mpml(node_num, element_num, element_order, element_node, node_xy, pml_nx, pml_ny, edge_size, xmin, xmax, ymin, ymax, vp_max,
                            use_mpml_xmin, use_mpml_xmax, use_mpml_ymin, use_mpml_ymax, mpml_dx, mpml_dy, mpml_dxx, mpml_dyy, mpml_dxx_pyx, mpml_dyy_pxy);

    fp_u = fopen("u.dat", "w");
    fp_w = fopen("w.dat", "w");
    fp_Energy_u = fopen("Energy_u.dat", "w");
    fp_Energy_w = fopen("Energy_w.dat", "w");

    /********************************************************
        assemble matrices, first in coo then call 
        fortran subroutines to convert coo to csr. 
        do not need use & to get the address of the pointers.
     *********************************************************/
    // mass lump
    mass_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, rho, mass_coo_i, mass_coo_j, mass_coo_x, 1);
    for (i = 0; i < node_num; i++)
        mass_lump[i] = mass_coo_x[i];
    // mass
    mass_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, rho, mass_coo_i, mass_coo_j, mass_coo_x, 0);
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &mass_csr_size, mass_coo_i, mass_coo_j, mass_coo_x, mass_csr_p, mass_csr_j_temp, mass_csr_x_temp);
    mass_csr_j = (int *)malloc(mass_csr_size * sizeof(int));
    mass_csr_x = (double *)malloc(mass_csr_size * sizeof(double));
    for (i = 0; i < mass_csr_size; i++)
    {
        mass_csr_j[i] = mass_csr_j_temp[i];
        mass_csr_x[i] = mass_csr_x_temp[i];
    }
    // stiffness1: dphidx * dphidx
    stif_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, stif_coo_i, stif_coo_j, stif_coo_x, 1);
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &stif1_csr_size, stif_coo_i, stif_coo_j, stif_coo_x, stif1_csr_p, stif_csr_j_temp, stif_csr_x_temp);
    stif1_csr_j = (int *)malloc(stif1_csr_size * sizeof(int));
    stif1_csr_x = (double *)malloc(stif1_csr_size * sizeof(double));
    for (i = 0; i < stif1_csr_size; i++)
    {
        stif1_csr_j[i] = stif_csr_j_temp[i];
        stif1_csr_x[i] = stif_csr_x_temp[i];
    }
    // stiffness2: dphidy * dphidy
    stif_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, stif_coo_i, stif_coo_j, stif_coo_x, 2);
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &stif2_csr_size, stif_coo_i, stif_coo_j, stif_coo_x, stif2_csr_p, stif_csr_j_temp, stif_csr_x_temp);
    stif2_csr_j = (int *)malloc(stif2_csr_size * sizeof(int));
    stif2_csr_x = (double *)malloc(stif2_csr_size * sizeof(double));
    for (i = 0; i < stif2_csr_size; i++)
    {
        stif2_csr_j[i] = stif_csr_j_temp[i];
        stif2_csr_x[i] = stif_csr_x_temp[i];
    }
    // stiffness3: dphidx * dphidy
    stif_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, stif_coo_i, stif_coo_j, stif_coo_x, 3);
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &stif3_csr_size, stif_coo_i, stif_coo_j, stif_coo_x, stif3_csr_p, stif_csr_j_temp, stif_csr_x_temp);
    stif3_csr_j = (int *)malloc(stif3_csr_size * sizeof(int));
    stif3_csr_x = (double *)malloc(stif3_csr_size * sizeof(double));
    for (i = 0; i < stif3_csr_size; i++)
    {
        stif3_csr_j[i] = stif_csr_j_temp[i];
        stif3_csr_x[i] = stif_csr_x_temp[i];
    }
    // stiffness4: dphidy * dphidx
    stif_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, stif_coo_i, stif_coo_j, stif_coo_x, 4);
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &stif4_csr_size, stif_coo_i, stif_coo_j, stif_coo_x, stif4_csr_p, stif_csr_j_temp, stif_csr_x_temp);
    stif4_csr_j = (int *)malloc(stif4_csr_size * sizeof(int));
    stif4_csr_x = (double *)malloc(stif4_csr_size * sizeof(double));
    for (i = 0; i < stif4_csr_size; i++)
    {
        stif4_csr_j[i] = stif_csr_j_temp[i];
        stif4_csr_x[i] = stif_csr_x_temp[i];
    }
    // stiffness5: phi * dphidx
    stif_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, stif_coo_i, stif_coo_j, stif_coo_x, 5);
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &stif5_csr_size, stif_coo_i, stif_coo_j, stif_coo_x, stif5_csr_p, stif_csr_j_temp, stif_csr_x_temp);
    stif5_csr_j = (int *)malloc(stif5_csr_size * sizeof(int));
    stif5_csr_x = (double *)malloc(stif5_csr_size * sizeof(double));
    for (i = 0; i < stif5_csr_size; i++)
    {
        stif5_csr_j[i] = stif_csr_j_temp[i];
        stif5_csr_x[i] = stif_csr_x_temp[i];
    }
    // stiffness6: phi * dphidy
    stif_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, stif_coo_i, stif_coo_j, stif_coo_x, 6);
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &stif6_csr_size, stif_coo_i, stif_coo_j, stif_coo_x, stif6_csr_p, stif_csr_j_temp, stif_csr_x_temp);
    stif6_csr_j = (int *)malloc(stif6_csr_size * sizeof(int));
    stif6_csr_x = (double *)malloc(stif6_csr_size * sizeof(double));
    for (i = 0; i < stif6_csr_size; i++)
    {
        stif6_csr_j[i] = stif_csr_j_temp[i];
        stif6_csr_x[i] = stif_csr_x_temp[i];
    }
    
    // write first two step values: u_old[node_num], u_now[node_num], energy[0], energy[1]
    for (i = 0; i < node_num; i++)
    {
        U_now[i] = 0.0;
        W_now[i] = 0.0;
        U1_now[i] = 0.0;
        U2_now[i] = 0.0;
        U3_now[i] = 0.0;
        W1_now[i] = 0.0;
        W2_now[i] = 0.0;
        W3_now[i] = 0.0;
        U1t_now[i] = 0.0;
        U2t_now[i] = 0.0;
        U3t_now[i] = 0.0;
        W1t_now[i] = 0.0;
        W2t_now[i] = 0.0;
        W3t_now[i] = 0.0;
        U1tt_now[i] = 0.0;
        U2tt_now[i] = 0.0;
        U3tt_now[i] = 0.0;
        W1tt_now[i] = 0.0;
        W2tt_now[i] = 0.0;
        W3tt_now[i] = 0.0;
        U1tt_new[i] = 0.0;
        U2tt_new[i] = 0.0;
        U3tt_new[i] = 0.0;
        W1tt_new[i] = 0.0;
        W2tt_new[i] = 0.0;
        W3tt_new[i] = 0.0;
        Lx1_now[i] = 0.0;
        Lx2_now[i] = 0.0;
        Lx3_now[i] = 0.0;
        Lx4_now[i] = 0.0;
        Ly1_now[i] = 0.0;
        Ly2_now[i] = 0.0;
        Ly3_now[i] = 0.0;
        Ly4_now[i] = 0.0;

        fprintf(fp_u, "%f	", U_now[i]);
        fprintf(fp_w, "%f	", W_now[i]);
    }
    fprintf(fp_u, "\n");
    fprintf(fp_w, "\n");

    Energy_u[0] = 0.0;
    Energy_u[1] = 0.0;
    Energy_w[0] = 0.0;
    Energy_w[1] = 0.0;

    fprintf(fp_Energy_u, "%f\n%f\n", Energy_u[0], Energy_u[1]);
    fprintf(fp_Energy_w, "%f\n%f\n", Energy_w[0], Energy_w[1]);

    // begin iteration: from 0 to step-1, time = (step + 1) * dt
    printf("\nTime iteration begin:\n");
    for (it = 2; it < step; it++)
    {
        time = (it + 1) * dt;
        if ((it + 1) % 100 == 0)
            printf("\n Iteration step: %-d, time: %-f s\n ", it + 1, time);
        Energy_u[it] = 0.0;
        Energy_w[it] = 0.0;
        point_source = seismic_source(f0, t0, 1.0e10, time);

        /********************************************************************************************************************************************
         Equation u1:
          mass * U1tt_new = - c11 * dphidx * dphidx * U_now - 2.0 * mpml_dx * phi * phi * U1t_now - mpml_dx * mpml_dx * phi * phi * U1_now 
                            + phi * phi * Lx1_now + phi * phi * Lx2_now + Source_x
        *********************************************************************************************************************************************/
        csr_matvec(csr_p_size, stif1_csr_p, stif1_csr_j, stif1_csr_x, U_now, equ1_1); // dphidx * dphidx * U_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, U1t_now, equ1_2);  // phi    * phi    * U1t_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, U1_now, equ1_3);   // phi    * phi    * U1_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, Lx1_now, equ1_4);  // phi    * phi    * Lx1_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, Lx2_now, equ1_5);  // phi    * phi    * Lx2_now
        /********************************************************************************************************************************************
         Equation u2:
          mass * U2tt_new = - c13 * dphidx * dphidy * W_now - c44 * dphidy * dphidx * W_now - mpml_dx * phi * phi * U2t_now 
                            - mpml_dy * phi * phi * U2t_now - mpml_dx * mpml_dy * phi * phi * U2_now
        *********************************************************************************************************************************************/
        csr_matvec(csr_p_size, stif3_csr_p, stif3_csr_j, stif3_csr_x, W_now, equ2_1); // dphidx * dphidy * W_now
        csr_matvec(csr_p_size, stif4_csr_p, stif4_csr_j, stif4_csr_x, W_now, equ2_2); // dphidey * dphidx * W_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, U2t_now, equ2_3);  // phi    * phi    * U2t_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, U2_now, equ2_4);   // phi    * phi    * U2_now
        /********************************************************************************************************************************************
         Equation u3:
          mass * U3tt_new = - c44 * dphidy * dphidy * U_now - 2.0 * mpml_dy * phi * phi * U3t_now - mpml_dy * mpml_dy * phi * phi * U3_now 
                            + phi * phi * Lx3_now + phi * phi * Lx4_now
        *********************************************************************************************************************************************/
        csr_matvec(csr_p_size, stif2_csr_p, stif2_csr_j, stif2_csr_x, U_now, equ3_1); // dphidy * dphidy * U_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, U3t_now, equ3_2);  // phi    * phi    * U3t_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, U3_now, equ3_3);   // phi    * phi    * U3_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, Lx3_now, equ3_4);  // phi    * phi    * Lx3_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, Lx4_now, equ3_5);  // phi    * phi    * Lx4_now
        /********************************************************************************************************************************************
         Equation u4:
          mass * Lx1_new = - dt * c11 * mpml_dxx * phi * dphidx * U_now - dt * mpml_dx * phi * phi * Lx1_now + phi * phi * Lx1_now
        *********************************************************************************************************************************************/
        csr_matvec(csr_p_size, stif5_csr_p, stif5_csr_j, stif5_csr_x, U_now, equ4_1); // phi    * dphidx * U_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, Lx1_now, equ4_2);  // phi    * phi    * Lx1_now
        /********************************************************************************************************************************************
         Equation u5:
          mass * Lx2_new = - dt * c44 * mpml_dyy_pxy * phi * dphidx * W_now - dt * mpml_dy * phi * phi * Lx2_now + phi * phi * Lx2_now
        *********************************************************************************************************************************************/
        csr_matvec(csr_p_size, stif5_csr_p, stif5_csr_j, stif5_csr_x, W_now, equ5_1); // phi    * dphidx * W_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, Lx2_now, equ5_2);  // phi    * phi    * Lx2_now
        /********************************************************************************************************************************************
         Equation u6:
          mass * Lx3_new = - dt * c13 * mpml_dxx_pyx * phi * dphidy * W_now - dt * mpml_dx * phi * phi * Lx3_now + phi * phi * Lx3_now
        *********************************************************************************************************************************************/
        csr_matvec(csr_p_size, stif6_csr_p, stif6_csr_j, stif6_csr_x, W_now, equ6_1); // phi    * dphidy * W_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, Lx3_now, equ6_2);  // phi    * phi    * Lx3_now
        /********************************************************************************************************************************************
         Equation u7:
         mass * Lx4_new = - dt * c44 * mpml_dyy * phi * dphidy * U_now - dt * mpml_dy * phi * phi * Lx4_now + phi * phi * Lx4_now
        *********************************************************************************************************************************************/
        csr_matvec(csr_p_size, stif6_csr_p, stif6_csr_j, stif6_csr_x, U_now, equ7_1); // phi    * dphidy * U_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, Lx4_now, equ7_2);  // phi    * phi    * Lx4_now

        #pragma omp parallel for private(i)
        for (i = 0; i < node_num; i++)
        {
            rhs_u1[i] = -c[0][i] * equ1_1[i] - 2.0 * mpml_dx[i] * equ1_2[i] - mpml_dx[i] * mpml_dx[i] * equ1_3[i] + equ1_4[i] + equ1_5[i] + (i == source_node) * point_source * sin(Angle_force * pi / 180.0);
            rhs_u2[i] = -c[1][i] * equ2_1[i] - c[3][i] * equ2_2[i] - mpml_dx[i] * equ2_3[i] - mpml_dy[i] * equ2_3[i] - mpml_dx[i] * mpml_dy[i] * equ2_4[i];
            rhs_u3[i] = -c[3][i] * equ3_1[i] - 2.0 * mpml_dy[i] * equ3_2[i] - mpml_dy[i] * mpml_dy[i] * equ3_3[i] + equ3_4[i] + equ3_5[i];
            rhs_u4[i] = -dt * c[0][i] * mpml_dxx[i] * equ4_1[i] - dt * mpml_dx[i] * equ4_2[i] + equ4_2[i];
            rhs_u5[i] = -dt * c[3][i] * mpml_dyy_pxy[i] * equ5_1[i] - dt * mpml_dy[i] * equ5_2[i] + equ5_2[i];
            rhs_u6[i] = -dt * c[1][i] * mpml_dxx_pyx[i] * equ6_1[i] - dt * mpml_dx[i] * equ6_2[i] + equ6_2[i];
            rhs_u7[i] = -dt * c[3][i] * mpml_dyy[i] * equ7_1[i] - dt * mpml_dy[i] * equ7_2[i] + equ7_2[i];
        }

        /********************************************************************************************************************************************
         Equation w1:
          mass * W1tt_new = - c44 * dphidx * dphidx * W_now - 2.0 * mpml_dx * phi * phi * W1t_now - mpml_dx * mpml_dx * phi * phi * W1_now 
                            + phi * phi * Ly1_now + phi * phi * Ly2_now + Source_y
        *********************************************************************************************************************************************/
        csr_matvec(csr_p_size, stif1_csr_p, stif1_csr_j, stif1_csr_x, W_now, equ1_1); // dphidx * dphidx * W_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, W1t_now, equ1_2);  // phi    * phi    * W1t_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, W1_now, equ1_3);   // phi    * phi    * W1_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, Ly1_now, equ1_4);  // phi    * phi    * Ly1_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, Ly2_now, equ1_5);  // phi    * phi    * Ly2_now
        /********************************************************************************************************************************************
         Equation w2:
          mass * W2tt_new = - c44 * dphidx * dphidy * U_now - c13 * dphidy * dphidx * U_now - mpml_dx * phi * phi * W2t_now 
                            - mpml_dy * phi * phi * W2t_now - mpml_dx * mpml_dy * phi * phi * W2_now
        *********************************************************************************************************************************************/
        csr_matvec(csr_p_size, stif3_csr_p, stif3_csr_j, stif3_csr_x, U_now, equ2_1); // dphidx * dphidy * U_now
        csr_matvec(csr_p_size, stif4_csr_p, stif4_csr_j, stif4_csr_x, U_now, equ2_2); // dphidy * dphidx * U_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, W2t_now, equ2_3);  // phi    * phi    * W2t_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, W2_now, equ2_4);   // phi    * phi    * W2_now
        /********************************************************************************************************************************************
         Equation w3:
          mass * W3tt_new = - c33 * dphidy * dphidy * W_now - 2.0 * mpml_dy * phi * phi * W3t_now - mpml_dy * mpml_dy * phi * phi * W3_now 
                            + phi * phi * Ly3_now + phi * phi * Ly4_now
        *********************************************************************************************************************************************/
        csr_matvec(csr_p_size, stif2_csr_p, stif2_csr_j, stif2_csr_x, W_now, equ3_1); // dphidy * dphidy * W_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, W3t_now, equ3_2);  // phi    * phi    * W3t_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, W3_now, equ3_3);   // phi    * phi    * W3_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, Ly3_now, equ3_4);  // phi    * phi    * Ly3_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, Ly4_now, equ3_5);  // phi    * phi    * Ly4_now
        /********************************************************************************************************************************************
         Equation 4:
          mass * Ly1_new = - dt * c44 * mpml_dxx * phi * dphidx * W_now - dt * mpml_dx * phi * phi * Ly1_now + phi * phi * Ly1_now
        *********************************************************************************************************************************************/
        csr_matvec(csr_p_size, stif5_csr_p, stif5_csr_j, stif5_csr_x, W_now, equ4_1); // phi    * dphidx * W_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, Ly1_now, equ4_2);  // phi    * phi    * Ly1_now
        /********************************************************************************************************************************************
         Equation w5:
          mass * Ly2_new = - dt * c13 * mpml_dyy_pxy * phi * dphidx * U_now - dt * mpml_dy * phi * phi * Ly2_now + phi * phi * Ly2_now
        *********************************************************************************************************************************************/
        csr_matvec(csr_p_size, stif5_csr_p, stif5_csr_j, stif5_csr_x, U_now, equ5_1); // phi    * dphidx * U_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, Ly2_now, equ5_2);  // phi    * phi    * Ly2_now
        /********************************************************************************************************************************************
         Equation w6:
          mass * Ly3_new = - dt * c44 * mpml_dxx_pyx * phi * dphidy * U_now - dt * mpml_dx * phi * phi * Ly3_now + phi * phi * Ly3_now
        *********************************************************************************************************************************************/
        csr_matvec(csr_p_size, stif6_csr_p, stif6_csr_j, stif6_csr_x, U_now, equ6_1); // phi    * dphidy * U_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, Ly3_now, equ6_2);  // phi    * phi    * Ly3_now
        /********************************************************************************************************************************************
         Equation w7:
          mass * Ly4_new = - dt * c33 * mpml_dyy * phi * dphidy * W_now - dt * mpml_dy * phi * phi * Ly4_now + phi * phi * Ly4_now
        *********************************************************************************************************************************************/
        csr_matvec(csr_p_size, stif6_csr_p, stif6_csr_j, stif6_csr_x, W_now, equ7_1); // phi    * dphidy * W_now
        csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, Ly4_now, equ7_2);  // phi    * phi    * Ly4_now

        #pragma omp parallel for private(i)
        for (i = 0; i < node_num; i++)
        {
            rhs_w1[i] = -c[3][i] * equ1_1[i] - 2.0 * mpml_dx[i] * equ1_2[i] - mpml_dx[i] * mpml_dx[i] * equ1_3[i] + equ1_4[i] + equ1_5[i] + (i == source_node) * point_source * cos(Angle_force * pi / 180.0);
            rhs_w2[i] = -c[3][i] * equ2_1[i] - c[1][i] * equ2_2[i] - mpml_dx[i] * equ2_3[i] - mpml_dy[i] * equ2_3[i] - mpml_dx[i] * mpml_dy[i] * equ2_4[i];
            rhs_w3[i] = -c[2][i] * equ3_1[i] - 2.0 * mpml_dy[i] * equ3_2[i] - mpml_dy[i] * mpml_dy[i] * equ3_3[i] + equ3_4[i] + equ3_5[i];
            rhs_w4[i] = -dt * c[3][i] * mpml_dxx[i] * equ4_1[i] - dt * mpml_dx[i] * equ4_2[i] + equ4_2[i];
            rhs_w5[i] = -dt * c[1][i] * mpml_dyy_pxy[i] * equ5_1[i] - dt * mpml_dy[i] * equ5_2[i] + equ5_2[i];
            rhs_w6[i] = -dt * c[3][i] * mpml_dxx_pyx[i] * equ6_1[i] - dt * mpml_dx[i] * equ6_2[i] + equ6_2[i];
            rhs_w7[i] = -dt * c[2][i] * mpml_dyy[i] * equ7_1[i] - dt * mpml_dy[i] * equ7_2[i] + equ7_2[i];
        }

        /***********************************
                 solve liner system
        ************************************/
        if (strcmp(solver, "masslump") == 0)
        {
            #pragma omp parallel for private(i)
            for (i = 0; i < node_num; i++)
            {
                if (mass_lump[i] != 0.0)
                {
                    U1tt_new[i] = rhs_u1[i] / (1.0 * mass_lump[i]);
                    U2tt_new[i] = rhs_u2[i] / (1.0 * mass_lump[i]);
                    U3tt_new[i] = rhs_u3[i] / (1.0 * mass_lump[i]);
                    Lx1_now[i] = rhs_u4[i] / (1.0 * mass_lump[i]);
                    Lx2_now[i] = rhs_u5[i] / (1.0 * mass_lump[i]);
                    Lx3_now[i] = rhs_u6[i] / (1.0 * mass_lump[i]);
                    Lx4_now[i] = rhs_u7[i] / (1.0 * mass_lump[i]);
                    W1tt_new[i] = rhs_w1[i] / (1.0 * mass_lump[i]);
                    W2tt_new[i] = rhs_w2[i] / (1.0 * mass_lump[i]);
                    W3tt_new[i] = rhs_w3[i] / (1.0 * mass_lump[i]);
                    Ly1_now[i] = rhs_w4[i] / (1.0 * mass_lump[i]);
                    Ly2_now[i] = rhs_w5[i] / (1.0 * mass_lump[i]);
                    Ly3_now[i] = rhs_w6[i] / (1.0 * mass_lump[i]);
                    Ly4_now[i] = rhs_w7[i] / (1.0 * mass_lump[i]);
                }
                else
                    printf("ELASTIC_MPML - Fatal Error: zero value in the mass_csr_lumped"); // in case zero value
            }
        }
        else if (strcmp(solver, "mgmres") == 0)
        {
            pmgmres_ilu_cr(node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, U1tt_new, rhs_u1, itr_max, mr, tol_abs, tol_rel);
            pmgmres_ilu_cr(node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, U2tt_new, rhs_u2, itr_max, mr, tol_abs, tol_rel);
            pmgmres_ilu_cr(node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, U3tt_new, rhs_u3, itr_max, mr, tol_abs, tol_rel);
            pmgmres_ilu_cr(node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, Lx1_now, rhs_u4, itr_max, mr, tol_abs, tol_rel);
            pmgmres_ilu_cr(node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, Lx2_now, rhs_u5, itr_max, mr, tol_abs, tol_rel);
            pmgmres_ilu_cr(node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, Lx3_now, rhs_u6, itr_max, mr, tol_abs, tol_rel);
            pmgmres_ilu_cr(node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, Lx4_now, rhs_u7, itr_max, mr, tol_abs, tol_rel);
            pmgmres_ilu_cr(node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, W1tt_new, rhs_w1, itr_max, mr, tol_abs, tol_rel);
            pmgmres_ilu_cr(node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, W2tt_new, rhs_w2, itr_max, mr, tol_abs, tol_rel);
            pmgmres_ilu_cr(node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, W3tt_new, rhs_w3, itr_max, mr, tol_abs, tol_rel);
            pmgmres_ilu_cr(node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, Ly1_now, rhs_w4, itr_max, mr, tol_abs, tol_rel);
            pmgmres_ilu_cr(node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, Ly2_now, rhs_w5, itr_max, mr, tol_abs, tol_rel);
            pmgmres_ilu_cr(node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, Ly3_now, rhs_w6, itr_max, mr, tol_abs, tol_rel);
            pmgmres_ilu_cr(node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, Ly4_now, rhs_w7, itr_max, mr, tol_abs, tol_rel);
        }
        else
        {
            fprintf(stderr, "\n");
            fprintf(stderr, "ELASTIC_MPML- Fatal error!\n");
            fprintf(stderr, "  Solver type is not set = \"%s\".\n", solver);
            exit(1);
        }
        #pragma omp parallel for private(i)
        for (i = 0; i < node_num; i++)
        {
            U1_now[i] = U1_now[i] + U1t_now[i] * dt + ((0.5 - alpha) * U1tt_now[i] + alpha * U1tt_new[i]) * dt * dt;
            U2_now[i] = U2_now[i] + U2t_now[i] * dt + ((0.5 - alpha) * U2tt_now[i] + alpha * U2tt_new[i]) * dt * dt;
            U3_now[i] = U3_now[i] + U3t_now[i] * dt + ((0.5 - alpha) * U3tt_now[i] + alpha * U3tt_new[i]) * dt * dt;
            W1_now[i] = W1_now[i] + W1t_now[i] * dt + ((0.5 - alpha) * W1tt_now[i] + alpha * W1tt_new[i]) * dt * dt;
            W2_now[i] = W2_now[i] + W2t_now[i] * dt + ((0.5 - alpha) * W2tt_now[i] + alpha * W2tt_new[i]) * dt * dt;
            W3_now[i] = W3_now[i] + W3t_now[i] * dt + ((0.5 - alpha) * W3tt_now[i] + alpha * W3tt_new[i]) * dt * dt;
            U1t_now[i] = U1t_now[i] + ((1 - delta) * U1tt_now[i] + delta * U1tt_new[i]) * dt;
            U2t_now[i] = U2t_now[i] + ((1 - delta) * U2tt_now[i] + delta * U2tt_new[i]) * dt;
            U3t_now[i] = U3t_now[i] + ((1 - delta) * U3tt_now[i] + delta * U3tt_new[i]) * dt;
            W1t_now[i] = W1t_now[i] + ((1 - delta) * W1tt_now[i] + delta * W1tt_new[i]) * dt;
            W2t_now[i] = W2t_now[i] + ((1 - delta) * W2tt_now[i] + delta * W2tt_new[i]) * dt;
            W3t_now[i] = W3t_now[i] + ((1 - delta) * W3tt_now[i] + delta * W3tt_new[i]) * dt;
            U1tt_now[i] = U1tt_new[i];
            U2tt_now[i] = U2tt_new[i];
            U3tt_now[i] = U3tt_new[i];
            W1tt_now[i] = W1tt_new[i];
            W2tt_now[i] = W2tt_new[i];
            W3tt_now[i] = W3tt_new[i];
            U_now[i] = U1_now[i] + U2_now[i] + U3_now[i];
            W_now[i] = W1_now[i] + W2_now[i] + W3_now[i];
            Energy_u[it] += U_now[i] * U_now[i];
            Energy_w[it] += W_now[i] * W_now[i];
            if (Energy_u[it] > 10e6 || Energy_w[it] > 10e6)
            {
                fprintf(stderr, "\n");
                fprintf(stderr, "ELASTIC_MPML - Fatal error!\n");
                fprintf(stderr, "Energy exceeds maximum value!\n");
                exit(1);
            }
        }
        if ((it + 1) % 50 == 0)
        {
            for (i = 0; i < node_num; i++)
            {
                fprintf(fp_u, "%f	", U_now[i]);
                fprintf(fp_w, "%f	", W_now[i]);
            }
            fprintf(fp_u, "\n");
            fprintf(fp_w, "\n");
        }
        fprintf(fp_Energy_u, "%f\n", Energy_u[it]);
        fprintf(fp_Energy_w, "%f\n", Energy_w[it]);
    }

    printf("\nTime iteration end!\n");
    free(mass_lump);
    free(mass_coo_i);
    free(mass_coo_j);
    free(mass_coo_x);
    free(stif_coo_i);
    free(stif_coo_j);
    free(stif_coo_x);
    free(mass_csr_p);
    free(mass_csr_j);
    free(mass_csr_x);
    free(stif1_csr_p);
    free(stif1_csr_j);
    free(stif1_csr_x);
    free(stif2_csr_p);
    free(stif2_csr_j);
    free(stif2_csr_x);
    free(stif3_csr_p);
    free(stif3_csr_j);
    free(stif3_csr_x);
    free(stif4_csr_p);
    free(stif4_csr_j);
    free(stif4_csr_x);
    free(stif5_csr_p);
    free(stif5_csr_j);
    free(stif5_csr_x);
    free(stif6_csr_p);
    free(stif6_csr_j);
    free(stif6_csr_x);
    free(mass_csr_j_temp);
    free(mass_csr_x_temp);
    free(stif_csr_j_temp);
    free(stif_csr_x_temp);
    free(rho);
    free(vp);
    free(vs);
    for (i = 0; i < 4; i++)
        free(c[i]);
    free(c);
    free(mpml_dx);
    free(mpml_dy);
    free(mpml_dxx);
    free(mpml_dyy);
    free(mpml_dxx_pyx);
    free(mpml_dyy_pxy);
    free(Energy_u);
    free(Energy_w);
    free(U_now);
    free(W_now);
    free(U1_now);
    free(U2_now);
    free(U3_now);
    free(W1_now);
    free(W2_now);
    free(W3_now);
    free(U1t_now);
    free(U2t_now);
    free(U3t_now);
    free(W1t_now);
    free(W2t_now);
    free(W3t_now);
    free(U1tt_now);
    free(U2tt_now);
    free(U3tt_now);
    free(W1tt_now);
    free(W2tt_now);
    free(W3tt_now);
    free(U1tt_new);
    free(U2tt_new);
    free(U3tt_new);
    free(W1tt_new);
    free(W2tt_new);
    free(W3tt_new);
    free(Lx1_now);
    free(Lx2_now);
    free(Lx3_now);
    free(Lx4_now);
    free(Ly1_now);
    free(Ly2_now);
    free(Ly3_now);
    free(Ly4_now);
    free(equ1_1);
    free(equ1_2);
    free(equ1_3);
    free(equ1_4);
    free(equ1_5);
    free(equ2_1);
    free(equ2_2);
    free(equ2_3);
    free(equ2_4);
    free(equ3_1);
    free(equ3_2);
    free(equ3_3);
    free(equ3_4);
    free(equ3_5);
    free(equ4_1);
    free(equ4_2);
    free(equ5_1);
    free(equ5_2);
    free(equ6_1);
    free(equ6_2);
    free(equ7_1);
    free(equ7_2);
    free(rhs_u1);
    free(rhs_u2);
    free(rhs_u3);
    free(rhs_u4);
    free(rhs_u5);
    free(rhs_u6);
    free(rhs_u7);
    free(rhs_w1);
    free(rhs_w2);
    free(rhs_w3);
    free(rhs_w4);
    free(rhs_w5);
    free(rhs_w6);
    free(rhs_w7);
    fclose(fp_u);
    fclose(fp_w);
    fclose(fp_Energy_u);
    fclose(fp_Energy_w);
    
    printf("\n Elastic_mpml Normal End!\n");
}
