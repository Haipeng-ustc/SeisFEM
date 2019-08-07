#include "../solver/mgmres/mgmres.c"

void elastic_newmark_mpml(char *type, int node_num, int element_num, int element_order, int *element_node, double **node_xy, int nnz, int csr_p_size, int step, \
     double dt, double f0, double t0, double edge_size, double xmin, double xmax, double ymin, double ymax, int pml_nx, int pml_ny, int source_node, char *solver)
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
    int mass_csr_size  = 0; // csr size for mass : j and x
    int stif1_csr_size = 0; // csr size for stif1: j and x
    int stif2_csr_size = 0; // csr size for stif2: j and x
    int stif3_csr_size = 0; // csr size for stif3: j and x
    int stif4_csr_size = 0; // csr size for stif4: j and x
    int stif5_csr_size = 0; // csr size for stif5: j and x
    int stif6_csr_size = 0; // csr size for stif6: j and x
    int *mass_coo_i = NULL;
    int *mass_coo_j = NULL;
    int *mass_csr_p = NULL;
    int *mass_csr_j = NULL;
    int *mass_csr_j_temp = NULL;
    double *mass_coo_x = NULL;
    double *mass_csr_x = NULL;
    double *mass_csr_x_temp = NULL;
    double *mass_lump = NULL;
    int *stif_coo_i = NULL;
    int *stif_coo_j = NULL;
    int *stif_csr_j_temp = NULL;
    double *stif_coo_x = NULL;
    double *stif_csr_x_temp = NULL;
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
    double *vp  = NULL;
    double *vs  = NULL;
    double **c  = NULL;
    double Angle_force = 135.0;
    double pi = 3.1415926535898;
    /***************************************
             time evolution parameters
     ****************************************/
	int i, it;
	double time;
    double vp_max;
	double point_source;
    double *U_now    = NULL;  double *W_now    = NULL;  double *energy   = NULL;
    double *U1_now   = NULL;  double *U2_now   = NULL;  double *U3_now   = NULL;  double *W1_now   = NULL;  double *W2_now   = NULL;  double *W3_now   = NULL;
    double *U1t_now  = NULL;  double *U2t_now  = NULL;  double *U3t_now  = NULL;  double *W1t_now  = NULL;  double *W2t_now  = NULL;  double *W3t_now  = NULL;
    double *U1tt_now = NULL;  double *U2tt_now = NULL;  double *U3tt_now = NULL;  double *W1tt_now = NULL;  double *W2tt_now = NULL;  double *W3tt_now = NULL; 
    double *Lx1_now  = NULL;  double *Lx2_now  = NULL;  double *Lx3_now  = NULL;  double *Lx4_now  = NULL; 
    double *Ly1_now  = NULL;  double *Ly2_now  = NULL;  double *Ly3_now  = NULL;  double *Ly4_now  = NULL;  
	double *rhs_u1   = NULL;  double *rhs_u2   = NULL;  double *rhs_u3   = NULL;  double *rhs_u4   = NULL;  double *rhs_u5   = NULL;  
    double *rhs_u6   = NULL;  double *rhs_u7   = NULL;  double *rhs_w1   = NULL;  double *rhs_w2   = NULL;  double *rhs_w3   = NULL; 
    double *rhs_w4   = NULL;  double *rhs_w5   = NULL;  double *rhs_w6   = NULL;  double *rhs_w7   = NULL;
	double *y1       = NULL;  double *y2       = NULL;  double *y3       = NULL;  double *y4       = NULL;  double *y5       = NULL;

    double tol_abs = 1.0e-08; // a relative tolerance comparing the current residual to the initial residual.
    double tol_rel = 1.0e-08; // an absolute tolerance applied to the current residual.
    int    itr_max = 100;     // the maximum number of (outer) iterations to take.
    int    mr      = 50;      // the maximum number of (inner) iterations to take.
    
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
    FILE *fp_u, *fp_w, *fp_energy;
    FILE *fp_stif_p, *fp_stif_j, *fp_stif_x;

    /***************************************
          allocate the dynamic arrays
     ****************************************/
    mass_lump       = (double *)malloc(node_num   * sizeof(double));
    mass_coo_i      = (int    *)malloc(nnz        * sizeof(int   ));
    mass_coo_j      = (int    *)malloc(nnz        * sizeof(int   ));
    mass_coo_x      = (double *)malloc(nnz        * sizeof(double));
    stif_coo_i      = (int    *)malloc(nnz        * sizeof(int   ));
    stif_coo_j      = (int    *)malloc(nnz        * sizeof(int   ));
    stif_coo_x      = (double *)malloc(nnz        * sizeof(double));
    mass_csr_p      = (int    *)malloc(csr_p_size * sizeof(int   ));
    mass_csr_j_temp = (int    *)malloc(nnz        * sizeof(int   ));
    mass_csr_x_temp = (double *)malloc(nnz        * sizeof(double));
    stif1_csr_p     = (int    *)malloc(csr_p_size * sizeof(int   ));
    stif2_csr_p     = (int    *)malloc(csr_p_size * sizeof(int   ));
    stif3_csr_p     = (int    *)malloc(csr_p_size * sizeof(int   ));
    stif4_csr_p     = (int    *)malloc(csr_p_size * sizeof(int   ));
    stif_csr_j_temp = (int    *)malloc(nnz        * sizeof(int   ));
    stif_csr_x_temp = (double *)malloc(nnz        * sizeof(double));
    rho             = (double *)malloc(node_num   * sizeof(double));
    vp              = (double *)malloc(node_num   * sizeof(double));
    vs              = (double *)malloc(node_num   * sizeof(double));
    c               = (double**)malloc(4          * sizeof(double));
    for(i = 0; i < 4; i++) c[i] = malloc(node_num * sizeof(double));
    mpml_dx         = (double *)malloc(node_num   * sizeof(double));
    mpml_dy         = (double *)malloc(node_num   * sizeof(double));
    mpml_dxx        = (double *)malloc(node_num   * sizeof(double));
    mpml_dyy        = (double *)malloc(node_num   * sizeof(double));
    mpml_dxx_pyx    = (double *)malloc(node_num   * sizeof(double));
    mpml_dyy_pxy    = (double *)malloc(node_num   * sizeof(double));
    // time evolution
    energy   = (double *)malloc(    step * sizeof(double));
    U_now    = (double *)malloc(node_num * sizeof(double));
    W_now    = (double *)malloc(node_num * sizeof(double));
    U1_now   = (double *)malloc(node_num * sizeof(double));
    U2_now   = (double *)malloc(node_num * sizeof(double));
    U3_now   = (double *)malloc(node_num * sizeof(double));
    U1_now   = (double *)malloc(node_num * sizeof(double));
    W2_now   = (double *)malloc(node_num * sizeof(double));
    W3_now   = (double *)malloc(node_num * sizeof(double));
    U1t_now  = (double *)malloc(node_num * sizeof(double));
    U2t_now  = (double *)malloc(node_num * sizeof(double));
    U3t_now  = (double *)malloc(node_num * sizeof(double));
    W1t_now  = (double *)malloc(node_num * sizeof(double));
    W2t_now  = (double *)malloc(node_num * sizeof(double));
    W3t_now  = (double *)malloc(node_num * sizeof(double));
    U1tt_now = (double *)malloc(node_num * sizeof(double));
    U2tt_now = (double *)malloc(node_num * sizeof(double));
    U3tt_now = (double *)malloc(node_num * sizeof(double));
    W1tt_now = (double *)malloc(node_num * sizeof(double));
    W2tt_now = (double *)malloc(node_num * sizeof(double));
    W3tt_now = (double *)malloc(node_num * sizeof(double));
    Lx1_now  = (double *)malloc(node_num * sizeof(double));
    Lx2_now  = (double *)malloc(node_num * sizeof(double));
    Lx3_now  = (double *)malloc(node_num * sizeof(double));
    Lx4_now  = (double *)malloc(node_num * sizeof(double));
    Ly1_now  = (double *)malloc(node_num * sizeof(double));
    Ly2_now  = (double *)malloc(node_num * sizeof(double));
    Ly3_now  = (double *)malloc(node_num * sizeof(double));
    Ly4_now  = (double *)malloc(node_num * sizeof(double));
    y1       = (double *)malloc(node_num * sizeof(double));
    y2       = (double *)malloc(node_num * sizeof(double));
    y3       = (double *)malloc(node_num * sizeof(double));
    y4       = (double *)malloc(node_num * sizeof(double));
    y5       = (double *)malloc(node_num * sizeof(double));
    rhs_u1   = (double *)malloc(node_num * sizeof(double));
    rhs_u2   = (double *)malloc(node_num * sizeof(double));
    rhs_u3   = (double *)malloc(node_num * sizeof(double));
    rhs_u4   = (double *)malloc(node_num * sizeof(double));

    
    elastic_model(node_num, element_num, element_order, element_node, node_xy, rho, vp, vs, c);
    vp_max = vp[0];
    for( i = 0; i < node_num; i++ )
    {
        if(vp[i] > vp_max) vp_max = vp[i];
    }

    absorbing_boundary_mpml(node_num, element_num, element_order, element_node, node_xy, pml_nx, pml_ny, edge_size, xmin, xmax, ymin, ymax, vp_max, \
                          use_mpml_xmin, use_mpml_xmax, use_mpml_ymin, use_mpml_ymax, mpml_dx, mpml_dy, mpml_dxx, mpml_dyy, mpml_dxx_pyx, mpml_dyy_pxy);
 
    mass_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, mass_coo_i, mass_coo_j, mass_coo_x, 1);
    for (i = 0; i < node_num; i++)  mass_lump[i] = mass_coo_x[i];
    
    /********************************************************
        assemble matrices, first in coo then call 
        fortran subroutines to convert coo to csr. 
        do not need use & to get the address of the pointers.
     *********************************************************/
    // mass
    mass_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, mass_coo_i, mass_coo_j, mass_coo_x, 0);
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size,  &mass_csr_size,  mass_coo_i,  mass_coo_j,  mass_coo_x,  mass_csr_p,  mass_csr_j_temp,  mass_csr_x_temp);
    mass_csr_j  = (   int *)malloc(mass_csr_size * sizeof(int)   );
    mass_csr_x  = (double *)malloc(mass_csr_size * sizeof(double));
    for (i = 0; i < mass_csr_size; i++)
    {
        mass_csr_j[i] = mass_csr_j_temp[i];
        mass_csr_x[i] = mass_csr_x_temp[i];
    }
    // stiffness1: dphidx * dphidx
    stif_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, stif_coo_i, stif_coo_j, stif_coo_x, 1); 
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &stif1_csr_size, stif_coo_i, stif_coo_j, stif_coo_x, stif1_csr_p, stif_csr_j_temp, stif_csr_x_temp);
    stif1_csr_j = (   int *)malloc(stif1_csr_size * sizeof(int)   );
    stif1_csr_x = (double *)malloc(stif1_csr_size * sizeof(double));
    for (i = 0; i < stif1_csr_size; i++)
    {
        stif1_csr_j[i] = stif_csr_j_temp[i];
        stif1_csr_x[i] = stif_csr_x_temp[i];
    }
    // stiffness2: dphidy * dphidy
    stif_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, stif_coo_i, stif_coo_j, stif_coo_x, 2);
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &stif2_csr_size, stif_coo_i, stif_coo_j, stif_coo_x, stif2_csr_p, stif_csr_j_temp, stif_csr_x_temp);
    stif2_csr_j = (   int *)malloc(stif2_csr_size * sizeof(int)   );
    stif2_csr_x = (double *)malloc(stif2_csr_size * sizeof(double));
    for (i = 0; i < stif2_csr_size; i++)
    {
        stif2_csr_j[i] = stif_csr_j_temp[i];
        stif2_csr_x[i] = stif_csr_x_temp[i];
    }
    // stiffness3: dphidx * dphidy
    stif_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, stif_coo_i, stif_coo_j, stif_coo_x, 3);
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &stif3_csr_size, stif_coo_i, stif_coo_j, stif_coo_x, stif3_csr_p, stif_csr_j_temp, stif_csr_x_temp);
    stif3_csr_j = (   int *)malloc(stif3_csr_size * sizeof(int)   );
    stif3_csr_x = (double *)malloc(stif3_csr_size * sizeof(double));
    for (i = 0; i < stif3_csr_size; i++)
    {
        stif3_csr_j[i] = stif_csr_j_temp[i];
        stif3_csr_x[i] = stif_csr_x_temp[i];
    }
    // stiffness4: dphidy * dphidx
    stif_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, stif_coo_i, stif_coo_j, stif_coo_x, 4);
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &stif4_csr_size, stif_coo_i, stif_coo_j, stif_coo_x, stif4_csr_p, stif_csr_j_temp, stif_csr_x_temp);
    stif4_csr_j = (   int *)malloc(stif4_csr_size * sizeof(int)   );
    stif4_csr_x = (double *)malloc(stif4_csr_size * sizeof(double));
    for (i = 0; i < stif4_csr_size; i++)
    {
        stif4_csr_j[i] = stif_csr_j_temp[i];
        stif4_csr_x[i] = stif_csr_x_temp[i];
    }
    // stiffness5: phi * dphidx
    stif_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, stif_coo_i, stif_coo_j, stif_coo_x, 5);
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &stif5_csr_size, stif_coo_i, stif_coo_j, stif_coo_x, stif5_csr_p, stif_csr_j_temp, stif_csr_x_temp);
    stif5_csr_j = (   int *)malloc(stif5_csr_size * sizeof(int));
    stif5_csr_x = (double *)malloc(stif5_csr_size * sizeof(double));
    for (i = 0; i < stif5_csr_size; i++)
    {
        stif5_csr_j[i] = stif_csr_j_temp[i];
        stif5_csr_x[i] = stif_csr_x_temp[i];
    }
    // stiffness6: phi * dphidy
    stif_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, stif_coo_i, stif_coo_j, stif_coo_x, 6);
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &stif6_csr_size, stif_coo_i, stif_coo_j, stif_coo_x, stif6_csr_p, stif_csr_j_temp, stif_csr_x_temp);
    stif6_csr_j = (   int *)malloc(stif6_csr_size * sizeof(int));
    stif6_csr_x = (double *)malloc(stif6_csr_size * sizeof(double));
    for (i = 0; i < stif6_csr_size; i++)
    {
        stif6_csr_j[i] = stif_csr_j_temp[i];
        stif6_csr_x[i] = stif_csr_x_temp[i];
    }
    

    fp_u      = fopen("u.dat", "w");
    fp_w      = fopen("w.dat", "w");
	fp_energy = fopen("energy.dat", "w");
	fp_stif_p = fopen("stif_p.dat", "w");
	fp_stif_j = fopen("stif_j.dat", "w");
    fp_stif_x = fopen("stif_x.dat", "w");


    for(i=0;i<csr_p_size;i++)
    {
        fprintf(fp_stif_p,"%d\n",stif1_csr_p[i]);
    }
    for(i=0;i<stif4_csr_size;i++)
    {
        fprintf(fp_stif_j,"%d\n",stif1_csr_j[i]);
        fprintf(fp_stif_x,"%f\n",stif1_csr_x[i]);
    }

	// write first two step values: u_old[node_num], u_now[node_num], energy[0], energy[1]
	 for (i = 0; i < node_num; i++)
	{
		U1_now[i]   = 0.0;   U2_now[i]   = 0.0;   U3_now[i]   = 0.0;
        U1t_now[i]  = 0.0;   U2t_now[i]  = 0.0;   U3t_now[i]  = 0.0;
        U1tt_now[i] = 0.0;   U2tt_now[i] = 0.0;   U3tt_now[i] = 0.0;
        W1_now[i]   = 0.0;   W2_now[i]   = 0.0;   W3_now[i]   = 0.0;
        W1t_now[i]  = 0.0;   W2t_now[i]  = 0.0;   W3t_now[i]  = 0.0;
        W1tt_now[i] = 0.0;   W2tt_now[i] = 0.0;   W3tt_now[i] = 0.0;
        Lx1_now[i]  = 0.0;   Lx2_now[i]  = 0.0;   Lx3_now[i]  = 0.0;   Lx4_now[i]  = 0.0;
        Ly1_now[i]  = 0.0;   Ly1_now[i]  = 0.0;   Ly1_now[i]  = 0.0;   Ly4_now[i]  = 0.0;
        U_now[i]    = 0.0;   W_now[i]    = 0.0;
		fprintf(fp_u, "%f	", U_now[i]);
		fprintf(fp_w, "%f	", U_now[i]);
	}
	fprintf(fp_u, "\n");
	fprintf(fp_w, "\n");

	energy[0] = 0.0;
	energy[1] = 0.0;
	fprintf(fp_energy, "%f\n", energy[0]);
	fprintf(fp_energy, "%f\n", energy[1]);

	// begin iteration: from 0 to step-1, time = (step + 1) * dt
    printf("\nTime iteration begin:\n");
	for (it = 2; it < step; it++)
	{
        time = (it + 1) * dt;
        if ((it + 1) % 100 == 0) printf("\n Iteration step: %-d, time: %-f\n ", it + 1, time);
		energy[it] = 0.0;
        point_source = seismic_source(f0, t0, 1.0, time);
        for (i = 0; i < node_num; i++)
        {
            rhs1[i] = 0.0; rhs2[i] = 0.0; rhs3[i] = 0.0; rhs4[i] = 0.0;
        }
        /* 
         c[0] c11
         c[1] c13
         c[2] c33
         c[3] c44
        */

       /********************************************************************************************************************************************
        Equation u1:
         mass * U1tt_new = - c11 * dphidx * dphidx * U_now - 2.0 * rho * mpml_dx * phi * phi * U1t_now - rho * mpml_dx * mpml_dx * phi * phi * U1_now 
                           + phi * phi * Lx1_now + phi * phi * Lx2_now + Source_x
        *********************************************************************************************************************************************/
		for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; y3[i] = 0.0; y4[i] = 0.0; y5[i] = 0.0;
        }
        csr_matvec(csr_p_size, stif1_csr_p, stif1_csr_j, stif1_csr_x, U_now  , y1); // dphidx * dphidx * U_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, U1t_now, y2); // phi    * phi    * U1t_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, U1_now , y3); // phi    * phi    * U1_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Lx1_now, y4); // phi    * phi    * Lx1_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Lx2_now, y5); // phi    * phi    * Lx2_now
        for (i = 0; i < node_num; i++)
        {
            rhs_u1[i]  =  - c[0][i] * y1[i] - 2.0 * rho[i] * mpml_dx[i] * y2[i] -  rho[i] * mpml_dx[i] * mpml_dx[i] * y3[i] 
                        + y4[i] + y5[i] + (i == source_node) * point_source * sin( Angle_force * pi / 180.0 );
        }
       
       /********************************************************************************************************************************************
        Equation u2:
         mass * U2tt_new = - c13 * dphidx * dphidy * W_now - c44 * dphidy * dphidx * W_now - rho * mpml_dx * phi * phi * U2t_now 
                           - rho * mpml_dy * phi * phi * U2t_now - rho * mpml_dx * mpml_dy * phi * phi * U2_now
        *********************************************************************************************************************************************/
        for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; y3[i] = 0.0; y4[i] = 0.0;
        }
        csr_matvec(csr_p_size, stif3_csr_p, stif3_csr_j, stif3_csr_x, W_now  , y1); // dphidx * dphidy * W_now
        csr_matvec(csr_p_size, stif4_csr_p, stif4_csr_j, stif4_csr_x, W_now  , y2); // dphidy * dphidx * W_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, U2t_now, y3); // phi    * phi    * U2t_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, U2_now , y4); // phi    * phi    * U2_now
        for (i = 0; i < node_num; i++)
        {
            rhs_u2[i]  =  - c[1][i] * y1[i] - c[3][i] * y2[i] - rho[i] * mpml_dx[i] * y3[i] - rho[i] * mpml_dy[i] * y3[i] - rho[i] * mpml_dx[i] * mpml_dy[i] * y4[i]; 
        }

       /********************************************************************************************************************************************
        Equation u3:
         mass * U3tt_new = - c44 * dphidy * dphidy * U_now - 2.0 * rho * mpml_dy * phi * phi * U3t_now - rho * mpml_dy * mpml_dy * phi * phi * U3_now 
                           + phi * phi * Lx3_now + phi * phi * Lx4_now
        *********************************************************************************************************************************************/
        for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; y3[i] = 0.0; y4[i] = 0.0; y5[i] = 0.0;
        }
        csr_matvec(csr_p_size, stif2_csr_p, stif2_csr_j, stif2_csr_x, U_now  , y1); // dphidy * dphidy * U_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, U3t_now, y2); // phi    * phi    * U3t_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, U3_now , y3); // phi    * phi    * U3_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Lx3_now, y4); // phi    * phi    * Lx3_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Lx4_now, y5); // phi    * phi    * Lx4_now
        for (i = 0; i < node_num; i++)
        {
            rhs_u3[i]  =  - c[3][i] * y1[i] - 2.0 * rho[i] * mpml_dy[i] * y2[i] -  rho[i] * mpml_dy[i] * mpml_dy[i] * y3[i] + y4[i] + y5[i];
        }

        /********************************************************************************************************************************************
         Equation u4:
          mass * Lx1_new = - dt * c11 * mpml_dxx * phi * dphidx * U_now - dt * rho * mpml_dx * phi * phi * Lx1_now - phi * phi * Lx1_now
        *********************************************************************************************************************************************/
        for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; 
        }
        csr_matvec(csr_p_size, stif5_csr_p, stif5_csr_j, stif5_csr_x, U_now  , y1); // phi    * dphidx * U_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Lx1_now, y2); // phi    * phi    * Lx1_now
        for (i = 0; i < node_num; i++)
        {
            rhs_u4[i]  =  - dt * c[0][i] * mpml_dxx[i] * y1[i] - dt * rho[i] * mpml_dx[i] * y2[i] - y2[i];
        }

       /********************************************************************************************************************************************
         Equation u5:
          mass * Lx2_new = - dt * c44 * mpml_dyy_pxy * phi * dphidx * W_now - dt * rho * mpml_dy * phi * phi * Lx2_now - phi * phi * Lx2_now
        *********************************************************************************************************************************************/
        for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; 
        }
        csr_matvec(csr_p_size, stif5_csr_p, stif5_csr_j, stif5_csr_x, W_now  , y1); // phi    * dphidx * W_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Lx2_now, y2); // phi    * phi    * Lx2_now
        for (i = 0; i < node_num; i++)
        {
            rhs_u5[i]  =  - dt * c[3][i] * mpml_dyy_pxy[i] * y1[i] - dt * rho[i] * mpml_dx[i] * y2[i] - y2[i];
        }
        
        /********************************************************************************************************************************************
         Equation u6:
          mass * Lx3_new = - dt * c13 * mpml_dxx_pyx * phi * dphidy * W_now - dt * rho * mpml_dx * phi * phi * Lx3_now - phi * phi * Lx3_now
        *********************************************************************************************************************************************/
        for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; 
        }
        csr_matvec(csr_p_size, stif6_csr_p, stif6_csr_j, stif6_csr_x, W_now  , y1); // phi    * dphidy * W_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Lx3_now, y2); // phi    * phi    * Lx3_now
        for (i = 0; i < node_num; i++)
        {
            rhs_u6[i]  =  - dt * c[1][i] * mpml_dxx_pyx[i] * y1[i] - dt * rho[i] * mpml_dx[i] * y2[i] - y2[i];
        }

       /********************************************************************************************************************************************
         Equation u7:
         mass * Lx4_new = - dt * c44 * mpml_dyy * phi * dphidy * U_now - dt * rho * mpml_dy * phi * phi * Lx4_now - phi * phi * Lx4_now
        *********************************************************************************************************************************************/
        for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; 
        }
        csr_matvec(csr_p_size, stif6_csr_p, stif6_csr_j, stif6_csr_x, U_now  , y1); // phi    * dphidy * U_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Lx4_now, y2); // phi    * phi    * Lx4_now
        for (i = 0; i < node_num; i++)
        {
            rhs_u7[i]  =  - dt * c[3][i] * mpml_dyy[i] * y1[i] - dt * rho[i] * mpml_dx[i] * y2[i] - y2[i];
        }
        
       /********************************************************************************************************************************************
        Equation w1:
         mass * W1tt_new = - c44 * dphidx * dphidx * W_now - 2.0 * rho * mpml_dx * phi * phi * W1t_now - rho * mpml_dx * mpml_dx * phi * phi * W1_now 
                           + phi * phi * Ly1_now + phi * phi * Ly2_now + Source_y
        *********************************************************************************************************************************************/
		for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; y3[i] = 0.0; y4[i] = 0.0; y5[i] = 0.0;
        }
        csr_matvec(csr_p_size, stif1_csr_p, stif1_csr_j, stif1_csr_x, W_now  , y1); // dphidx * dphidx * W_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, W1t_now, y2); // phi    * phi    * W1t_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, W1_now , y3); // phi    * phi    * W1_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Ly1_now, y4); // phi    * phi    * Ly1_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Ly2_now, y5); // phi    * phi    * Ly2_now
        for (i = 0; i < node_num; i++)
        {
            rhs_w1[i]  =  - c[3][i] * y1[i] - 2.0 * rho[i] * mpml_dx[i] * y2[i] - rho[i] * mpml_dx[i] * mpml_dx[i] * y3[i] 
                        + y4[i] + y5[i] + (i == source_node) * point_source * cos( Angle_force * pi / 180.0 );
        }
       
       /********************************************************************************************************************************************
        Equation w2:
         mass * W2tt_new = - c44 * dphidx * dphidy * U_now - c13 * dphidy * dphidx * U_now - rho * mpml_dx * phi * phi * W2t_now 
                           - rho * mpml_dy * phi * phi * W2t_now - rho * mpml_dx * mpml_dy * phi * phi * W2_now
        *********************************************************************************************************************************************/
        for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; y3[i] = 0.0; y4[i] = 0.0;
        }
        csr_matvec(csr_p_size, stif3_csr_p, stif3_csr_j, stif3_csr_x, U_now  , y1); // dphidx * dphidy * U_now
        csr_matvec(csr_p_size, stif4_csr_p, stif4_csr_j, stif4_csr_x, U_now  , y2); // dphidy * dphidx * U_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, W2t_now, y3); // phi    * phi    * W2t_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, W2_now , y4); // phi    * phi    * W2_now
        for (i = 0; i < node_num; i++)
        {
            rhs_w2[i]  =  - c[3][i] * y1[i] - c[1][i] *  y2[i] - rho[i] * mpml_dx[i] * y3[i] - rho[i] * mpml_dy[i] * y3[i] - rho[i] * mpml_dx[i] * mpml_dy[i] * y4[i]; 
        }

       /********************************************************************************************************************************************
        Equation w3:
         mass * W3tt_new = - c33 * dphidy * dphidy * W_now - 2.0 * rho * mpml_dy * phi * phi * W3t_now - rho * mpml_dy * mpml_dy * phi * phi * W3_now 
                           + phi * phi * Ly3_now + phi * phi * Ly4_now
        *********************************************************************************************************************************************/
        for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; y3[i] = 0.0; y4[i] = 0.0; y5[i] = 0.0;
        }
        csr_matvec(csr_p_size, stif2_csr_p, stif2_csr_j, stif2_csr_x, W_now  , y1); // dphidy * dphidy * W_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, W3t_now, y2); // phi    * phi    * W3t_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, W3_now , y3); // phi    * phi    * W3_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Ly3_now, y4); // phi    * phi    * Ly3_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Ly4_now, y5); // phi    * phi    * Ly4_now
        for (i = 0; i < node_num; i++)
        {
            rhs_w3[i]  =  - c[2][i] * y1[i] - 2.0 * rho[i] * mpml_dy[i] * y2[i] -  rho[i] * mpml_dy[i] * mpml_dy[i] * y3[i] + y4[i] + y5[i];
        }

        /********************************************************************************************************************************************
         Equation 4:
          mass * Ly1_new = - dt * c44 * mpml_dxx * phi * dphidx * W_now - dt * rho * mpml_dx * phi * phi * Ly1_now - phi * phi * Ly1_now
        *********************************************************************************************************************************************/
        for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; 
        }
        csr_matvec(csr_p_size, stif5_csr_p, stif5_csr_j, stif5_csr_x, W_now  , y1); // phi    * dphidx * W_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Ly1_now, y2); // phi    * phi    * Ly1_now
        for (i = 0; i < node_num; i++)
        {
            rhs_w4[i]  =  - dt * c[3][i] * mpml_dxx[i] * y1[i] - dt * rho[i] * mpml_dx[i] * y2[i] - y2[i];
        }

       /********************************************************************************************************************************************
         Equation w5:
          mass * Ly2_new = - dt * c13 * mpml_dyy_pxy * phi * dphidx * U_now - dt * rho * mpml_dy * phi * phi * Ly2_now - phi * phi * Ly2_now
        *********************************************************************************************************************************************/
        for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; 
        }
        csr_matvec(csr_p_size, stif5_csr_p, stif5_csr_j, stif5_csr_x, U_now  , y1); // phi    * dphidx * U_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Ly2_now, y2); // phi    * phi    * Ly2_now
        for (i = 0; i < node_num; i++)
        {
            rhs_w5[i]  =  - dt * c[1][i] * mpml_dyy_pxy[i] * y1[i] - dt * rho[i] * mpml_dx[i] * y2[i] - y2[i];
        }
        
        /********************************************************************************************************************************************
         Equation w6:
          mass * Ly3_new = - dt * c44 * mpml_dxx_pyx * phi * dphidy * U_now - dt * rho * mpml_dx * phi * phi * Ly3_now - phi * phi * Ly3_now
        *********************************************************************************************************************************************/
        for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; 
        }
        csr_matvec(csr_p_size, stif6_csr_p, stif6_csr_j, stif6_csr_x, U_now  , y1); // phi    * dphidy * U_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Ly3_now, y2); // phi    * phi    * Ly3_now
        for (i = 0; i < node_num; i++)
        {
            rhs_w6[i]  =  - dt * c[3][i] * mpml_dxx_pyx[i] * y1[i] - dt * rho[i] * mpml_dx[i] * y2[i] - y2[i];
        }

       /********************************************************************************************************************************************
         Equation w7:
         mass * Ly4_new = - dt * c33 * mpml_dyy * phi * dphidy * W_now - dt * rho * mpml_dy * phi * phi * Ly4_now - phi * phi * Ly4_now
        *********************************************************************************************************************************************/
        for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; 
        }
        csr_matvec(csr_p_size, stif6_csr_p, stif6_csr_j, stif6_csr_x, W_now  , y1); // phi    * dphidy * W_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, Ly4_now, y2); // phi    * phi    * Ly4_now
        for (i = 0; i < node_num; i++)
        {
            rhs_w7[i]  =  - dt * c[2][i] * mpml_dyy[i] * y1[i] - dt * rho[i] * mpml_dx[i] * y2[i] - y2[i];
        }
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        // here can use openmp
        if( strcmp(solver, "masslump") == 0 )
        {
            for (i = 0; i < node_num; i++)
            {
                if (mass_lump[i] != 0.0)
                {
                    u1_now[i] = rhs1[i] / (1.0 * mass_lump[i]);
                    p1_now[i] = rhs2[i] / (1.0 * mass_lump[i]);
                    u2_now[i] = rhs3[i] / (1.0 * mass_lump[i]);
                    p2_now[i] = rhs4[i] / (1.0 * mass_lump[i]);
                }
                else  printf("Fatal Error: zero value in the mass_csr_lumped"); // in case zero value
            }
        }
        else if (strcmp(solver, "mgmres") == 0)
        {
            pmgmres_ilu_cr ( node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, u1_now, rhs1, itr_max, mr, tol_abs, tol_rel );
            pmgmres_ilu_cr ( node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, p1_now, rhs2, itr_max, mr, tol_abs, tol_rel );
            pmgmres_ilu_cr ( node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, u2_now, rhs3, itr_max, mr, tol_abs, tol_rel );
            pmgmres_ilu_cr ( node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, p2_now, rhs4, itr_max, mr, tol_abs, tol_rel );
        }
        else
        {
            fprintf(stderr, "\n");
            fprintf(stderr, "ACOUSTIC_PML - Fatal error!\n");
            fprintf(stderr, "  Solver type is not set = \"%s\".\n", solver);
            exit(1);
        }
       
        for (i = 0; i < node_num; i++)
        {
            u_now[i] = u1_now[i] + u2_now[i];
            energy[it] += u_now[i]*u_now[i];
            if( energy[it] > 10e6)
            {
                fprintf(stderr, "\n");
                fprintf(stderr, "ACOUSTIC_PML - Fatal error!\n");
                fprintf(stderr, "Energy exceeds maximum value!\n");
                exit(1);
            }
        }
		if ((it + 1) % 20 == 0)
		{
			for (i = 0; i < node_num; i++)
			{
				fprintf(fp_u, "%f	", u_now[i]);
			}
			fprintf(fp_u, "\n");
		}
		fprintf(fp_energy, "%f\n", energy[it]);
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
    free(mass_csr_j_temp);
    free(mass_csr_x_temp);
    free(stif_csr_j_temp);
    free(stif_csr_x_temp);
    free(rho);
    free(vp);
    free(mpml_dx);
    free(mpml_dy);
    free(mpml_dxx);
    free(mpml_dyy);
    free(mpml_dxx_pyx);
    free(mpml_dyy_pxy);
    free(energy);
    free(u_now);
	free(u1_old);
	free(u2_old);
	free(u1_now);
	free(u2_now);
	free(p1_now);
	free(p2_now);
	free(y1);
	free(y2);
	free(y3);
    free(y4);
	free(rhs1);
    free(rhs2);
	free(rhs3);
	free(rhs4);

	fclose(fp_u);
	fclose(fp_energy);
	fclose(fp_stif_p);
    fclose(fp_stif_j);
	fclose(fp_stif_x);


	printf("\n Acoutsic Normal End!\n");
}
