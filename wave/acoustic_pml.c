#include "../solver/mgmres/mgmres.c"

void acoustic_pml(char *type, int node_num, int element_num, int element_order, int *element_node, double **node_xy, int nnz, int csr_p_size, int step, \
     double dt, double f0, double t0, double edge_size, double xmin, double xmax, double ymin, double ymax, int pml_nx, int pml_ny, int source_node, char *solver)
{
    
    /***************************************
     coo(i,j,x) and csr(p,j,x) arrays
     ****************************************/
    int  mass_csr_size = 0; // csr size for mass : j and x
    int stif1_csr_size = 0; // csr size for stif1: j and x
    int stif2_csr_size = 0; // csr size for stif2: j and x
    int stif3_csr_size = 0; // csr size for stif3: j and x
    int stif4_csr_size = 0; // csr size for stif4: j and x
    int *mass_coo_i = NULL;
    int *mass_coo_j = NULL;
    int *mass_csr_p = NULL;
    int *mass_csr_j = NULL;
    int *mass_csr_j_temp = NULL;
    int *stif1_coo_i = NULL;
    int *stif1_coo_j = NULL;
    int *stif2_coo_i = NULL;
    int *stif2_coo_j = NULL;
    int *stif3_coo_i = NULL;
    int *stif3_coo_j = NULL;
    int *stif4_coo_i = NULL;
    int *stif4_coo_j = NULL;
    int *stif1_csr_p = NULL;
    int *stif1_csr_j = NULL;
    int *stif2_csr_p = NULL;
    int *stif2_csr_j = NULL;
    int *stif3_csr_p = NULL;
    int *stif3_csr_j = NULL;
    int *stif4_csr_p = NULL;
    int *stif4_csr_j = NULL;
    int *stif1_csr_j_temp = NULL;
    int *stif2_csr_j_temp = NULL;
    int *stif3_csr_j_temp = NULL;
    int *stif4_csr_j_temp = NULL;
    double *stif1_coo_x = NULL;
    double *stif1_csr_x = NULL;
    double *stif2_coo_x = NULL;
    double *stif2_csr_x = NULL;
    double *stif3_coo_x = NULL;
    double *stif3_csr_x = NULL;
    double *stif4_coo_x = NULL;
    double *stif4_csr_x = NULL;
    double *stif1_csr_x_temp = NULL;
    double *stif2_csr_x_temp = NULL;
    double *stif3_csr_x_temp = NULL;
    double *stif4_csr_x_temp = NULL;
    double *mass_coo_x = NULL;
    double *mass_csr_x = NULL;
    double *mass_csr_x_temp = NULL;
    double *mass_lump = NULL;
    
    /***************************************
     model density and velocity parameters
     ****************************************/
    double *rho = NULL;
    double *vp = NULL;
    
    /***************************************
             time evolution parameters
     ****************************************/
	int i, it;
	double time;
    double vp_max;
	double point_source;
	double *rhs1   = NULL;
    double *rhs2   = NULL;
    double *rhs3   = NULL;
    double *rhs4   = NULL;
    double *u_now  = NULL;
	double *u1_old = NULL;
	double *u1_now = NULL;
    double *u2_old = NULL;
    double *u2_now = NULL;
    double *p1_now = NULL;
    double *p2_now = NULL;
	double *y1     = NULL;
	double *y2     = NULL;
	double *y3     = NULL;
    double *y4     = NULL;
	double *energy = NULL;
    double tol_abs = 1.0e-08; // a relative tolerance comparing the current residual to the initial residual.
    double tol_rel = 1.0e-08; // an absolute tolerance applied to the current residual.
    int    itr_max = 100;     // the maximum number of (outer) iterations to take.
    int    mr      = 50;      // the maximum number of (inner) iterations to take.
    
    /***************************************
     abosorbing bc, pml and mpml parameters
     ****************************************/
    int use_pml_xmin = 1;
    int use_pml_xmax = 1;
    int use_pml_ymin = 1;
    int use_pml_ymax = 1; // use_pml_ymax = 0, free surface on ymax.
    double *pml_dx  = NULL;
    double *pml_dy  = NULL;
    double *pml_dxx = NULL;
    double *pml_dyy = NULL;
    
    /***************************************
                  file pointers
     ****************************************/
    FILE *fp_u, *fp_energy;
    FILE *fp_stif_p, *fp_stif_j, *fp_stif_x;

    /***************************************
     allocate the dynamic arrays
     ****************************************/
    // mass lump
    mass_lump = (double *)malloc(node_num * sizeof(double));
    // coo (i,j,x), mass and stif matrices
    mass_coo_i = (int *)malloc(nnz * sizeof(int));
    mass_coo_j = (int *)malloc(nnz * sizeof(int));
    mass_coo_x = (double *)malloc(nnz * sizeof(double));
    stif1_coo_i = (int *)malloc(nnz * sizeof(int));
    stif1_coo_j = (int *)malloc(nnz * sizeof(int));
    stif2_coo_i = (int *)malloc(nnz * sizeof(int));
    stif2_coo_j = (int *)malloc(nnz * sizeof(int));
    stif3_coo_i = (int *)malloc(nnz * sizeof(int));
    stif3_coo_j = (int *)malloc(nnz * sizeof(int));
    stif4_coo_i = (int *)malloc(nnz * sizeof(int));
    stif4_coo_j = (int *)malloc(nnz * sizeof(int));
    stif1_coo_x = (double *)malloc(nnz * sizeof(double));
    stif2_coo_x = (double *)malloc(nnz * sizeof(double));
    stif3_coo_x = (double *)malloc(nnz * sizeof(double));
    stif4_coo_x = (double *)malloc(nnz * sizeof(double));

    // csr (p,j,x), mass and stif matrices
    mass_csr_p      = (int *)malloc(csr_p_size * sizeof(int));
    mass_csr_j_temp = (int *)malloc(nnz * sizeof(int));
    mass_csr_x_temp = (double *)malloc(nnz * sizeof(double));
    stif1_csr_p     = (int *)malloc(csr_p_size * sizeof(int));
    stif2_csr_p     = (int *)malloc(csr_p_size * sizeof(int));
    stif3_csr_p     = (int *)malloc(csr_p_size * sizeof(int));
    stif4_csr_p     = (int *)malloc(csr_p_size * sizeof(int));
    stif1_csr_j_temp = (int *)malloc(nnz * sizeof(int));
    stif2_csr_j_temp = (int *)malloc(nnz * sizeof(int));
    stif3_csr_j_temp = (int *)malloc(nnz * sizeof(int));
    stif4_csr_j_temp = (int *)malloc(nnz * sizeof(int));
    stif1_csr_x_temp = (double *)malloc(nnz * sizeof(double));
    stif2_csr_x_temp = (double *)malloc(nnz * sizeof(double));
    stif3_csr_x_temp = (double *)malloc(nnz * sizeof(double));
    stif4_csr_x_temp = (double *)malloc(nnz * sizeof(double));
    // model density and velocity matrices
    rho = (double *)malloc(node_num * sizeof(double));
    vp  = (double *)malloc(node_num * sizeof(double));
    // pml and mpml matrices
    pml_dx  = (double *)malloc(node_num * sizeof(double));
    pml_dy  = (double *)malloc(node_num * sizeof(double));
    pml_dxx = (double *)malloc(node_num * sizeof(double));
    pml_dyy = (double *)malloc(node_num * sizeof(double));
    // time evolution
    energy = (double *)malloc(step * sizeof(double));
    u_now  = (double *)malloc(node_num * sizeof(double));
    u1_old = (double *)malloc(node_num * sizeof(double));
    u1_now = (double *)malloc(node_num * sizeof(double));
    u2_old = (double *)malloc(node_num * sizeof(double));
    u2_now = (double *)malloc(node_num * sizeof(double));
    p1_now = (double *)malloc(node_num * sizeof(double));
    p2_now = (double *)malloc(node_num * sizeof(double));
    y1     = (double *)malloc(node_num * sizeof(double));
    y2     = (double *)malloc(node_num * sizeof(double));
    y3     = (double *)malloc(node_num * sizeof(double));
    y4     = (double *)malloc(node_num * sizeof(double));
    rhs1   = (double *)malloc(node_num * sizeof(double));
    rhs2   = (double *)malloc(node_num * sizeof(double));
    rhs3   = (double *)malloc(node_num * sizeof(double));
    rhs4   = (double *)malloc(node_num * sizeof(double));

    
    acoustic_model(node_num, element_num, element_order, element_node, node_xy, rho, vp);
    vp_max = vp[0];
    for( i = 0; i < node_num; i++ )
    {
        if(vp[i] > vp_max) vp_max = vp[i];
    }
    absorbing_boundary_pml(node_num, element_num, element_order, element_node, node_xy, pml_nx, pml_ny, edge_size, xmin, xmax, ymin, ymax, vp_max,
                           use_pml_xmin, use_pml_xmax, use_pml_ymin, use_pml_ymax, pml_dx, pml_dy, pml_dxx, pml_dyy);
    // get lump mass matrix
    mass_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, mass_coo_i, mass_coo_j, mass_coo_x, 1);
    for (i = 0; i < node_num; i++)  mass_lump[i] = mass_coo_x[i];
    // get mass and stif matrices in coo format
    mass_sparse_all(type, node_num, element_num, element_order, element_node, node_xy,  mass_coo_i,  mass_coo_j,  mass_coo_x, 0);
    stif_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, stif1_coo_i, stif1_coo_j, stif1_coo_x, 2); // dphidx * dphidx
    stif_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, stif2_coo_i, stif2_coo_j, stif2_coo_x, 4); //    phi * dphidx
    stif_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, stif3_coo_i, stif3_coo_j, stif3_coo_x, 3); // dphidy * dphidy
    stif_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, stif4_coo_i, stif4_coo_j, stif4_coo_x, 5); //    phi * dphidy

    // call fortran subroutines to convert coo to csr. Attention: mi, mj, mass, Mp_temp, Mj_temp, Mass_temp are addresses,
    // do not need use & to get the address of the pointers.
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size,  &mass_csr_size,  mass_coo_i,  mass_coo_j,  mass_coo_x,  mass_csr_p,  mass_csr_j_temp,  mass_csr_x_temp);
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &stif1_csr_size, stif1_coo_i, stif1_coo_j, stif1_coo_x, stif1_csr_p, stif1_csr_j_temp, stif1_csr_x_temp);
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &stif2_csr_size, stif2_coo_i, stif2_coo_j, stif2_coo_x, stif2_csr_p, stif2_csr_j_temp, stif2_csr_x_temp);
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &stif3_csr_size, stif3_coo_i, stif3_coo_j, stif3_coo_x, stif3_csr_p, stif3_csr_j_temp, stif3_csr_x_temp);
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &stif4_csr_size, stif4_coo_i, stif4_coo_j, stif4_coo_x, stif4_csr_p, stif4_csr_j_temp, stif4_csr_x_temp);

    // allocate the real size mass and stif csr matrices
    mass_csr_j  = (int *)malloc(mass_csr_size * sizeof(int));
    mass_csr_x  = (double *)malloc(mass_csr_size * sizeof(double));
    stif1_csr_j = (int *)malloc(stif1_csr_size * sizeof(int));
    stif2_csr_j = (int *)malloc(stif2_csr_size * sizeof(int));
    stif3_csr_j = (int *)malloc(stif3_csr_size * sizeof(int));
    stif4_csr_j = (int *)malloc(stif4_csr_size * sizeof(int));
    stif1_csr_x = (double *)malloc(stif1_csr_size * sizeof(double));
    stif2_csr_x = (double *)malloc(stif2_csr_size * sizeof(double));
    stif3_csr_x = (double *)malloc(stif3_csr_size * sizeof(double));
    stif4_csr_x = (double *)malloc(stif4_csr_size * sizeof(double));

    for (i = 0; i < mass_csr_size; i++)
    {
        mass_csr_j[i] = mass_csr_j_temp[i];
        mass_csr_x[i] = mass_csr_x_temp[i];
    }
    for (i = 0; i < stif1_csr_size; i++)
    {
        stif1_csr_j[i] = stif1_csr_j_temp[i];
        stif1_csr_x[i] = stif1_csr_x_temp[i];
    }
    for (i = 0; i < stif2_csr_size; i++)
    {
        stif2_csr_j[i] = stif2_csr_j_temp[i];
        stif2_csr_x[i] = stif2_csr_x_temp[i];
    }
    for (i = 0; i < stif3_csr_size; i++)
    {
        stif3_csr_j[i] = stif3_csr_j_temp[i];
        stif3_csr_x[i] = stif3_csr_x_temp[i];
    }
    for (i = 0; i < stif4_csr_size; i++)
    {
        stif4_csr_j[i] = stif4_csr_j_temp[i];
        stif4_csr_x[i] = stif4_csr_x_temp[i];
    }

    fp_u      = fopen("u.dat", "w");
	fp_energy = fopen("energy.dat", "w");
	fp_stif_p = fopen("stif_p.dat", "w");
	fp_stif_j = fopen("stif_j.dat", "w");
    fp_stif_x = fopen("stif_x.dat", "w");


    for(i=0;i<csr_p_size;i++)
    {
        fprintf(fp_stif_p,"%d\n",stif4_csr_p[i]);
    }
    for(i=0;i<stif4_csr_size;i++)
    {
        fprintf(fp_stif_j,"%d\n",stif4_csr_j[i]);
        fprintf(fp_stif_x,"%f\n",stif4_csr_x[i]);
    }

	// write first two step values: u_old[node_num], u_now[node_num], energy[0], energy[1]
	 for (i = 0; i < node_num; i++)
	{
		u1_old[i] = 0.0;
        u2_old[i] = 0.0;
		u1_now[i] = 0.0;
        u2_now[i] = 0.0;
        p1_now[i] = 0.0;
        p2_now[i] = 0.0;
        u_now[i]  = 0.0;
		fprintf(fp_u, "%f	", u_now[i]);
	}
	fprintf(fp_u, "\n");

	energy[0] = 0.0;
	energy[1] = 0.0;
	fprintf(fp_energy, "%f\n", energy[0]);
	fprintf(fp_energy, "%f\n", energy[1]);

	// begin iteration: from 0 to step-1, time = (step + 1) * dt,  with total iteration number = step
	// give u(0) ans u(1)
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
        
        // equation 1:
        // mass * u1_new = 2*mass*u1_now - mass*u1_old - 2*dt*pml_dx*mass*u1_now + 2*dt*pml_dx*mass*u1_old
        //               - dt*dt*pml_dx*pml_dx*mass*u1_now + vp*vp*dt*dt*stif1*u_now + vp*vp*dt*dt*mass*p1_now
        // stif1: dphidx * dphidx
		for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0; y2[i] = 0.0; y3[i] = 0.0; y4[i] = 0.0;
		}
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, u1_now, y1); // mass  * u1_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, u1_old, y2); // mass  * u1_old
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, p1_now, y3); // mass  * p1_now
		csr_matvec(csr_p_size, stif1_csr_p, stif1_csr_j, stif1_csr_x,  u_now, y4); // stif1 *  u_now
        for (i = 0; i < node_num; i++)
        {
            rhs1[i]  =  2.0 * y1[i] - y2[i] - 2.0 * dt * pml_dx[i] * y1[i] + 2 * dt * pml_dx[i] * y2[i]  \
                     -  dt * dt * pml_dx[i] * pml_dx[i] * y1[i] + vp[i] * vp[i] * dt * dt * y4[i] + vp[i] * vp[i] * dt * dt * y3[i]  \
                     +  (i == source_node) * point_source ;
        }
       
        // equation 2
        // mass * p1_new = mass*p1_now - dt*pml_dx*mass*p1_now - dt*pml_dxx*stif2*u_now
        // stif2: phi * dphidx
        for (i = 0; i < node_num; i++)
        {
            y1[i] = 0.0; y2[i] = 0.0; 
        }
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, p1_now, y1); // mass  * p1_now
        csr_matvec(csr_p_size, stif2_csr_p, stif2_csr_j, stif2_csr_x,  u_now, y2); // stif2 *  u_now
        for (i = 0; i < node_num; i++)
        {
            rhs2[i]  =  y1[i] - dt * pml_dx[i] * y1[i] - dt * pml_dxx[i] * y2[i];
        }
        
        // equation 3:
        // mass * u2_new = 2*mass*u2_now - mass*u2_old - 2*dt*pml_dy*mass*u2_now + 2*dt*pml_dy*mass*u2_old
        //               - dt*dt*pml_dy*pml_dy*mass*u2_now + vp*vp*dt*dt*stif3*u_now + vp*vp*dt*dt*mass*p2_now
        // stif3: dphidy * dphidy
        for (i = 0; i < node_num; i++)
        {
            y1[i] = 0.0; y2[i] = 0.0; y3[i] = 0.0; y4[i] = 0.0;
        }
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, u2_now, y1); // mass  * u2_now
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, u2_old, y2); // mass  * u2_old
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, p2_now, y3); // mass  * p2_now
        csr_matvec(csr_p_size, stif3_csr_p, stif3_csr_j, stif3_csr_x,  u_now, y4); // stif3 *  u_now
        for (i = 0; i < node_num; i++)
        {
            rhs3[i]  =  2.0 * y1[i] - y2[i] - 2.0 * dt * pml_dy[i] * y1[i] + 2 * dt * pml_dy[i] * y2[i]  \
                     -  dt * dt * pml_dy[i] * pml_dy[i] * y1[i] + vp[i] * vp[i] * dt * dt * y4[i] + vp[i] * vp[i] * dt * dt * y3[i];
        }
        
        // equation 4
        // mass * p2_new = mass*p2_now - dt*pml_dy*mass*p2_now - dt*pml_dyy*stif4*u_now
        // stif4: phi * dphidy
        for (i = 0; i < node_num; i++)
        {
            y1[i] = 0.0; y2[i] = 0.0; 
        }
        csr_matvec(csr_p_size,  mass_csr_p,  mass_csr_j,  mass_csr_x, p2_now, y1); // mass  * p2_now
        csr_matvec(csr_p_size, stif4_csr_p, stif4_csr_j, stif4_csr_x,  u_now, y2); // stif4 *  u_now
        for (i = 0; i < node_num; i++)
        {
            rhs4[i]  =  y1[i] - dt * pml_dy[i] * y1[i] - dt * pml_dyy[i] * y2[i];
        }
  
        for (i = 0; i < node_num; i++)
        {
            u1_old[i] = u1_now[i];
            u2_old[i] = u2_now[i];
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
    free(stif1_coo_i);
    free(stif1_coo_j);
    free(stif1_coo_x);
    free(stif2_coo_i);
    free(stif2_coo_j);
    free(stif2_coo_x);
    free(stif3_coo_i);
    free(stif3_coo_j);
    free(stif3_coo_x);
    free(stif4_coo_i);
    free(stif4_coo_j);
    free(stif4_coo_x);
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
    free(stif1_csr_j_temp);
    free(stif1_csr_x_temp);
    free(stif2_csr_j_temp);
    free(stif2_csr_x_temp);
    free(stif3_csr_j_temp);
    free(stif3_csr_x_temp);
    free(stif4_csr_j_temp);
    free(stif4_csr_x_temp);
    free(rho);
    free(vp);
    free(pml_dx);
    free(pml_dy);
    free(pml_dxx);
    free(pml_dyy);
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
