
void acoustic(char *type, int node_num, int element_num, int element_order, int *element_node, double **node_xy, int nnz, int csr_p_size, int step, double dt, double f0, double t0, int source_node, char *solver)
{
    
    /***************************************
     coo(i,j,x) and csr(p,j,x) arrays
     ****************************************/
    int mass_csr_size = 0; // csr size for mass: j and x
    int stif_csr_size = 0; // csr size for stif: j and x
    int *mass_coo_i = NULL;
    int *mass_coo_j = NULL;
    int *stif_coo_i = NULL;
    int *stif_coo_j = NULL;
    int *mass_csr_p = NULL;
    int *mass_csr_j = NULL;
    int *stif_csr_p = NULL;
    int *stif_csr_j = NULL;
    int *mass_csr_j_temp = NULL;
    int *stif_csr_j_temp = NULL;
    double *mass_coo_x = NULL;
    double *stif_coo_x = NULL;
    double *mass_csr_x = NULL;
    double *stif_csr_x = NULL;
    double *mass_csr_x_temp = NULL;
    double *stif_csr_x_temp = NULL;
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
	double point_source;
	double *rhs   = NULL;
	double *u_old = NULL;
	double *u_now = NULL;
	double *y1    = NULL;
	double *y2    = NULL;
	double *y3    = NULL;
	double *energy= NULL;
    double tol_abs = 1.0e-08; // a relative tolerance comparing the current residual to the initial residual.
    double tol_rel = 1.0e-08; // an absolute tolerance applied to the current residual.
    int    itr_max = 100;     // the maximum number of (outer) iterations to take.
    int    mr      = 50;      // the maximum number of (inner) iterations to take.
    /***************************************
                  file pointers
     ****************************************/
    FILE *fp_u, *fp_energy;

    /***************************************
     allocate the dynamic arrays
     ****************************************/
    
    // mass lump
    mass_lump = (double *)malloc(node_num * sizeof(double));
    // coo (i,j,x), mass and stif matrices
    mass_coo_i = (int *)malloc(nnz * sizeof(int));
    stif_coo_i = (int *)malloc(nnz * sizeof(int));
    mass_coo_j = (int *)malloc(nnz * sizeof(int));
    stif_coo_j = (int *)malloc(nnz * sizeof(int));
    mass_coo_x = (double *)malloc(nnz * sizeof(double));
    stif_coo_x = (double *)malloc(nnz * sizeof(double));
    // csr (p,j,x), mass and stif matrices
    mass_csr_p = (int *)malloc(csr_p_size * sizeof(int));
    stif_csr_p = (int *)malloc(csr_p_size * sizeof(int));
    mass_csr_j_temp = (int *)malloc(nnz * sizeof(int));
    stif_csr_j_temp = (int *)malloc(nnz * sizeof(int));
    mass_csr_x_temp = (double *)malloc(nnz * sizeof(double));
    stif_csr_x_temp = (double *)malloc(nnz * sizeof(double));
    // model density and velocity matrices
    rho = (double *)malloc(node_num * sizeof(double));
    vp  = (double *)malloc(node_num * sizeof(double));
    rhs = (double *)malloc(node_num * sizeof(double));
    // time evolution
    u_old = (double *)malloc(node_num * sizeof(double));
    u_now = (double *)malloc(node_num * sizeof(double));
    y1    = (double *)malloc(node_num * sizeof(double));
    y2    = (double *)malloc(node_num * sizeof(double));
    y3    = (double *)malloc(node_num * sizeof(double));
    energy= (double *)malloc(node_num * sizeof(double));
    
    acoustic_model(node_num, element_num, element_order, element_node, node_xy, rho, vp);
    
    // get lump mass matrix
    mass_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, mass_coo_i, mass_coo_j, mass_coo_x, 1);
    for (i = 0; i < node_num; i++)  mass_lump[i] = mass_coo_x[i];
    // get mass and stif matrices in coo format
    mass_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, mass_coo_i, mass_coo_j, mass_coo_x, 0);
    stif_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, stif_coo_i, stif_coo_j, stif_coo_x, 1);
    // call fortran subroutines to convert coo to csr. Attention: mi, mj, mass, Mp_temp, Mj_temp, Mass_temp are addresses,
    // do not need use & to get the address of the pointers.
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &mass_csr_size, mass_coo_i, mass_coo_j, mass_coo_x, mass_csr_p, mass_csr_j_temp, mass_csr_x_temp);
    __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &stif_csr_size, stif_coo_i, stif_coo_j, stif_coo_x, stif_csr_p, stif_csr_j_temp, stif_csr_x_temp);
    // allocate the real size mass and stif csr matrices
    mass_csr_j = (int *)malloc(mass_csr_size * sizeof(int));
    mass_csr_x = (double *)malloc(mass_csr_size * sizeof(double));
    stif_csr_j = (int *)malloc(stif_csr_size * sizeof(int));
    stif_csr_x = (double *)malloc(stif_csr_size * sizeof(double));
    for (i = 0; i < mass_csr_size; i++)
    {
        mass_csr_j[i] = mass_csr_j_temp[i];
        mass_csr_x[i] = mass_csr_x_temp[i];
    }
    for (i = 0; i < stif_csr_size; i++)
    {
        stif_csr_j[i] = stif_csr_j_temp[i];
        stif_csr_x[i] = stif_csr_x_temp[i];
    }

	fp_u = fopen("u.dat", "w");
	fp_energy = fopen("energy.dat", "w");

	// write first two step values: u_old[node_num], u_now[node_num], energy[0], energy[1]
	for (i = 0; i < node_num; i++)
	{
		u_old[i] = 0.0;
		u_now[i] = 0.0;
		fprintf(fp_u, "%f	", u_old[i]);
	}
	fprintf(fp_u, "\n");

	energy[0] = 0.0;
	energy[1] = 0.0;
	fprintf(fp_energy, "%f\n", energy[0]);
	fprintf(fp_energy, "%f\n", energy[1]);

	// begin iteration: from 0 to step-1, time = (step + 1) * dt,  with total iteration number = step
	// give u(0) ans u(1)
    
	for (it = 2; it < step; it++)
	{
		time = (it + 1) * dt;
		energy[it] = 0.0;
		if ((it + 1) % 100 == 0) printf("iteration step: %-d, time: %-f\n ", it + 1, time);

		for (i = 0; i < node_num; i++)
		{
			y1[i] = 0.0;
			y2[i] = 0.0;
			y3[i] = 0.0;
			rhs[i] = 0.0;
		}
		point_source = seismic_source(f0, t0, 1.0, time);
		csr_matvec(csr_p_size, stif_csr_p, stif_csr_j, stif_csr_x, u_now, y1); // stif * u_now
		csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, u_now, y2); // mass * u_now
		csr_matvec(csr_p_size, mass_csr_p, mass_csr_j, mass_csr_x, u_old, y3); // mass * u_old
        
        for (i = 0; i < node_num; i++)
        {
            rhs[i]   = (i == source_node) * point_source -  vp[i] * vp[i] * dt * dt * y1[i] + 2.0 * y2[i] - y3[i];
            u_old[i] = u_now[i];
        }
		// here can use openmp
        if( strcmp(solver, "masslump") == 0 )
        {
            for (i = 0; i < node_num; i++)
            {
                if (mass_lump[i] != 0.0)  u_now[i] = rhs[i] / (1.0 * mass_lump[i]);
                else  printf("Fatal Error: zero value in the mass_csr_lumped"); // in case zero value
            }
        }
        else if (strcmp(solver, "mgmres") == 0)
        {
            pmgmres_ilu_cr ( node_num, mass_csr_size, mass_csr_p, mass_csr_j, mass_csr_x, u_now, rhs, itr_max, mr, tol_abs, tol_rel );
        }
        else
        {
            fprintf(stderr, "\n");
            fprintf(stderr, "ACOUSTIC - Fatal error!\n");
            fprintf(stderr, "  Solver type is not set = \"%s\".\n", solver);
            exit(1);
        }
       
        for (i = 0; i < node_num; i++)
        {
            energy[it] += u_now[i]*u_now[i];
        }
		if ((it + 1) % 10 == 0)
		{
			for (i = 0; i < node_num; i++)
			{
				fprintf(fp_u, "%f	", u_now[i]);
			}
			fprintf(fp_u, "\n");
		}
		fprintf(fp_energy, "%f\n", energy[it]);
	}
    
    free(mass_coo_i);
    free(mass_coo_j);
    free(mass_coo_x);
    free(stif_coo_i);
    free(stif_coo_j);
    free(stif_coo_x);
    free(mass_csr_p);
    free(mass_csr_j);
    free(mass_csr_x);
    free(stif_csr_p);
    free(stif_csr_j);
    free(stif_csr_x);
    free(mass_csr_j_temp);
    free(mass_csr_x_temp);
    free(stif_csr_j_temp);
    free(stif_csr_x_temp);
    free(rho);
    free(vp);
	free(rhs);
	free(u_old);
	free(u_now);
	free(y1);
	free(y2);
	free(y3);
	free(energy);
	fclose(fp_u);
	fclose(fp_energy);
	printf("\n Acoutsic Normal End!\n");
}
