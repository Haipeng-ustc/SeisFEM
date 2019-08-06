#include "../solver/mgmres/mgmres.c"

void acoustic(int node_num, int step, double dt, double f0, double t0, int source_node, double *rho, double *vp, double *mass_lump, int mass_csr_size,
			  int *mass_csr_p, int *mass_csr_j, double *mass_csr_x, int *stif_csr_p, int *stif_csr_j, double *stif_csr_x, char *solver)

{
	int i, it, csr_p_size;
	double time;
	double point_source;
	double *rhs = NULL;
	double *u_old = NULL;
	double *u_now = NULL;
    double *u_now_temp = NULL;

	double *y1 = NULL;
	double *y2 = NULL;
	double *y3 = NULL;
	double *energy = NULL;
    // for mgmres solver
    double tol_abs = 1.0e-08; // a relative tolerance comparing the current residual to the initial residual.
    double tol_rel = 1.0e-08; // an absolute tolerance applied to the current residual.
    int    itr_max = 100;     // the maximum number of (outer) iterations to take.
    int    mr      = 50;      // the maximum number of (inner) iterations to take.
    double sum = 0.0;
	FILE *fp_u, *fp_energy;

    
    
    
    
    
    
    
    
	csr_p_size = node_num + 1;
	rhs = (double *)malloc(node_num * sizeof(double));
	u_old = (double *)malloc(node_num * sizeof(double));
	u_now = (double *)malloc(node_num * sizeof(double));
    u_now_temp = (double *)malloc(node_num * sizeof(double));

	y1 = (double *)malloc(node_num * sizeof(double));
	y2 = (double *)malloc(node_num * sizeof(double));
	y3 = (double *)malloc(node_num * sizeof(double));
	energy = (double *)malloc(node_num * sizeof(double));


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
	free(rhs);
	free(u_old);
	free(u_now);
    free(u_now_temp);

	free(y1);
	free(y2);
	free(y3);
	free(energy);
	fclose(fp_u);
	fclose(fp_energy);
	printf("\n Acoutsic Normal End!\n");
}
