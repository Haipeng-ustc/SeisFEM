#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#include <string.h>

int main()
{

    int i, j;
    char *type;
    char *solver;
    int type_code = 0;
    int solver_code = 0;
    int nelemx = 0;
    int nelemy = 0;
    double edge_size = 0.0;
    int pml_nx = 20;
    int pml_ny = 20;
    int node_num = 0;
    int element_num = 0;
    int element_order = 0;
    int nnz = 0;
    int csr_p_size = 0;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    int *element_node = NULL;
    double **node_xy = NULL;
    double program_start_time, program_run_time;
    int use_exterior_mesh = 0;
    /***************************************
        seismic source and receiver	
  ****************************************/
    double f0 = 0.0;
    double t0 = 0.0;
    int rec_num;
    int src_num;
    int *src_node;
    int *rec_node;
    double src_x;
    double src_y;
    double rec_x_first;
    double rec_y_first;
    double rec_x_last;
    double rec_y_last;
     double src_x_first;
  double src_y_first;
  double src_x_last;
  double src_y_last;
  
    double *rec_x;
    double *rec_y;
    /***************************************
              time evolution
  ****************************************/
    int step = 0;
    double dt = 0.0;
    FILE *fp_par;
    fp_par = fopen("./par/par.txt", "r");
    fscanf(fp_par, "## Mesh parameters\n");
    fscanf(fp_par, "use_exterior_mesh = %d\n", &use_exterior_mesh);
    fscanf(fp_par, "type_code = %d\n", &type_code);
    fscanf(fp_par, "nelemx = %d\n", &nelemx);
    fscanf(fp_par, "nelemy = %d\n", &nelemy);
    fscanf(fp_par, "edge_size = %lf\n", &edge_size);
    fscanf(fp_par, "## source parameters\n");
    fscanf(fp_par, "f0 = %lf\n", &f0);
    fscanf(fp_par, "t0 = %lf\n", &t0);
    fscanf(fp_par, "src_num = %d\n", &src_num);
    fscanf(fp_par, "src_x_first = %lf\n", &src_x_first);
    fscanf(fp_par, "src_y_first = %lf\n", &src_y_first);
    fscanf(fp_par, "src_x_last  = %lf\n", &src_x_last);
    fscanf(fp_par, "src_y_last  = %lf\n", &src_y_last);
    fscanf(fp_par, "## receiver parameters\n");
    fscanf(fp_par, "rec_num = %d\n", &rec_num);
    fscanf(fp_par, "rec_x_first = %lf\n", &rec_x_first);
    fscanf(fp_par, "rec_y_first = %lf\n", &rec_y_first);
    fscanf(fp_par, "rec_x_last  = %lf\n", &rec_x_last);
    fscanf(fp_par, "rec_y_last  = %lf\n", &rec_y_last);
    fscanf(fp_par, "## time evolution parameters\n");
    fscanf(fp_par, "dt = %lf\n", &dt);
    fscanf(fp_par, "step = %d\n", &step);
    fscanf(fp_par, "solver_code = %d\n", &solver_code);
    fclose(fp_par);

    printf("%d\n", use_exterior_mesh);
    printf("%d\n", type_code);
    printf("%d\n", nelemx);
    printf("%d\n", nelemy);
    printf("%f\n", edge_size);
    printf("%f\n", f0);
    printf("%f\n", t0);
    printf("%d\n", src_num);
    printf("%f\n", src_x_first);
    printf("%f\n", src_y_first);
    printf("%f\n", src_x_last);
    printf("%f\n", src_y_last);

    printf("%d\n", rec_num);
    printf("%f\n", rec_x_first);
    printf("%f\n", rec_y_first);
    printf("%f\n", rec_x_last);
    printf("%f\n", rec_y_last);

    printf("%f\n", dt);
    printf("%d\n", step);
    printf("%d\n", solver_code);
}
