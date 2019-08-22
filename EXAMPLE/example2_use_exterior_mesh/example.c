#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "../../mesh/element_type.c"
#include "../../mesh/mesh_element_order.c"
#include "../../mesh/mesh_node_num.c"
#include "../../mesh/mesh_element_num.c"
#include "../../mesh/mesh_element.c"
#include "../../mesh/mesh_xy.c"
#include "../../assemble/mass_sparse_all.c"
#include "../../assemble/stif_sparse_all.c"
#include "../../model_elastic_parameter/model_elastic_parameter.c"
#include "../../pml/abc_mpml.c"
#include "../../sparse_matrix/csr_matvec.c"
#include "../../source_receiver/node_location.c"
#include "../../source_receiver/set_receiver_node.c"
#include "../../source_receiver/set_source_node.c"
#include "../../source_receiver/seismic_source.c"
#include "../../source_receiver/set_rec_and_src.c"
#include "../../solver/solver_type.c"
#include "../../solver/pardiso/pardiso_unsym.c"
#include "../../solver/mgmres/mgmres.c"
#include "../../time_evolution/elastic_wave.c"

int main()
/******************************************************************************/
/*
 *       TEST_ALL tests the all subroutines.
 *
 
    ELEMENT_TYPE List:
    I  ELEMENT_TYPE   Definition
    -  ------------   ----------
    1  T3             3 node linear triangle;
    2  T6             6 node quadratic triangle;
    3  T10            10 node cubic triangle.
    4  Q4             4 node linear Lagrange/serendipity quadrilateral;
    5  Q9             9 node quadratic Lagrange quadrilateral;
    6  Q16            16 node cubic Lagrange quadrilateral.
 
    SOLVER_TYPE List:
    I  SOLVER_TYPE   Definition
    -  ------------   ----------
    1  pardiso        require license and only valid for username haipeng;
    2  mgmres         Generalized Minimum Residual (GMRES) algorithm, CSR format;
    3  masslump       masslump.
*/
{

  /***************************************
          model mesh parameters  		
  ****************************************/
  int i, j;
  char *type;
  char *solver;
  int type_code = 0;
  int solver_code = 0;
  int nelemx = 0;
  int nelemy = 0;
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
  double edge_size = 0.0;
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
  double src_x_first;
  double src_y_first;
  double src_x_last;
  double src_y_last;
  double rec_x_first;
  double rec_y_first;
  double rec_x_last;
  double rec_y_last;
  double *src_x;
  double *src_y;
  double *rec_x;
  double *rec_y;
  double x_temp, y_temp;
  /***************************************
              time evolution
  ****************************************/
  int step = 0;
  double dt = 0.0;

  /***************************************
              file pointers
  ****************************************/
  FILE *fp_par;
  FILE *fp_element_node;
  FILE *fp_node_xy;
  FILE *fp_num;
  /***************************************
            prepare parameters
  ****************************************/
  program_start_time = omp_get_wtime();

  if ((fp_par = fopen("par.txt", "r")) == NULL)
    printf("\n can not open par.txt file\n");
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

  /***************************************
    select the element type and mesh model
  ****************************************/
  type = element_type(type_code);
  solver = solver_type(solver_code);
  element_order = mesh_element_order(type);

  if (use_exterior_mesh == 1)
  {
    // use exterior mesh and read the mesh scheme


    if ((fp_num = fopen("./mesh_exterior/element_and_node_num.txt", "r")) == NULL)
      printf("\n element_and_node_num.txt file cannot open\n");
    if ((fp_element_node = fopen("./mesh_exterior/element_nodes.txt", "r")) == NULL)
      printf("\n element file cannot open\n");
    if ((fp_node_xy = fopen("./mesh_exterior/nodes_xy.txt", "r")) == NULL)
      printf("\n node_xy file cannot open\n");
    
    fscanf(fp_num, "element_num = %d\n", &element_num);
    fscanf(fp_num, "node_num = %d\n", &node_num);
    fclose(fp_num);

    element_node = (int *)malloc(element_order * element_num * sizeof(int));
    node_xy = (double **)malloc(2 * sizeof(double));
    for (i = 0; i < 2; i++)
      node_xy[i] = malloc(sizeof(double) * node_num);

    // read element_node
    for (i = 0; i < element_num * element_order; i = i + element_order)
    {
      for (j = 0; j < element_order; j++)
      {
        fscanf(fp_element_node, "%d ", &element_node[i + j]);
      }
      fscanf(fp_element_node, "\n");
    }
    // read node_xy
    for (i = 0; i < node_num; i++)
    {
      fscanf(fp_node_xy, "%lf %lf\n", &x_temp, &y_temp);
      node_xy[0][i] = x_temp;
      node_xy[1][i] = y_temp;
    }
    fclose(fp_element_node);
    fclose(fp_node_xy);
  }
  else
  {
    // use internal mesh and save the mesh scheme
    xmin = 0.0;
    ymin = 0.0;
    xmax = edge_size * nelemx;
    ymax = edge_size * nelemy;
    node_num = mesh_node_num(type, nelemx, nelemy);
    element_num = mesh_element_num(type, nelemx, nelemy);
    
    element_node = (int *)malloc(element_order * element_num * sizeof(int));
    node_xy = (double **)malloc(2 * sizeof(double));
    for (i = 0; i < 2; i++)
      node_xy[i] = malloc(sizeof(double) * node_num);

    mesh_element(type, nelemx, nelemy, element_node);
    mesh_xy(type, nelemx, nelemy, node_num, xmin, xmax, ymin, ymax, node_xy);
    fp_element_node = fopen("./mesh_internal/element_node.txt", "w");
    fp_node_xy = fopen("./mesh_internal/node_xy.txt", "w");
    for (i = 0; i < element_num * element_order; i = i + element_order)
    {
      for (j = 0; j < element_order; j++)
        fprintf(fp_element_node, "%d  ", element_node[i + j]);
      fprintf(fp_element_node, "\n");
    }

    for (i = 0; i < node_num; i++)
    {
      fprintf(fp_node_xy, "%f  %f\n", node_xy[0][i], node_xy[1][i]);
    }
    fclose(fp_element_node);
    fclose(fp_node_xy);
  }

  csr_p_size = node_num + 1;
  nnz = element_num * element_order * element_order;

  printf("\n Test all \n");
  printf("\n Element type is     %s\n", type);
  printf("\n Node number is      %d\n", node_num);
  printf("\n Element number is   %d\n", element_num);
  printf("\n Element order is    %d\n", element_order);
  printf("\n csr_p_size is       %d\n", csr_p_size);
  printf("\n None zero number is %d\n", nnz);
  printf("\n solver is           %s\n", solver);
  printf("\n Element type is     %s\n", type);
  printf("\n Node number is      %d\n", node_num);

  rec_node = (int *)malloc(rec_num * sizeof(int));
  rec_x = (double *)malloc(rec_num * sizeof(double));
  rec_y = (double *)malloc(rec_num * sizeof(double));
  src_node = (int *)malloc(src_num * sizeof(int));
  src_x = (double *)malloc(src_num * sizeof(double));
  src_y = (double *)malloc(src_num * sizeof(double));

  set_rec_and_src(rec_num, rec_x_first, rec_y_first, rec_x_last, rec_y_last, rec_x, rec_y, src_num, src_x_first, src_y_first, src_x_last, src_y_last, src_x, src_y);
  set_receiver_node(rec_num, node_num, edge_size, rec_x, rec_y, node_xy, rec_node);
  set_source_node(src_num, node_num, edge_size, src_x, src_y, node_xy, src_node);

  elastic_wave(type, node_num, element_num, element_order, element_node, node_xy, nnz, csr_p_size, step, dt, f0, t0, edge_size, \
               xmin, xmax, ymin, ymax, pml_nx, pml_ny, src_num, src_node, rec_num, rec_node, solver, use_exterior_mesh);

  /***************************************
              free memory
  ****************************************/
  fclose(fp_element_node);
  fclose(fp_node_xy);
  free(type);
  free(solver);
  for (i = 0; i < 2; i++)
    free(node_xy[i]);
  free(node_xy);
  free(rec_node);
  free(rec_x);
  free(rec_y);
  free(src_node);
  free(src_x);
  free(src_y);
  program_run_time = omp_get_wtime() - program_start_time;
  printf("\n Total run time is: %f\n", program_run_time);
  printf("\n Example Normal End!\n");
  return 0;
}
