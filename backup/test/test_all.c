#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "../mesh/element_type.c"
#include "../mesh/mesh_element_order.c"
#include "../mesh/mesh_node_num.c"
#include "../mesh/mesh_node_num_read.c"
#include "../mesh/mesh_element_num.c"
#include "../mesh/mesh_element_num_read.c"
#include "../mesh/mesh_element.c"
#include "../mesh/mesh_element_read.c"
#include "../mesh/mesh_xy.c"
#include "../mesh/mesh_xy_read.c"
#include "../assemble/mass_sparse_all.c"
#include "../assemble/stif_sparse_all.c"
#include "../model/elastic_model.c"
#include "../absorbing_boundary/absorbing_boundary_mpml.c"
#include "../sparse_matrix/csr_matvec.c"
#include "../source_receiver/node_location.c"
#include "../source_receiver/seismic_source.c"
#include "../solver/solver_type.c"
#include "../solver/mgmres/mgmres.c"
#include "../wave/elastic_mpml.c"

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
    1  mgmres         Generalized Minimum Residual (GMRES) algorithm, CSR format;
    2  superlu        LU decomposition to solver linear system;
    3  masslump       masslump.
*/
{

  /***************************************
          model mesh parameters  		
  ****************************************/
  int i, j;
  char *type;
  char *solver;
  int type_code = 1;
  int solver_code = 3;
  int nelemx = 400;
  int nelemy = 400;
  int edge_size = 5;
  int pml_nx = 10;
  int pml_ny = 10;
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

  /***************************************
              seismic source		
  ****************************************/
  double source_x;
  double source_y;
  double f0 = 10.0;
  double t0 = 1.2 / 10.0;
  int source_node;

  /***************************************
              time evolution
  ****************************************/
  int step = 2000;
  double dt = 0.0005;

  /***************************************
              file pointers
  ****************************************/
  FILE *fp_element_node;
  FILE *fp_node_xy;
  char element_name[30] = "elements_400x400.txt";
  char node_xy_name[30] = "nodes_400x400.txt";
  char mesh_par_name[30] = "mesh_par_400x400.txt";
  /***************************************
            prepare parameters
  ****************************************/
  program_start_time = omp_get_wtime();
  xmin = 0.0;
  ymin = 0.0;
  xmax = edge_size * nelemx;
  ymax = edge_size * nelemy;
  source_x = xmax / 2.0;
  source_y = ymax / 2.0;

  /***************************************
    select the element type and mesh model
  ****************************************/
  type = element_type(type_code);
  solver = solver_type(solver_code);
  element_order = mesh_element_order(type);
  //node_num = mesh_node_num(type, nelemx, nelemy);
  //element_num = mesh_element_num(type, nelemx, nelemy);
  node_num = mesh_node_num_read(mesh_par_name);
  element_num = mesh_element_num_read(mesh_par_name);
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
  element_node = (int *)malloc(element_order * element_num * sizeof(int));
  node_xy = (double **)malloc(2 * sizeof(double));
  for (i = 0; i < 2; i++)
    node_xy[i] = malloc(sizeof(double) * node_num);

  mesh_element_read(element_name, element_num, element_order, element_node);
  mesh_xy_read(node_xy_name, node_num, node_xy);
  //mesh_element(type, nelemx, nelemy, element_node);
  //mesh_xy(type, nelemx, nelemy, node_num, xmin, xmax, ymin, ymax, node_xy);

  source_node = node_location(node_num, edge_size, node_xy, source_x, source_y);

  elastic_mpml(type, node_num, element_num, element_order, element_node, node_xy, nnz, csr_p_size, step, dt, f0, t0, edge_size, xmin, xmax, ymin, ymax, pml_nx, pml_ny, source_node, solver);

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

  printf("\n Free all dynamic arrays!\n");
  program_run_time = omp_get_wtime() - program_start_time;
  printf("\n Total run time is: %f\n", program_run_time);
  printf("\n Normal End!\n");
  return 0;
}
