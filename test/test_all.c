#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "../mesh/element_type.c"
#include "../mesh/mesh_node_num.c"
#include "../mesh/mesh_element.c"
#include "../mesh/mesh_element_num.c"
#include "../mesh/mesh_element_order.c"
#include "../mesh/mesh_xy.c"
#include "../assemble/mass_sparse_all.c"
#include "../assemble/stif_sparse_all.c"
#include "../model/acoustic_model.c"
#include "../model/elastic_model.c"
#include "../absorbing_boundary/absorbing_boundary_pml.c"
#include "../absorbing_boundary/absorbing_boundary_mpml.c"
#include "../sparse_matrix/csr_matvec.c"
#include "../source_receiver/node_location.c"
#include "../source_receiver/seismic_source.c"
//#include "../wave/acoustic.c"
#include "../wave/acoustic_pml.c"
#include "../solver/solver_type.c"
#include "timestamp.c"

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
    6  Q16            16 node cubic Lagrange quadrilateral;
 
    SOLVER_TYPE List:
    I  SOLVER_TYPE   Definition
    -  ------------   ----------
    1  mgmres         Generalized Minimum Residual (GMRES) algorithm,
                      using compressed row sparse matrix format
    2  superlu        LU decomposition to solver linear system;
    3  masslump       masslump.
*/
{

  /***************************************
          model mesh parameters  		
  ****************************************/
  int element, order, i, j;
  char *type;
  char *solver;
  int type_code = 4;
  int solver_code = 3;
  int nelemx = 200;
  int nelemy = 200;
  int edge_size = 10;
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
  /***************************************
              seismic source		
  ****************************************/
  double source_x;
  double source_y;
  double f0;
  double t0;
  int source_node;

  /***************************************
              time evolution
  ****************************************/
  int step;
  double dt;

  /***************************************
              file pointers
  ****************************************/
  FILE *fp_node_xy;

  /***************************************
            prepare parameters
  ****************************************/
  timestamp ( );

  xmin = 0.0;
  ymin = 0.0;
  xmax = edge_size * nelemx;
  ymax = edge_size * nelemy;
  source_x = xmax / 2.0;
  source_y = ymax / 2.0;
  f0 = 20;
  t0 = 1.2 / f0;
  dt = 0.0005;
  step = 400;

  /***************************************
    select the element type and mesh model
  ****************************************/
  type          = element_type(type_code);
  solver        = solver_type(solver_code);
  node_num      = mesh_node_num(type, nelemx, nelemy);
  element_num   = mesh_element_num(type, nelemx, nelemy);
  element_node  = mesh_element(type, nelemx, nelemy);
  element_order = mesh_element_order(type);
  csr_p_size    = node_num + 1;
  nnz           = element_num * element_order * element_order;
  printf("\n Test all \n");
  printf("\n Element type is     %s\n", type);
  printf("\n Node number is      %d\n", node_num);
  printf("\n Element number is   %d\n", element_num);
  printf("\n Element order is    %d\n", element_order);
  printf("\n csr_p_size is       %d\n", csr_p_size);
  printf("\n None zero number is %d\n", nnz);
  printf("\n solver is           %s\n", solver);

  /***************************************
            call functions
  ****************************************/
  // set node_xy
  node_xy = (double **)malloc(2 * sizeof(double));
  for (i = 0; i < 2; i++)  node_xy[i] = malloc(sizeof(double) * node_num);

  mesh_xy(type, nelemx, nelemy, node_num, xmin, xmax, ymin, ymax, node_xy);
  fp_node_xy = fopen("node_xy.dat", "w");

  for (i = 0; i < node_num; i++) fprintf(fp_node_xy, "%f	%f\n", node_xy[0][i], node_xy[1][i]);
  // set source_node
  source_node = node_location(node_num, edge_size, node_xy, source_x, source_y);
  printf("\n Source node is      %d\n", source_node);

  //acoustic(type, node_num, element_num, element_order, element_node, node_xy, nnz, csr_p_size, step, dt, f0, t0, source_node, solver);
  acoustic_pml(type, node_num, element_num, element_order, element_node, node_xy, nnz, csr_p_size, step, dt, f0, t0, edge_size, \
               xmin, xmax, ymin, ymax, pml_nx, pml_ny, source_node, solver);

    
  /***************************************
              free memory
  ****************************************/
  free(type);
  free(element_node);
  for (i = 0; i < 2; i++) free(node_xy[i]);
  free(node_xy);

  fclose(fp_node_xy);

  printf("\n Free all dynamic arrays!\n");
  printf("\n Normal End!\n");
  timestamp ( );
  return 0;
}
