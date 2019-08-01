#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../mesh/element_type.c"
#include "../mesh/mesh_node_num.c"
#include "../mesh/mesh_element.c"
#include "../mesh/mesh_element_num.c"
#include "../mesh/mesh_element_order.c"
#include "../mesh/mesh_xy.c"
#include "../assemble_sparse/mass_sparse_all.c"
#include "../assemble_sparse/stif_sparse_all.c"
#include "../model/acoustic_model.c"
#include "../model/elastic_model.c"
#include "../absorbing_boundary/absorbing_boundary_pml.c"
#include "../absorbing_boundary/absorbing_boundary_mpml.c"
#include "../coo2csr_lib/csr_matvec.c"
#include "../source_receiver/node_location.c"
#include "../source_receiver/seismic_source.c"

#include "../wave/acoustic.c"
int main()
/******************************************************************************/
/*
 *       MESH_TEST tests the all subroutines.
 *
  List:

    I  ELEMENT_TYPE   Definition
    -  ------------   ----------
    1  T3             3 node linear triangle;
    2  T6             6 node quadratic triangle;
    3  T10            10 node cubic triangle.
    4  Q4             4 node linear Lagrange/serendipity quadrilateral;
    5  Q9             9 node quadratic Lagrange quadrilateral;
    6  Q16            16 node cubic Lagrange quadrilateral;

*/
{
  
  /***************************************
          model mesh parameters  		
  ****************************************/
  int element,order,i,j;
  char *type;
  int type_code = 1;
  int nelemx    = 100;
  int nelemy    = 100;
  int edge_size = 10;
  int node_num      = 0;
  int element_num   = 0;
  int element_order = 0;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  int    *element_node = NULL;
  double **node_xy     = NULL;
 /***************************************
              seismic source		
  ****************************************/
  double source_x;
  double source_y;
  double f0;
  double t0;
  int 	 source_node;
  
  /***************************************
              time evolution
  ****************************************/
  int  step;
  double dt;
    
  /***************************************
      coo(i,j,x) and csr(p,j,x) arrays  		
  ****************************************/
  int     nnz        = 0;     // nonezero number
  int     csr_p_size = 0;     // csr size for p
  int     mass_csr_size = 0;  // csr size for mass: j and x
  int     stif_csr_size = 0;  // csr size for stif: j and x
  int     lumpflag   = 1;     // lump mass 
  int    *mass_coo_i = NULL;
  int    *mass_coo_j = NULL;
  double *mass_coo_x = NULL;
  int    *stif_coo_i = NULL;
  int    *stif_coo_j = NULL;
  double *stif_coo_x = NULL;
  int    *mass_csr_p = NULL;
  int 	 *mass_csr_j = NULL;
  double *mass_csr_x = NULL;
  int    *stif_csr_p = NULL;
  int 	 *stif_csr_j = NULL;
  double *stif_csr_x = NULL;
  int    *mass_csr_j_temp = NULL;
  double *mass_csr_x_temp = NULL;
  int    *stif_csr_j_temp = NULL;
  double *stif_csr_x_temp = NULL;
  double *mass_lump = NULL;

  /***************************************
    model density and velocity parameters  		
  ****************************************/
  double vp_max;
  double *rho = NULL;
  double *vp  = NULL;
  double *vs  = NULL;
  double **c  = NULL;

  /***************************************
    abosorbing bc, pml and mpml parameters		
  ****************************************/
  int pml_nx;
  int pml_ny;
  int use_pml_xmin = 1;
  int use_pml_xmax = 1;
  int use_pml_ymin = 1;
  int use_pml_ymax = 1;  // use_pml_ymax = 0, free surface on ymax.
  double *pml_dx  = NULL;
  double *pml_dy  = NULL;
  double *pml_dxx = NULL;
  double *pml_dyy = NULL;
  double *mpml_dx = NULL;
  double *mpml_dy = NULL;
  double *mpml_dxx = NULL;
  double *mpml_dyy = NULL;
  double *mpml_dxx_pyx = NULL;
  double *mpml_dyy_pxy = NULL;
  
  /***************************************
              file pointers
  ****************************************/
  FILE  *fp_node_xy;
  
  /***************************************
            prepare parameters
  ****************************************/
  xmin = 0.0;  
  ymin = 0.0;
  xmax = edge_size * nelemx;
  ymax = edge_size * nelemy; 
  source_x = xmax / 2.0;
  source_y = ymax / 2.0;
  f0   = 20;
  t0   = 1.2 / f0;
  dt   = 0.0005;
  step = 1000;
  
  pml_nx = 20;
  pml_ny = 20;
  vp_max = 2000.0;
  
  /***************************************
    select the element type and mesh model
  ****************************************/
  type 		= element_type 	    ( type_code            );
  node_num      = mesh_node_num     ( type, nelemx, nelemy );
  element_num   = mesh_element_num  ( type, nelemx, nelemy );
  element_node  = mesh_element      ( type, nelemx, nelemy );
  element_order = mesh_element_order( type                 );
  csr_p_size 	= node_num + 1;
  nnz		= element_num * element_order * element_order;
  printf ( "\n Test all \n");
  printf ( "\n Element type is     %s\n", type          );
  printf ( "\n Node number is      %d\n", node_num      );
  printf ( "\n Element number is   %d\n", element_num   );
  printf ( "\n Element order is    %d\n", element_order );
  printf ( "\n csr_p_size is       %d\n", csr_p_size    );
  printf ( "\n None zero number is %d\n", nnz 	        );
  
  /***************************************
        allocate the dynamic arrays
  ****************************************/
  printf ( "\n Allocate the dynamic arrays. \n");
  node_xy   = (double **)malloc( 2   * sizeof( double ) );  
  for(i=0;i<2;i++) node_xy[i]=malloc(sizeof(double)*node_num);
  // mass lump
  mass_lump       = (double  *)malloc(node_num * sizeof( double ) );
  // coo (i,j,x), mass and stif matrices 
  mass_coo_i = (int     *)malloc( nnz * sizeof( int    ) );
  mass_coo_j = (int     *)malloc( nnz * sizeof( int    ) );
  mass_coo_x = (double  *)malloc( nnz * sizeof( double ) );
  stif_coo_i = (int     *)malloc( nnz * sizeof( int    ) );
  stif_coo_j = (int     *)malloc( nnz * sizeof( int    ) );
  stif_coo_x = (double  *)malloc( nnz * sizeof( double ) );
  // csr (p,j,x), mass and stif matrices 
  mass_csr_p = (int     *)malloc( csr_p_size * sizeof( int    ) );
  stif_csr_p = (int     *)malloc( csr_p_size * sizeof( int    ) );
  mass_csr_j_temp = (int     *)malloc( nnz * sizeof( int    ) );
  mass_csr_x_temp = (double  *)malloc( nnz * sizeof( double ) );
  stif_csr_j_temp = (int     *)malloc( nnz * sizeof( int    ) );
  stif_csr_x_temp = (double  *)malloc( nnz * sizeof( double ) );
  // model density and velocity matrices
  rho      = (double  *)malloc( node_num * sizeof( double ) );
  vp       = (double  *)malloc( node_num * sizeof( double ) );
  vs       = (double  *)malloc( node_num * sizeof( double ) );
  c        = (double **)malloc( 4   * sizeof( double ) );  
  for(i = 0; i < 4; i++) c[i] = malloc( node_num * sizeof(double) );
  
  // pml and mpml matrices
  pml_dx   = (double  *)malloc( node_num * sizeof( double ) );
  pml_dy   = (double  *)malloc( node_num * sizeof( double ) ); 
  pml_dxx  = (double  *)malloc( node_num * sizeof( double ) );
  pml_dyy  = (double  *)malloc( node_num * sizeof( double ) );
  mpml_dx  = (double  *)malloc( node_num * sizeof( double ) );
  mpml_dy  = (double  *)malloc( node_num * sizeof( double ) ); 
  mpml_dxx = (double  *)malloc( node_num * sizeof( double ) );
  mpml_dyy = (double  *)malloc( node_num * sizeof( double ) );
  mpml_dxx_pyx = (double  *)malloc( node_num * sizeof( double ) );
  mpml_dyy_pxy = (double  *)malloc( node_num * sizeof( double ) );
  
  /***************************************
            call functions
  ****************************************/
  // set (x,y) for every node 
  mesh_xy( type, nelemx, nelemy, node_num, xmin, xmax, ymin, ymax, node_xy);
  fp_node_xy = fopen("node_xy.dat", "w");
  for(i=0;i<node_num;i++) fprintf( fp_node_xy, "%f	%f\n", node_xy[0][i], node_xy[1][i]);
  
  // set acoustic model 
  acoustic_model( node_num, element_num, element_order, element_node, node_xy, rho, vp);
  // set elastic model
  elastic_model ( node_num, element_num, element_order, element_node, node_xy, rho, c );
  // set pml abosorbing profiles
  absorbing_boundary_pml ( node_num, element_num, element_order, element_node, node_xy, pml_nx, pml_ny, edge_size, xmin, xmax, ymin, ymax, vp_max, \
			   use_pml_xmin, use_pml_xmax, use_pml_ymin, use_pml_ymax, pml_dx, pml_dy, pml_dxx, pml_dyy );
  // set mpml abosorbing profiles
  absorbing_boundary_mpml( node_num, element_num, element_order, element_node, node_xy, pml_nx, pml_ny, edge_size, xmin, xmax, ymin, ymax, vp_max, \
			   use_pml_xmin, use_pml_xmax, use_pml_ymin, use_pml_ymax, mpml_dx, mpml_dy, mpml_dxx, mpml_dyy, mpml_dxx_pyx, mpml_dyy_pxy );			  
  // get lump mass matrix 
  mass_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, rho, mass_coo_i, mass_coo_j, mass_coo_x, 1);
  for( i = 0; i < node_num; i++ ) mass_lump[i] = mass_coo_x[i];
  // get mass and stif matrices in coo format 
  mass_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, rho, mass_coo_i, mass_coo_j, mass_coo_x, 0);			  
  stif_sparse_all(type, node_num, element_num, element_order, element_node, node_xy, vp , stif_coo_i, stif_coo_j, stif_coo_x   );
  // call fortran subroutines to convert coo to csr. Attention: mi, mj, mass, Mp_temp, Mj_temp, Mass_temp arethe address, 
  // do not need use & to get the address of the pointers.
  __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &mass_csr_size, mass_coo_i, mass_coo_j, mass_coo_x, mass_csr_p, mass_csr_j_temp, mass_csr_x_temp ); 
  __coo2csr_lib_MOD_coo2csr_canonical(&nnz, &csr_p_size, &stif_csr_size, stif_coo_i, stif_coo_j, stif_coo_x, stif_csr_p, stif_csr_j_temp, stif_csr_x_temp ); 
  // allocate the real size mass and stif csr matrices
  mass_csr_j = (int      *)malloc( mass_csr_size * sizeof( int    ) );
  mass_csr_x = (double   *)malloc( mass_csr_size * sizeof( double ) );
  stif_csr_j = (int      *)malloc( stif_csr_size * sizeof( int    ) );
  stif_csr_x = (double   *)malloc( stif_csr_size * sizeof( double ) );
  
  for( i = 0; i < mass_csr_size; i++)
  {
     mass_csr_j[i] = mass_csr_j_temp[i];
     mass_csr_x[i] = mass_csr_x_temp[i];
  }
  for( i = 0; i < stif_csr_size; i++)
  {
     stif_csr_j[i] = stif_csr_j_temp[i];
     stif_csr_x[i] = stif_csr_x_temp[i];
  }   
  
  source_node = node_location( node_num, edge_size, node_xy, source_x, source_y );
  printf("\n Source node is %d\n",source_node);
  acoustic( node_num, step, dt, f0, t0, source_node, mass_lump, mass_csr_p, mass_csr_j, mass_csr_x, stif_csr_p, stif_csr_j, stif_csr_x);
  
  
  
  /***************************************
              free memory
  ****************************************/
  free ( type );
  free ( element_node );
  for(i=0;i<2;i++) free(node_xy[i]); 
  free(node_xy); 
  free(mass_coo_i);
  free(mass_coo_j);
  free(mass_coo_x);
  free(stif_coo_i);
  free(stif_coo_j);
  free(stif_coo_x);
  free(mass_csr_p);
  free(mass_csr_j_temp);
  free(mass_csr_x_temp);
  free(mass_csr_j);
  free(mass_csr_x);
  free(stif_csr_p);
  free(stif_csr_j_temp);
  free(stif_csr_x_temp);
  free(stif_csr_j);
  free(stif_csr_x);
  free(rho);
  free(vp);
  free(vs);
  for(i = 0; i < 4; i++) free(c[i]); 
  free(c); 
  
  free(pml_dx);
  free(pml_dy);
  free(pml_dxx);
  free(pml_dyy);
  free(mpml_dx);
  free(mpml_dy);
  free(mpml_dxx);
  free(mpml_dyy);
  free(mpml_dxx_pyx);
  free(mpml_dyy_pxy);
  
  fclose(fp_node_xy);

  printf("\n Free all dynamic arrays!\n");
  printf("\n Normal End!\n");

  return 0;

}
