#include "../shape/shape_t3.c"
#include "../shape/shape_t6.c"
#include "../shape/shape_t10.c"
#include "../shape/shape_q4.c"
#include "../shape/shape_q9.c"
#include "../shape/shape_q16.c"
#include "mass_sparse_t3.c"
#include "mass_sparse_t6.c"
#include "mass_sparse_t10.c"
#include "mass_sparse_q4.c"
#include "mass_sparse_q9.c"
#include "mass_sparse_q16.c"
void mass_sparse_all(char *type, int node_num, int element_num, int element_order, int *element_node, double **node_xy,
                     int *mass_coo_i, int *mass_coo_j, double *mass_coo_x, int lumpflag)

/******************************************************************************/
/*
  Purpose:

   mass_sparse_all computes the mass matrix, store the matrix in the coo format (i, j, x), according to the element type

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

  if (strcmp(type, "T3") == 0)
  {
    mass_sparse_t3(node_num, element_num, element_order, element_node, node_xy, mass_coo_i, mass_coo_j, mass_coo_x, lumpflag);
  }
  else if (strcmp(type, "T6") == 0)
  {
    mass_sparse_t6(node_num, element_num, element_order, element_node, node_xy, mass_coo_i, mass_coo_j, mass_coo_x, lumpflag);
  }
  else if (strcmp(type, "T10") == 0)
  {
    mass_sparse_t10(node_num, element_num, element_order, element_node, node_xy, mass_coo_i, mass_coo_j, mass_coo_x, lumpflag);
  }
  else if (strcmp(type, "Q4") == 0)
  {
    mass_sparse_q4(node_num, element_num, element_order, element_node, node_xy, mass_coo_i, mass_coo_j, mass_coo_x, lumpflag);
  }
  else if (strcmp(type, "Q9") == 0)
  {
    mass_sparse_q9(node_num, element_num, element_order, element_node, node_xy, mass_coo_i, mass_coo_j, mass_coo_x, lumpflag);
  }
  else if (strcmp(type, "Q16") == 0)
  {
    mass_sparse_q16(node_num, element_num, element_order, element_node, node_xy, mass_coo_i, mass_coo_j, mass_coo_x, lumpflag);
  }
  else
  {
    element_node = NULL;
    fprintf(stderr, "\n");
    fprintf(stderr, "MASS_SPARSE_ALL - Fatal error!\n");
    fprintf(stderr, "  Illegal value of type = \"%s\".\n", type);
    exit(1);
  }
}
