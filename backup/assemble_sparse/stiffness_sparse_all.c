#include "stiffness_sparse_t3.c"
#include "stiffness_sparse_t6.c"
#include "stiffness_sparse_t10.c"
#include "stiffness_sparse_q4.c"
#include "stiffness_sparse_q9.c"
#include "stiffness_sparse_q16.c"

void stiffness_sparse_all(char *type, int node_num, int element_num, int element_order,  int *element_node, double **node_xy, int *stiffi, int *stiffj, double *stiffness)

/******************************************************************************/
/*
  Purpose:

   stiffness_sparse_all computes the stiffness matrix, store the matrix in the coo format (i,j,value), according to the element type

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

  if ( strcmp(type,"T3") == 0 )
  {
    stiffness_sparse_t3(node_num, element_num, element_order, element_node, node_xy, stiffi, stiffj, stiffness);
  }
  else if ( strcmp(type,"T6") == 0 )
  {
    stiffness_sparse_t6(node_num, element_num, element_order, element_node, node_xy, stiffi, stiffj, stiffness);
  }
  else if ( strcmp(type,"T10") == 0 )
  {
    stiffness_sparse_t10(node_num, element_num, element_order, element_node, node_xy, stiffi, stiffj, stiffness);
  }
  else if ( strcmp(type,"Q4") == 0 )
  {
    stiffness_sparse_q4(node_num, element_num, element_order, element_node, node_xy, stiffi, stiffj, stiffness);
  }
  else if ( strcmp(type,"Q9") == 0 )
  {
    stiffness_sparse_q9(node_num, element_num, element_order, element_node, node_xy, stiffi, stiffj, stiffness);
  }
  else if ( strcmp(type,"Q16") == 0 )
  {
    stiffness_sparse_q16(node_num, element_num, element_order, element_node, node_xy, stiffi, stiffj, stiffness);
  }
  else
  {
    element_node = NULL;
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "STIFFNESS_SPARSE_ALL - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of type = \"%s\".\n", type );
    exit ( 1 );
  }

}
