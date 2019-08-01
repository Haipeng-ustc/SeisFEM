#include "mass_t3.c"
#include "mass_t6.c"
#include "mass_t10.c"
#include "mass_q4.c"
#include "mass_q9.c"
#include "mass_q16.c"

void mass_all(char *type, int node_num, int element_num, int element_order,  int *element_node, double **node_xy, double **mass)

/******************************************************************************/
/*
  Purpose:

   mass_all computes the mass matrix, according to the element type

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
    mass_t3(node_num, element_num, element_order, element_node, node_xy, mass);
  }
  else if ( strcmp(type,"T6") == 0 )
  {
    mass_t6(node_num, element_num, element_order, element_node, node_xy, mass);
  }
  else if ( strcmp(type,"T10") == 0 )
  {
    mass_t10(node_num, element_num, element_order, element_node, node_xy, mass);
  }
  else if ( strcmp(type,"Q4") == 0 )
  {
    mass_q4(node_num, element_num, element_order, element_node, node_xy, mass);
  }
  else if ( strcmp(type,"Q9") == 0 )
  {
    mass_q9(node_num, element_num, element_order, element_node, node_xy, mass);
  }
  else if ( strcmp(type,"Q16") == 0 )
  {
    mass_q16(node_num, element_num, element_order, element_node, node_xy, mass);
  }
  else
  {
    element_node = NULL;
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "MASS_ALL - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of type = \"%s\".\n", type );
    exit ( 1 );
  }


}
