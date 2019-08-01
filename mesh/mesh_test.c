#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "element_type.c"
#include "mesh_element_num.c"
#include "mesh_element_order.c"
#include "mesh_node_num.c"
#include "mesh_element.c"
void main()

/******************************************************************************/
/*
  Purpose: 

    MESH_TEST tests the mesh routines.

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
  char *type;
  int *element_node;
  int element_num = 0;
  int element_order = 0;
  int node_num = 0;
  int nelemx;
  int nelemy;
  int i, j;
  /*
  NODE is defined as a vector rather than a two dimensional array,
  so that we can handle the various cases using a single array.
*/
  type = element_type(1);
  nelemx = 3;
  nelemy = 2;

  element_num = mesh_element_num(type, nelemx, nelemy);
  element_order = mesh_element_order(type);
  node_num = mesh_node_num(type, nelemx, nelemy);
  element_node = mesh_element(type, nelemx, nelemy);

  printf("mesh_TEST: Test the mesh routine for element %s\n", type);
  printf("Element number is %d\n", element_num);
  printf("Node number is %d\n", node_num);
  printf("Element order is %d\n", element_order);
  for (i = 0; i < element_order * element_num; i = i + element_order)
  {
    for (j = 0; j < element_order; j++)
    {
      printf("%4d", element_node[i + j]);
    }
    printf("\n");
  }

  free(type);
  free(element_node);

  return;
}
