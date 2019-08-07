#include "mesh_element_t3.c"
#include "mesh_element_t6.c"
#include "mesh_element_t10.c"
#include "mesh_element_q4.c"
#include "mesh_element_q9.c"
#include "mesh_element_q16.c"

int mesh_element(char *type, int nelemx, int nelemy, int *element_node)

/******************************************************************************/
/*
  Purpose:

    MESH_ELEMENT returns the element mesh associated with any available element.

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
     mesh_t3_element(nelemx, nelemy, element_node);
  }
  else if (strcmp(type, "T6") == 0)
  {
     mesh_t6_element(nelemx, nelemy, element_node);
  }
  else if (strcmp(type, "T10") == 0)
  {
     mesh_t10_element(nelemx, nelemy, element_node);
  }
  else if (strcmp(type, "Q4") == 0)
  {
     mesh_q4_element(nelemx, nelemy, element_node);
  }
  else if (strcmp(type, "Q9") == 0)
  {
     mesh_q9_element(nelemx, nelemy, element_node);
  }
  else if (strcmp(type, "Q16") == 0)
  {
     mesh_q16_element(nelemx, nelemy, element_node);
  }
  else
  {
    element_node = NULL;
    fprintf(stderr, "\n");
    fprintf(stderr, "MESH_ELEMENT - Fatal error!\n");
    fprintf(stderr, "  Illegal value of type = \"%s\".\n", type);
    exit(1);
  }
}
