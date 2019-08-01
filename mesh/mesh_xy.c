#include "mesh_xy_t3.c"
#include "mesh_xy_t6.c"
#include "mesh_xy_t10.c"
#include "mesh_xy_q4.c"
#include "mesh_xy_q9.c"
#include "mesh_xy_q16.c"

int mesh_xy(char *type, int nelemx, int nelemy, int node_num, double xmin, double xmax, double ymin, double ymax, double **node_xy)

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
  int *element_node;

  if (strcmp(type, "T3") == 0)
  {
    mesh_xy_t3(nelemx, nelemy, node_num, xmin, xmax, ymin, ymax, node_xy);
  }
  else if (strcmp(type, "T6") == 0)
  {
    mesh_xy_t6(nelemx, nelemy, node_num, xmin, xmax, ymin, ymax, node_xy);
  }
  else if (strcmp(type, "T10") == 0)
  {
    mesh_xy_t10(nelemx, nelemy, node_num, xmin, xmax, ymin, ymax, node_xy);
  }
  else if (strcmp(type, "Q4") == 0)
  {
    mesh_xy_q4(nelemx, nelemy, node_num, xmin, xmax, ymin, ymax, node_xy);
  }
  else if (strcmp(type, "Q9") == 0)
  {
    mesh_xy_q9(nelemx, nelemy, node_num, xmin, xmax, ymin, ymax, node_xy);
  }
  else if (strcmp(type, "Q16") == 0)
  {
    mesh_xy_q16(nelemx, nelemy, node_num, xmin, xmax, ymin, ymax, node_xy);
  }
  else
  {
    node_xy = NULL;
    fprintf(stderr, "\n");
    fprintf(stderr, "MESH_XY - Fatal error!\n");
    fprintf(stderr, "  Illegal value of type = \"%s\".\n", type);
    exit(1);
  }
}
