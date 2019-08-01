
int mesh_element_order(char *code)
/******************************************************************************/
/*
  Purpose:

    mesh_element_order returns the element order.
   
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
  int element_order;
  if (strcmp(code, "T3") == 0)
  {
    element_order = 3;
  }
  else if (strcmp(code, "T6") == 0)
  {
    element_order = 6;
  }
  else if (strcmp(code, "T10") == 0)
  {
    element_order = 10;
  }
  else if (strcmp(code, "Q4") == 0)
  {
    element_order = 4;
  }
  else if (strcmp(code, "Q9") == 0)
  {
    element_order = 9;
  }
  else if (strcmp(code, "Q16") == 0)
  {
    element_order = 16;
  }
  else
  {
    fprintf(stderr, "\n");
    fprintf(stderr, "MESH_ELEMENT_ORDER - Fatal error!\n");
    fprintf(stderr, "  Illegal value of CODE = \"%s\".\n", code);
    element_order = -1;
    exit(1);
  }

  return element_order;
}
