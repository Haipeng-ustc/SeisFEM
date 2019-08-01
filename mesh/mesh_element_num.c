
int mesh_element_num ( char *code, int nelemx, int nelemy)

/******************************************************************************/
/*
  Purpose:

    GRID_ELEMENT_NUM returns the number of elements.

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
  int element_num;

  if ( strcmp(code,"T3") == 0 )
  {
    element_num = 2 * nelemx * nelemy;
  }
  else if ( strcmp(code,"T6") == 0 )
  {
    element_num = 2 * nelemx * nelemy;
  }
  else if ( strcmp(code,"T10") == 0 )
  {
    element_num = 2 * nelemx * nelemy;
  }
  else if ( strcmp(code,"Q4") == 0 )
  {
    element_num = nelemx * nelemy;
  }
  else if ( strcmp(code,"Q9") == 0 )
  {
    element_num = nelemx * nelemy;
  }
  else if ( strcmp(code,"Q16") == 0 )
  {
    element_num = nelemx * nelemy;
  }
  else
  {
    printf ("\n" );
    printf ("MESH_ELEMENT_NUM - Fatal error!\n" );
    printf ("  Illegal value of CODE = \"%s\".\n", code );
    element_num = -1;
    exit ( 1 );
  }

  return element_num;
}
