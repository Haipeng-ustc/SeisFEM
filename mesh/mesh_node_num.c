
int mesh_node_num ( char *code, int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose:

    MESH_NODE_NUM returns the number of nodes in a grid.

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
  int node_num;

  if ( strcmp(code,"T3") == 0 )
  {
    node_num = ( nelemx + 1 ) * ( nelemy + 1 );;
  }
  else if ( strcmp(code,"T6") == 0 )
  {
    node_num = ( 2 * nelemx + 1 ) * ( 2 * nelemy + 1 );
  }
  else if ( strcmp(code,"T10") == 0 )
  {
    node_num = ( 3 * nelemx + 1 ) * ( 3 * nelemy + 1 );
  }
  else if ( strcmp(code,"Q4") == 0 )
  {
    node_num = ( nelemx + 1 ) * ( nelemy + 1 );
  }
  else if ( strcmp(code,"Q9") == 0 )
  {
    node_num = ( 2 * nelemx + 1 ) * ( 2 * nelemy + 1 );
  }
  else if ( strcmp(code,"Q16") == 0 )
  {
    node_num = ( 3 * nelemx + 1 ) * ( 3 * nelemy + 1 );
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "MESH_NODE_NUM - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of CODE = \"%s\".\n", code );
    node_num = -1;
    exit ( 1 );
  }

  return node_num;
}
