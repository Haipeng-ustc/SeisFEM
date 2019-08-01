
int *mesh_q4_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    MESH_Q4_ELEMENT produces a grid of 4 node quadrilaterals.

    For each element, the nodes are listed in counter-clockwise order.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1, 2,  6,  5;
         2, 3,  7,  6;
         3, 4,  8,  7;
         5, 6, 10,  9;
         6, 7, 11, 10;
         7, 8, 12, 11.

  Grid:

    9---10---11---12
    |    |    |    |
    |    |    |    |
    |  4 |  5 |  6 |
    |    |    |    |
    5----6----7----8
    |    |    |    |
    |    |    |    |
    |  1 |  2 |  3 |
    |    |    |    |
    1----2----3----4

  Element Q4:

    |
    1  4-----3
    |  |     |
    |  |     |
    S  |     |
    |  |     |
    |  |     |
    0  1-----2
    |
    +--0--R--1-->

*/
{
  int element;
  int *element_node;
  int element_order = 4;
  int i;
  int j;
  int ne;
  int nw;
  int se;
  int sw;

  element_node = ( int * ) malloc ( element_order*nelemx*nelemy * sizeof ( int ) );
/*
  Node labeling:

    NW---NE
     |    |
    SW---SE
*/
  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = i     + ( j - 1 ) * ( nelemx + 1 );
      se = i + 1 + ( j - 1 ) * ( nelemx + 1 );
      nw = i     +   j       * ( nelemx + 1 );
      ne = i + 1 +   j       * ( nelemx + 1 );

      element_node[0 + element * element_order] = sw;
      element_node[1 + element * element_order] = se;
      element_node[2 + element * element_order] = ne;
      element_node[3 + element * element_order] = nw;

      element = element + 1;
    }
  }

  return element_node;
}
