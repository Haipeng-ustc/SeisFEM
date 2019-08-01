
int *mesh_t3_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    MESH_T3_ELEMENT produces a grid of pairs of 3 node triangles.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1,  2,  5;
         6,  5,  2;
         2,  3,  6;
         7,  6,  3;
         3,  4,  7;
         8,  7,  4;
         5,  6,  9;
        10,  9,  6;
         6,  7, 10;
        11, 10,  7;
         7,  8, 11;
        12, 11,  8.

  Grid:

    9---10---11---12
    |\ 8 |\10 |\12 |
    | \  | \  | \  |
    |  \ |  \ |  \ |
    |  7\|  9\| 11\|
    5----6----7----8
    |\ 2 |\ 4 |\ 6 |
    | \  | \  | \  |
    |  \ |  \ |  \ |
    |  1\|  3\|  5\|
    1----2----3----4

  Element T3:

    |
    1  3
    |  |\
    |  | \
    S  |  \
    |  |   \
    |  |    \
    0  1-----2
    |
    +--0--R--1-->

*/
{
  int element;
  int *element_node;
  int element_order = 3;
  int i;
  int j;
  int ne;
  int nw;
  int se;
  int sw;

  element_node = ( int * ) malloc ( element_order*2*nelemx*nelemy * sizeof ( int ) );
/*
  Node labeling:

    NW--NE
     |\ |
     | \|
    SW--SE
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

      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = nw;
      element = element + 1;

      element_node[0+element*element_order] = ne;
      element_node[1+element*element_order] = nw;
      element_node[2+element*element_order] = se;
      element = element + 1;
    }
  }

  return element_node;
}
