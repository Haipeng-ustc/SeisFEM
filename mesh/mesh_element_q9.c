
int *mesh_q9_element ( int nelemx, int nelemy )

/******************************************************************************/
/*
  Purpose: 

    GRID_Q9_ELEMENT produces a grid of 9 node quadrilaterals.

  Example:

    Input:

      NELEMX = 3, NELEMY = 2

    Output:

      ELEMENT_NODE = 
         1,  3, 17, 15,  2, 10, 16,  8,  9;
         3,  5, 19, 17,  4, 12, 18, 10, 11;
         5,  7, 21, 19,  6, 14, 20, 12, 13;
        15, 17, 31, 29, 16, 24, 30, 22, 23;
        17, 19, 33, 31, 18, 26, 32, 24, 25;
        19, 21, 35, 33, 20, 28, 34, 26, 27.

  Grid:

   29---30---31---32---33---34---35
    |    .    |    .    |    .    |
    |    .    |    .    |    .    |
   22 . 23 . 24 . 25 . 26 . 27 . 28
    |    .    |    .    |    .    |
    | 4  .    | 5  .    | 6  .    |
   15---16---17---18---19---20---21
    |    .    |    .    |    .    |
    |    .    |    .    |    .    |
    8 .  9 . 10 . 11 . 12 . 13 . 14
    |    .    |    .    |    .    |
    | 1  .    | 2  .    | 3  .    |
    1----2----3----4----5----6----7

  Element Q9:

    |
    1  4--7--3
    |  |     |
    |  |     |
    S  8  9  6
    |  |     |
    |  |     |
    0  1--5--2
    |
    +--0--R--1-->

*/
{
  int c;
  int e;
  int element;
  int *element_node;
  int element_order = 9;
  int i;
  int j;
  int n;
  int ne;
  int nw;
  int s;
  int se;
  int sw;
  int w;

  element_node = ( int * ) malloc ( element_order*nelemx*nelemy * sizeof ( int ) );
/*
  Node labeling:

    NW----N----NE
     |          |
     W    C     E
     |          |
    SW----S----SE
*/
  element = 0;

  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = 2 * ( j - 1 )  * ( 2 * nelemx + 1 ) + 2 * ( i - 1 ) + 1;
      w  = sw +               2 * nelemx + 1;
      nw = sw +         2 * ( 2 * nelemx + 1 );

      s  = sw + 1;
      c  = sw + 1 +               2 * nelemx + 1;
      n  = sw + 1 +         2 * ( 2 * nelemx + 1 );

      se = sw + 2;
      e  = sw + 2 +               2 * nelemx + 1;
      ne = sw + 2 +         2 * ( 2 * nelemx + 1 );

      element_node[0 + element * element_order] = sw;
      element_node[1 + element * element_order] = se;
      element_node[2 + element * element_order] = ne;
      element_node[3 + element * element_order] = nw;
      element_node[4 + element * element_order] = s;
      element_node[5 + element * element_order] = e;
      element_node[6 + element * element_order] = n;
      element_node[7 + element * element_order] = w;
      element_node[8 + element * element_order] = c;

      element = element + 1;
    }
  }

  return element_node;
}
