double mesh_xy_t3 ( int nx, int ny, int node_num, double xmin, double xmax, double ymin, double ymax, double **node_xy)

/******************************************************************************/
/*

 grid_t3_xy sets the XY coordinates of the nodes.


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
  int i;
  int j;

for(j=0; j<ny+1;j++)
{
    for(i=0; i<nx+1;i++)
    {
        node_xy[0][i+(j) * (nx + 1)] = ((nx - i) * xmin + (i) * xmax) / (1.0 * nx);
        node_xy[1][i+(j) * (nx + 1)] = ((ny - j) * ymin + (j) * ymax) / (1.0 * ny);
    }
}

}

