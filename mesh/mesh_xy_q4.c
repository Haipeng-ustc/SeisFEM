double mesh_xy_q4(int nx, int ny, int node_num, double xmin, double xmax, double ymin, double ymax, double **node_xy)

/******************************************************************************/
/*

 mesh_xy_q4 sets the XY coordinates of the nodes.

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
  int i;
  int j;

  for (j = 0; j < ny + 1; j++)
  {
    for (i = 0; i < nx + 1; i++)
    {
      node_xy[0][i + (j) * (nx + 1)] = ((nx - i) * xmin + (i)*xmax) / (1.0 * nx);
      node_xy[1][i + (j) * (nx + 1)] = ((ny - j) * ymin + (j)*ymax) / (1.0 * ny);
    }
  }
}
