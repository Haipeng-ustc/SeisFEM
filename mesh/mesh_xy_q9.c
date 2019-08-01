double mesh_xy_q9 ( int nx, int ny, int node_num, double xmin, double xmax, double ymin, double ymax, double **node_xy)

/******************************************************************************/
/*

 mesh_xy_q9 sets the XY coordinates of the nodes.



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
  int i;
  int j;

for(j=0; j < 2 * ny + 1;j++)
{
    for(i=0; i < 2 * nx + 1;i++)
    {
        node_xy[0][ i + j * (2 * nx + 1)] = i * ( xmax - xmin ) / (2.0 * nx);
        node_xy[1][ i + j * (2 * nx + 1)] = j * ( ymax - ymin ) / (2.0 * ny);
    }
}


}




