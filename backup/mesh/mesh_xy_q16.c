double mesh_xy_q16(int nx, int ny, int node_num, double xmin, double xmax, double ymin, double ymax, double **node_xy)

/******************************************************************************/
/*

 mesh_xy_q16 sets the XY coordinates of the nodes.


  Grid:

   43-44-45-46-47-48-49
    |        |        |
    |        |        |
   36 37 38 39 40 41 42
    |        |        |
    |        |        |
   29 30 31 32 33 34 35
    |        |        |
    | 3      | 4      |
   22-23-24-25-26-27-28
    |        |        |
    |        |        |
   15 16 17 18 19 20 21
    |        |        |
    |        |        |
    8  9 10 11 12 13 14
    |        |        |
    | 1      | 2      |
    1--2--3--4--5--6--7

*/
{
    int i;
    int j;

    for (j = 0; j < 3 * ny + 1; j++)
    {
        for (i = 0; i < 3 * nx + 1; i++)
        {
            node_xy[0][i + j * (3 * nx + 1)] = i * (xmax - xmin) / (3.0 * nx);
            node_xy[1][i + j * (3 * nx + 1)] = j * (ymax - ymin) / (3.0 * ny);
        }
    }
}
