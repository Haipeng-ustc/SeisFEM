void shape_t3 ( double r, double s, double t[], double dtdr[], double dtds[] )

/*
*****************************************************************************

 SHAPE_T3 evaluates shape functions for a 3 node reference triangle.

  Reference Element T3:

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
  	t[0] = 1.000 - r - s;
  	t[1] =         r    ;
  	t[2] =             s;

  	dtdr[0] = -1.000;
  	dtdr[1] =  1.000;
  	dtdr[2] =  0.000;

  	dtds[0] = -1.000;
  	dtds[1] =  0.000;
  	dtds[2] =  1.000;
}
