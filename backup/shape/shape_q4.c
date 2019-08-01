void shape_q4(double r, double s, double t[], double dtdr[], double dtds[] )

/*
*****************************************************************************

 SHAPE_Q4 evaluates shape functions for a 4 node reference quadrilateral.

  Reference Element Q4:

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
	t[0] = ( 1.000 - r ) * ( 1.000 - s );
  	t[1] =           r   * ( 1.000 - s );
  	t[2] =           r   *           s;
  	t[3] = ( 1.000 - r ) *           s;

  	dtdr[0] = - 1.000 + s;
  	dtdr[1] =   1.000 - s;    
  	dtdr[2] =           s;
  	dtdr[3] =         - s;

  	dtds[0] = - 1.000 + r;
  	dtds[1] =         - r;
  	dtds[2] =           r;
  	dtds[3] =   1.000 - r;

}
