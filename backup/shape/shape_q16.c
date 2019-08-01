void shape_q16(double r, double s, double t[], double dtdr[], double dtds[] )

/*
*****************************************************************************

 SHAPE_Q16 evaluates shape functions for a 16 node reference quadrilateral.

  Reference Element Q16:

    |
    1 13--14--15--16
    |  |   :   :   |
    |  |   :   :   |
    |  9..10..11..12
    S  |   :   :   |
    |  |   :   :   |
    |  5...6...7...8
    |  |   :   :   |
    |  |   :   :   |  
    0  1---2---3---4
    |
    +--0-----R-----1-->

*/

{
	double ra, rb, rc, rd;
	double sa, sb, sc, sd;
	double dbcd, dacd, dabd, dabc;

  	ra = r - 0.000;
  	rb = r - 1.000 / 3.000;
  	rc = r - 2.000 / 3.000;
  	rd = r - 1.000;

  	sa = s - 0.000;
  	sb = s - 1.000 / 3.000;
  	sc = s - 2.000 / 3.000;
  	sd = s - 1.000;

  	t[0]  =   (  81.000 / 4.000 ) * rb * rc * rd * sb * sc * sd;
  	t[1]  = - ( 243.000 / 4.000 ) * ra * rc * rd * sb * sc * sd;
  	t[2]  =   ( 243.000 / 4.000 ) * ra * rb * rd * sb * sc * sd;
  	t[3]  = - (  81.000 / 4.000 ) * ra * rb * rc * sb * sc * sd;
	
  	t[4]  = - ( 243.000 / 4.000 ) * rb * rc * rd * sa * sc * sd;
  	t[5]  =   ( 729.000 / 4.000 ) * ra * rc * rd * sa * sc * sd;
  	t[6]  = - ( 729.000 / 4.000 ) * ra * rb * rd * sa * sc * sd;
 	t[7]  =   ( 243.000 / 4.000 ) * ra * rb * rc * sa * sc * sd;

  	t[8]  =   ( 243.000 / 4.000 ) * rb * rc * rd * sa * sb * sd;
  	t[9]  = - ( 729.000 / 4.000 ) * ra * rc * rd * sa * sb * sd;
  	t[10] =   ( 729.000 / 4.000 ) * ra * rb * rd * sa * sb * sd;
  	t[11] = - ( 243.000 / 4.000 ) * ra * rb * rc * sa * sb * sd;

  	t[12] = - (  81.000 / 4.000 ) * rb * rc * rd * sa * sb * sc;
  	t[13] =   ( 243.000 / 4.000 ) * ra * rc * rd * sa * sb * sc;
  	t[14] = - ( 243.000 / 4.000 ) * ra * rb * rd * sa * sb * sc;
  	t[15] =   (  81.000 / 4.000 ) * ra * rb * rc * sa * sb * sc;

  	dbcd = 3.000 * r * r -  4.000 * r         + 11.000 / 9.000;
  	dacd = 3.000 * r * r - 10.000 * r / 3.000 +  2.000 / 3.000;
  	dabd = 3.000 * r * r -  8.000 * r / 3.000 +  1.000 / 3.000;
  	dabc = 3.000 * r * r -  2.000 * r         +  2.000 / 9.000;

  	dtdr[0]  =   (  81.000 / 4.000 ) * dbcd * sb * sc * sd;
  	dtdr[1]  = - ( 243.000 / 4.000 ) * dacd * sb * sc * sd;
  	dtdr[2]  =   ( 243.000 / 4.000 ) * dabd * sb * sc * sd;
  	dtdr[3]  = - (  81.000 / 4.000 ) * dabc * sb * sc * sd;
  	dtdr[4]  = - ( 243.000 / 4.000 ) * dbcd * sa * sc * sd;
  	dtdr[5]  =   ( 729.000 / 4.000 ) * dacd * sa * sc * sd;
	dtdr[6]  = - ( 729.000 / 4.000 ) * dabd * sa * sc * sd;
  	dtdr[7]  =   ( 243.000 / 4.000 ) * dabc * sa * sc * sd;
  	dtdr[8]  =   ( 243.000 / 4.000 ) * dbcd * sa * sb * sd;
 	dtdr[9]  = - ( 729.000 / 4.000 ) * dacd * sa * sb * sd;
  	dtdr[10] =   ( 729.000 / 4.000 ) * dabd * sa * sb * sd;
  	dtdr[11] = - ( 243.000 / 4.000 ) * dabc * sa * sb * sd;
  	dtdr[12] = - (  81.000 / 4.000 ) * dbcd * sa * sb * sc;
 	dtdr[13] =   ( 243.000 / 4.000 ) * dacd * sa * sb * sc;
 	dtdr[14] = - ( 243.000 / 4.000 ) * dabd * sa * sb * sc;
  	dtdr[15] =   (  81.000 / 4.000 ) * dabc * sa * sb * sc;

  	dbcd = 3.000 * s * s -  4.000 * s         + 11.000 / 9.000;
  	dacd = 3.000 * s * s - 10.000 * s / 3.000 +  2.000 / 3.000;
  	dabd = 3.000 * s * s -  8.000 * s / 3.000 +  1.000 / 3.000;
  	dabc = 3.000 * s * s -  2.000 * s         +  2.000 / 9.000;

  	dtds[0]  =   (  81.000 / 4.000 ) * rb * rc * rd * dbcd;
  	dtds[1]  = - ( 243.000 / 4.000 ) * ra * rc * rd * dbcd;
  	dtds[2]  =   ( 243.000 / 4.000 ) * ra * rb * rd * dbcd;
  	dtds[3]  = - (  81.000 / 4.000 ) * ra * rb * rc * dbcd;
  	dtds[4]  = - ( 243.000 / 4.000 ) * rb * rc * rd * dacd;
  	dtds[5]  =   ( 729.000 / 4.000 ) * ra * rc * rd * dacd;
  	dtds[6]  = - ( 729.000 / 4.000 ) * ra * rb * rd * dacd;
  	dtds[7]  =   ( 243.000 / 4.000 ) * ra * rb * rc * dacd;
  	dtds[8]  =   ( 243.000 / 4.000 ) * rb * rc * rd * dabd;
  	dtds[9]  = - ( 729.000 / 4.000 ) * ra * rc * rd * dabd;
  	dtds[10] =   ( 729.000 / 4.000 ) * ra * rb * rd * dabd;
  	dtds[11] = - ( 243.000 / 4.000 ) * ra * rb * rc * dabd;
  	dtds[12] = - (  81.000 / 4.000 ) * rb * rc * rd * dabc;
  	dtds[13] =   ( 243.000 / 4.000 ) * ra * rc * rd * dabc;
  	dtds[14] = - ( 243.000 / 4.000 ) * ra * rb * rd * dabc;
  	dtds[15] =   (  81.000 / 4.000 ) * ra * rb * rc * dabc;

}
