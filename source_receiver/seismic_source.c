double seismic_source(double f0, double t0, double magnitide, double t)
/******************************************************************************/
/*
  Purpose:

    seismic_source returns the point source.
    seismic point_source Ricker wavelet  
    f0: dominant frequency 
    t0: onset time
    magnitide : maxmium amplititude
    t : time
    point_source: the point_source value at time t

*/

{
  double point_source;
  double a;
  double pi = 3.141592653589793238462643;
  a = (pi * f0 * (t - t0)) * (pi * f0 * (t - t0));
  point_source = magnitide * (1 - 2 * a) * exp(-a);

  return point_source;
}
