function h = model_topo_fh( p )

  hmax = 5;
  hmin = 1;
%
%  Estimate the minimum and maximum distances.model_topo_fh
%
  dmax = max ( abs ( p24_fd ( p ) ) );
  dmin = min ( abs ( p24_fd ( p ) ) );
%
%  Assign the density of each point.

  h = ( ( dmax - abs ( p24_fd ( p ) )        ) * hmin   ...
      + (        abs ( p24_fd ( p ) ) - dmin ) * hmax ) ...
      / ( dmax                        - dmin );

  return
end