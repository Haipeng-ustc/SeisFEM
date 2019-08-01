function h = p24_fh ( p )

%*****************************************************************************80
%
%% P24_FH is a nonuniform mesh density function for problem 24.
%
%  Discussion:
%
%    We want small triangles near the boundaries of the domain, and
%    larger ones in the interior.  To achieve this, we assume that P
%    represents a large sampling of points in the region; we compute
%    the minimum and maximum distances of the points to the boundary,
%    and we assign mesh density values of HMIN to the closest points,
%    HMAX to the furthest ones, and linearly vary H between them. 
%
%    Note that the points inside the region have negative signed distance,
%    and those furthest from the boundary have the most negative value.
%    Thus, we take the absolute value of this distance to get the positive
%    distance we would prefer to work with.  The program does not expect to
%    receive input points which are actually outside the region.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    18 January 2016
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, real P(N,2), one or more points, where the mesh density function
%    is to be evaluated.
%
%    Output, real H(N), the value of the mesh density function H(P).
%

%
%  Set the minimum and maximum densities.
%
  hmax = 5;
  hmin = 1;
%
%  Estimate the minimum and maximum distances.
%
  dmax = max ( abs ( p24_fd ( p ) ) );
  dmin = min ( abs ( p24_fd ( p ) ) );
%
%  Assign the density of each point.
%
%   h = ( ( dmax - abs ( p24_fd ( p ) )        ) * hmin   ...
%       + (        abs ( p24_fd ( p ) ) - dmin ) * hmax ) ...
%       / ( dmax                        - dmin );
%   h = 0.5;
  np = size ( p, 1 );

h = 0.025 + 0.05 * sin ( pi * p(1:np,1) ) .* sin ( pi * p(1:np,2) );
  return
end

