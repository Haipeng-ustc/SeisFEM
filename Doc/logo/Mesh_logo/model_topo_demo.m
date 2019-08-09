function model_topo_demo ( iteration_max, h, fh )

  if ( nargin < 1 )
    iteration_max = 50;
    fprintf ( 1, '\n' );
    fprintf ( 1, 'P24_DEMO - Note:\n' );
    fprintf ( 1, '  No value of ITERATION_MAX was supplied.\n' );
    fprintf ( 1, '  The default value ITERATION_MAX = %d will be used.\n', ...
      iteration_max );
  end

  if ( nargin < 2 )
    h = 0.015;
    fprintf ( 1, '\n' );
    fprintf ( 1, 'P24_DEMO - Note:\n' );
    fprintf ( 1, '  No value of H was supplied.\n' );
    fprintf ( 1, '  The default value H = %g will be used.\n', h );
  end

  if ( nargin < 3 )
    fh = @p24_fh;
    fprintf ( 1, '\n' );
    fprintf ( 1, 'P24_DEMO - Note:\n' );
    fprintf ( 1, '  No value of FH was supplied.\n' );
    fprintf ( 1, '  The default variable density function will be used.\n', h );
  end
%
%  Put the random number generator into a fixed initial state.
%
  rand ( 'state', 111 );
%
%  Set the rendering method for the current figure to Z-buffering.
%
  set ( gcf, 'rend', 'z' );

  fprintf ( 1, '\n' );
  fprintf ( 1, 'Problem 24:\n' );
  fprintf ( 1, '  The mysterious hand, h = %f\n', h )

  fd = @model_topo_fd;
  box  =  [0 , 0 ; 1000 , 1200 ];

  fixed = load('model_topography.dat');

  [ p, t ] = distmesh_2d ( fd, fh, h, box, iteration_max, fixed );

  post_2d ( p, t, fh )
%
%  Write a PostScript image of the triangulation.
%
  [ node_num, junk ] = size ( p );
  [ tri_num , junk ] = size ( t );
  p = p';
  t = t';
  node_show = 0;
  triangle_show = 1;

  triangulation_order3_plot ( 'model_topo_mesh.eps', node_num, p, tri_num, ...
    t, node_show, triangle_show );
%
%  Write a text file containing the nodes.
%
  r8mat_write ( 'model_topo_nodes.txt', 2, node_num, p );
%
%  Write a text file containing the triangles.
%
  i4mat_write ( 'model_topo_elements.txt', 3, tri_num, t );

  return
end
