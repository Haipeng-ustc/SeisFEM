function distmesh_test ( )

%*****************************************************************************80
%
%% DISTMESH_TEST tests DISTMESH.
%
%  Licensing:
%
%    (C) 2004 Per-Olof Persson. 
%    See COPYRIGHT.TXT for details.
%
%  Modified:
%
%    09 January 2019
%
%  Reference:
%
%    Per-Olof Persson, Gilbert Strang,
%    A Simple Mesh Generator in MATLAB,
%    SIAM Review,
%    Volume 46, Number 2, June 2004, pages 329-345.
%
%  Local parameters:
%
%    Local, integer ITERATION_MAX, the maximum number of iterations that 
%    DISTMESH should take.  (The program might take fewer iterations if it 
%    detects convergence.)
%
%    Local, real H, the desired initial spacing.  If a uniform mesh density is used,
%    this is the approximate spacing throughout the region.
%
%    Local, pointer FH, an inline formula, or the name of an M file, which calculates
%    the mesh density function.
%
  addpath ( '../distmesh' );

  timestamp ( );
  fprintf ( 1, '\n' );
  fprintf ( 1, 'DISTMESH_TEST\n' );
  fprintf ( 1, '  Test DISTMESH.\n' );
%
%  Problem 1, the circle, with spacings H = 0.4, 0.2, 0.1.
%
  iteration_max = 200;
  h = 0.40;
  p01_demo ( iteration_max, h );
  filename = 'p01_mesh_h040.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );

  iteration_max = 200;
  h = 0.20;
  p01_demo ( iteration_max, h );
  filename = 'p01_mesh_h020.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );

  iteration_max = 200;
  h = 0.10;
  p01_demo ( iteration_max, h );
  filename = 'p01_mesh_h010.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 2, unit circle with a hole.
%
  iteration_max = 200;
  h = 0.10;
  p02_demo ( iteration_max, h )
  filename = 'p02_mesh.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 3, square with a hole, uniform density.
%
  iteration_max = 200;
  h = 0.15;
  fh = @p03_fh;
  p03_demo ( iteration_max, h, fh )
  filename = 'p03_mesh_h015.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 3, square with a hole, finer density near the hole.
%
  iteration_max = 300;
  h = 0.05;
  fh = inline ( 'min(4*sqrt(sum(p.^2,2))-1,2)', 'p' );
  p03_demo ( iteration_max, h, fh )
  filename = 'p03_mesh_h005.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 4, hexagon with a hexagonal hole.
%
  iteration_max = 200;
  h = 0.1;
  p04_demo ( iteration_max, h );
  filename = 'p04_mesh.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 5, the horn.
%
  iteration_max = 200;
  h = 0.020;
  p05_demo ( iteration_max, h );
  filename = 'p05_mesh.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 6, the superellipse.
%  Needs MAPLE to run...
%
  if ( false )
    iteration_max = 200;
    h = 0.08;
    p06_demo ( iteration_max, h );
    filename = 'p06_mesh.png';
    print ( '-dpng', filename );
    fprintf ( 1, '  Graphics saved as "%s"\n', filename );
  else
    fprintf ( 1, '\n' );
    fprintf ( 1, '  Skipping problem 6, the superellipse.\n' );
  end
%
%  Problem 7, the bicycle seat.
%  Needs MAPLE to run...
%
  if ( false )
    iteration_max = 200;
    h = 0.75;
    p07_demo ( iteration_max, h );
    filename = 'p07_mesh.png';
    print ( '-dpng', filename );
    fprintf ( 1, '  Graphics saved as "%s"\n', filename );
  else
    fprintf ( 1, '\n' );
    fprintf ( 1, '  Skipping problem 7, the bicycle seat.\n' );
  end
%
%  Problem 8, the holey pie slice, uniform density.
%
  iteration_max = 200;
  h = 0.025;
  fh = @huniform;
  p08_demo ( iteration_max, h, fh );
  filename = 'p08_mesh_h025.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 8, the holey pie slice, variable density.
%
  iteration_max = 200;
  h = 0.005;
  fh = @p08_fh;
  p08_demo ( iteration_max, h, fh );
  filename = 'p08_mesh_h005.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 9, Jeff Borggaard's square with two hexagonal holes.
%
  iteration_max = 200;
  h = 0.08;
  p09_demo ( iteration_max, h );
  filename = 'p09_mesh.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 10, a simple square.
%
  iteration_max = 200;
  h = 0.10;
  p10_demo ( iteration_max, h )
  filename = 'p10_mesh.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 11, the L-shaped region.
%
  iteration_max = 200;
  h = 0.10;
  p11_demo ( iteration_max, h )
  filename = 'p11_mesh.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 12, John Shadid's H-shaped region.
%
  iteration_max = 200;
  h = 0.05;
  p12_demo ( iteration_max, h )
  filename = 'p12_mesh.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 13, the Sandia Fork.
%
  iteration_max = 200;
  h = 0.025;
  p13_demo ( iteration_max, h )
  filename = 'p13_mesh.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 14, Marcus Garvie's Lake Alpha with Beta island.
%
  if ( true )
    iteration_max = 50;
    h = 20.0;
    fh = @huniform;
    p14_demo ( iteration_max, h, fh )
    filename = 'p14_mesh_h20.png';
    print ( '-dpng', filename );
    fprintf ( 1, '  Graphics saved as "%s"\n', filename );
  else
    fprintf ( 1, '\n' );
    fprintf ( 1, '  Skipping problem 14.\n' );
  end
%
%  Problem 14, Marcus Garvie's Lake Alpha with Beta island.
%  You may get endless warnings from DELAUNAYN about duplicate points.
%  Shut off this annoying drivel with "warning off".
%
  if ( true )
    warning off
    iteration_max = 50;
    h = 10.0;
    fh = @p14_fh;
    p14_demo ( iteration_max, h, fh )
    filename = 'p14_mesh_h10.png';
    print ( '-dpng', filename );
    fprintf ( 1, '  Graphics saved as "%s"\n', filename );
  else
    fprintf ( 1, '\n' );
    fprintf ( 1, '  Skipping problem 14.\n' );
  end
%
%  Problem 15, Sangbum Kim's forward step problem.
%
  iteration_max = 50;
  h = 0.2;
  p15_demo ( iteration_max, h )
  filename = 'p15_mesh.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 16, Kevin Pond's elbow.
%
  iteration_max = 50;
  h = 0.05;
  p16_demo ( iteration_max, h )
  filename = 'p16_mesh.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 17, Reuleaux triangle problem.
%
  iteration_max = 50;
  h = 0.05;
  p17_demo ( iteration_max, h )
  filename = 'p17_mesh.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 18, Dumbbell problem.
%
  iteration_max = 50;
  h = 0.100;
  p18_demo ( iteration_max, h )
  filename = 'p18_mesh.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 19, Dumbbell problem.
%
  iteration_max = 200;
  h = 0.025;
  p19_demo ( iteration_max, h )
  filename = 'p19_mesh.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 20, ICAM Wright House problem.
%
  iteration_max = 2;
  h = 2.0;
  p20_demo ( iteration_max, h )
  filename = 'p20_mesh.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 21, Zhu Wang's quarter round problem.
%
  iteration_max = 200;
  h = 0.100;
  p21_demo ( iteration_max, h )
  filename = 'p21_mesh.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 22, Hans-Werner van Wyk's Big C.
%
  iteration_max = 200;
  h = 0.100;
  p22_demo ( iteration_max, h )
  filename = 'p22_mesh.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 23, Mike Schneier's nonuniform square.
%
  iteration_max = 200;
  h = 0.025;
  p23_demo ( iteration_max, h )
  filename = 'p23_mesh.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Problem 24, the hand.
%
  iteration_max = 50;
  h = 0.015;
  fh = @p24_fh;
  p24_demo ( iteration_max, h, fh )
  filename = 'p24_mesh.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Graphics saved as "%s"\n', filename );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'DISTMESH_TEST\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( );

  rmpath ( '../distmesh' );

  return
end
