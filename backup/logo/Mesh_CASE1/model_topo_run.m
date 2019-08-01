 clear;
 clc;
 fd   =  @model_topo_fd;
 fh   =  @huniform;
 h0   =  5 ;  %  10 
 box  =  [0 , 0 ; 1000 , 2000 ];
 iteration_max = 800; % 500 
 pfix = load('m.dat');
 [ p , t ]=  distmesh_2d( fd, fh, h0, box, iteration_max ,pfix );
 
 post_2d ( p, t, fh );
 [ Node_num_test, junk ] = size ( p );
 [ Element_num_test , junk ] = size ( t );
 p = p';
 t = t';
 node_show = 0;
 triangle_show = 1;

 triangulation_order3_plot ( 'm_mesh.eps', Node_num_test, p, Element_num_test, ...
                            t, node_show, triangle_show );
%
%  Write a text file containing the nodes.
%
  r8mat_write ( 'm_nodes.txt', 2, Node_num_test, p );
%
%  Write a text file containing the triangles.
%
  i4mat_write ( 'm_elements.txt', 3, Element_num_test, t );
  
  triangulation_rcm ( 'm' );
  
  