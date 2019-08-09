function d = model_topo_fd ( p )

  v = load('m.dat');
  d = dpoly ( p, v );
  
%   d1 = dexpr ( p, 'y-cos(x)' );
% 
%   g1 = drectangle ( p, 0.0, 1.0, 0.0, 0.5  );
%   g2 = drectangle ( p, 0.0, 0.5, 0.0, 1.0 );
% 
%   d = dunion ( g1, g2 );

  
  return
end
