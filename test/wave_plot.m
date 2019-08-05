clear;
clc;
energy =  load('energy.dat');
u         =  load('u.dat');
node_xy = load('node_xy.dat');
x = node_xy( : , 1);
y = node_xy( : , 2);
[step, junk] = size(u);
for time = 1 : 1 : 20
    uu=u(time,:);
    scatter(x, y, 50, uu);
    pbaspect([ 1 1 1 ]);
    caxis([-1e-5 1e-5]);
    hold off
    drawnow
    pause(0.2)
end