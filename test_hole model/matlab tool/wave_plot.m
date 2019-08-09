clear;
clc;
Energy_u =  load('Energy_u.dat');
Energy_w =  load('Energy_w.dat');
u        =  load('u.dat');
w        =  load('w.dat');
node_xy  =  load('hole_nodes_200x400.txt');

x = node_xy( : , 1);
y = node_xy( : , 2);
[step, junk] = size(u);
figure(1)
caxis_value = 0.1;
for time = 1 : 1 : step
    uu=u(time,:);
    ww=w(time,:);
    subplot(1,2,1)
    scatter(x, y, 50, uu);
    pbaspect([ 1 0.5 1 ]);
    caxis([-caxis_value caxis_value]);
    subplot(1,2,2)
    scatter(x, y, 50, ww);
    pbaspect([ 1 0.5 1 ]);
    caxis([-caxis_value caxis_value]);
    hold off
    drawnow
    pause(0.2)
end