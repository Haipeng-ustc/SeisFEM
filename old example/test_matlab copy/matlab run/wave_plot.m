clear;
clc;
Energy_u =  load('energy_u_shot_1.dat');
Energy_w =  load('energy_w_shot_1.dat');
u        =  load('wavefield_u_shot_1.dat');
w        =  load('wavefield_w_shot_1.dat');
node_xy  =  load('hole_nodes_200x400.txt');

x = node_xy( : , 1);
y = node_xy( : , 2);
[step, node_num] = size(u);
nx = sqrt(Node_num);
ny = sqrt(Node_num);
figure(1)
caxis_value = 1e-12;
for time = 1 : 5 : step
%     uu=reshape(u(time,:),[nx, ny]);
%     ww=reshape(w(time,:),[nx, ny]);
    uu = U(time,:);
    ww = W(time,:);

    subplot(1,2,1)
    %imagesc(uu);
    scatter(x, y, 50, uu);
    pbaspect([ 1 1 1 ]);
    caxis([-caxis_value caxis_value]);
    subplot(1,2,2)
    %imagesc(ww);
    scatter(x, y, 50, ww);
    pbaspect([ 1 1 1 ]);
    caxis([-caxis_value caxis_value]);
    hold off
    drawnow
    pause(0.2)
end