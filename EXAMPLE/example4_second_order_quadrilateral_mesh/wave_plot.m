clear;
clc;
u        =  load('wavefield_u_shot_1.txt');
w        =  load('wavefield_w_shot_1.txt');
node_xy  =  load('../mesh/node_xy.txt');

x = node_xy( : , 1);
y = node_xy( : , 2);
[step, node_num] = size(u);
nx = sqrt(node_num);
ny = sqrt(node_num);
figure(1)
caxis_value = 0.1;
for time = 1 : 1 : step
%     uu=reshape(u(time,:),[nx, ny]);
%     ww=reshape(w(time,:),[nx, ny]);
    uu = u(time,:);
    ww = w(time,:);

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