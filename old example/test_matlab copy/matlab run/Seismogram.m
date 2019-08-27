function [Seism_U, Seism_W]= Seismogram( U, W, step, Rec_node)

Seism_U = zeros([step, 1]);
Seism_W = zeros([step, 1]);
Seism_U = U(:, Rec_node);
Seism_W = W(:, Rec_node);
end 