function [Seism_n_u, Seism_n_w, Rec_location] = multi_receiver(x, z, Rec_num, Rec_nx, Rec_nz, Node_num, ave_size, step, U, W)

Seism_n_u = zeros(step, Rec_num);
Seism_n_w = zeros(step, Rec_num);
Rec_location = zeros(Rec_num,2);
for i = 1 : Rec_num
    
    [Rec_node] = Receiver_Location( x, z, Rec_nx(i), Rec_nz(i), Node_num, ave_size);
    Rec_location(i,:) = [x(Rec_node), z(Rec_node)];    
    [Seism_n_u(:,i), Seism_n_w(:,i)]= Seismogram( U, W, step, Rec_node);
    
end

end 


