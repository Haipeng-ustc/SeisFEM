node_xy = load('node_xy.txt');
element_node = load('element_node.txt');
mass_matrix = load('mass_matrix.txt');
stiffness_matrix= load('stiffness_matrix.txt');
[node_num, junk] =  size(node_xy);
mass = sparse( mass_matrix( : , 1 ), mass_matrix( : , 2 ),  mass_matrix( : , 3 ), node_num, node_num  );
stiffness = sparse( stiffness_matrix( : , 1 ), stiffness_matrix( : , 2 ),  stiffness_matrix( : , 3 ), node_num, node_num  );
acoustic_model = load('acoustic_model.txt');
elastic_model = load('elastic_model.txt');
rho   = elastic_model( : ,1);
c11  = elastic_model( : ,2);
c13  = elastic_model( : ,3);
c33  = elastic_model( : ,4);
c44  = elastic_model( : ,5);
x = node_xy(:,1);
y  =node_xy(:,2);
pml = load('pml.txt');
pml_dx = pml(:,1);
pml_dy = pml(:,2);
pml_dxx = pml(:,3);
pml_dyy = pml(:,4);
figure(1)
subplot(2,2,1)
scatter(x,y,50, pml_dx);
subplot(2,2,2)
scatter(x,y,50, pml_dy);
subplot(2,2,3)
scatter(x,y,50, pml_dxx);
subplot(2,2,4)
scatter(x,y,50, pml_dyy);

mpml = load('mpml.txt');
mpml_dx = mpml(:,1);
mpml_dy = mpml(:,2);
mpml_dxx = mpml(:,3);
mpml_dyy = mpml(:,4);
mpml_dxx_pyx = mpml(:,5);
mpml_dyy_pxy = mpml(:,6);
figure(2)
subplot(3,2,1)
scatter(x,y,50, mpml_dx);
subplot(3,2,2)
scatter(x,y,50, mpml_dy);
subplot(3,2,3)
scatter(x,y,50, mpml_dxx);
subplot(3,2,4)
scatter(x,y,50, mpml_dyy);
subplot(3,2,5)
scatter(x,y,50, mpml_dxx_pyx);
subplot(3,2,6)
scatter(x,y,50, mpml_dyy_pxy);
 mpml_dx = reshape(mpml_dx, [101,101]);
 mpml_dy = reshape(mpml_dy, [101,101]);
 mpml_dxx = reshape(mpml_dxx, [101,101]);
 mpml_dyy = reshape(mpml_dyy, [101,101]);
 mpml_dxx_pyx = reshape(mpml_dxx_pyx, [101,101]);
 mpml_dyy_ = reshape(mpml_dyy_pxy, [101,101]);

mass_lump = load('mass_matrix_lump.txt');

 y = load('y.txt');
 csr = load('csr.txt');
 Bp = load('Bp.txt');
 
% subplot(1,2,1)
% imagesc(mass);
% subplot(1,2,2)
% imagesc(stiffness);


% showmesh
% figure(3)
% showmesh(node_xy,element_node);
% findnode(node_xy);
