Nx = 8001;
Nz = 4001;
x = linspace(0, 2000, Nx)';
z = linspace(0, 1000, Nz)';
vp = zeros(Nx,Nz);
vs = zeros(Nx,Nz);
for i = 1 : Nx
    for j = 1 : Nz
    x_val = x(i); 
    z_val = z(j); 
    [rho, c11, c13, c33, c44, Vp, Vs] = Elastic_Parameters( z_val, x_val);
    vp(i,j) = Vp;
    vs(i,j) = Vs;
    end 
end 
imagesc(x,z,flipud(vp'));
pbaspect([ 1 0.5 1 ]);
caxis([0 3000]);


%  subplot(1,2,1)
%  scatter(z,x,50,vp);
%  
%  subplot(1,2,2)
%  scatter(z,x,50,vs);
%  pbaspect([ 1 0.5 1 ]);
