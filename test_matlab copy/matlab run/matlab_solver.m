clear;clc;
% csr format matrices
mass_csr_p   = load('mass_p.dat');
stif1_csr_p = load('stif1_p.dat');
stif2_csr_p = load('stif2_p.dat');
stif3_csr_p = load('stif3_p.dat');
stif4_csr_p = load('stif4_p.dat');
stif5_csr_p = load('stif5_p.dat');
stif6_csr_p = load('stif6_p.dat');
mass_csr_j_x   = load('mass_j_x.dat');
stif1_csr_j_x = load('stif1_j_x.dat');
stif2_csr_j_x = load('stif2_j_x.dat');
stif3_csr_j_x = load('stif3_j_x.dat');
stif4_csr_j_x = load('stif4_j_x.dat');
stif5_csr_j_x = load('stif5_j_x.dat');
stif6_csr_j_x = load('stif6_j_x.dat');
mass  = csr2coo(mass_csr_p,   mass_csr_j_x(:,1), mass_csr_j_x(:,2) );  %phi*phi
stif1 = csr2coo(stif1_csr_p, stif1_csr_j_x(:,1), stif1_csr_j_x(:,2) ); %dphi_dx * dphi_dx
stif2 = csr2coo(stif2_csr_p, stif2_csr_j_x(:,1), stif2_csr_j_x(:,2) ); %dphi_dy * dphi_dy
stif3 = csr2coo(stif3_csr_p, stif3_csr_j_x(:,1), stif3_csr_j_x(:,2) ); %dphi_dx * dphi_dy
stif4 = csr2coo(stif4_csr_p, stif4_csr_j_x(:,1), stif4_csr_j_x(:,2) ); %dphi_dy * dphi_dx
stif5 = csr2coo(stif5_csr_p, stif5_csr_j_x(:,1), stif5_csr_j_x(:,2) ); %phi     * dphi_dx
stif6 = csr2coo(stif6_csr_p, stif6_csr_j_x(:,1), stif6_csr_j_x(:,2) ); %phi     * dphi_dy
clearvars mass_csr_p;
clearvars stif1_csr_p;
clearvars stif2_csr_p;
clearvars stif3_csr_p;
clearvars stif4_csr_p;
clearvars stif5_csr_p;
clearvars stif6_csr_p;
clearvars mass_csr_j_x;
clearvars stif1_csr_j_x;
clearvars stif2_csr_j_x;
clearvars stif3_csr_j_x;
clearvars stif4_csr_j_x;
clearvars stif5_csr_j_x;
clearvars stif6_csr_j_x;
% model_par rho vp vs c11 c13 c33 c44
model_par = load('model_par.dat');
rho = model_par(:,1);
vp  = model_par(:,2);
vs  = model_par(:,3);
c11 = model_par(:,4);
c13 = model_par(:,5);
c33 = model_par(:,6);
c44 = model_par(:,7);
clearvars model_par;
% node_xy and element_node
node_xy = load('hole_nodes_200x400.txt');
element_node = load('hole_elements_200x400.txt');
node_x = node_xy(:,1);
node_z = node_xy(:,2);
[Node_num,junk] = size(node_xy);
[element_num,element_order] = size(element_node);

% parameters for time evolution
dt = 0.0005;
f0 = 15.0;
t0 = 1.2/f0;
step = 2;
h    = 5;
source_x = max(node_x) / 2;
source_z = max(node_z) / 2;

Angle_Force = 90;
% mpml
[mpml_dx,mpml_dz, mpml_dxx, mpml_dzz, mpml_dxx_pzx, mpml_dzz_pxz] = MPML(max(vp), Node_num ,node_x, node_z, h);
% source location
[Source_node] = Source_Location( node_x, node_z, source_x, source_z, Node_num, h);

M   = mass .* rho;
Cx  = M    .* mpml_dx;  %  rho * phi * phi * mpml_dx
Cz  = M    .* mpml_dz;  %  rho * phi * phi * mpml_dz
Cxx = Cx   .* mpml_dx;  %  rho * phi * phi * mpml_dx * mpml_dx
Cxz = Cx   .* mpml_dz;  %  rho * phi * phi * mpml_dx * mpml_dz
Czz = Cx   .* mpml_dz;  %  rho * phi * phi * mpml_dz * mpml_dz
Kx1 = stif1 .* c11;
Kx2 = stif3 .* c13 + c44 .* stif4;
Kx3 = stif2 .* c44;
Kx4 = stif5 .* c11 .* mpml_dxx;
Kx5 = stif5 .* c44 .* mpml_dzz_pxz;
Kx6 = stif6 .* c13 .* mpml_dxx_pzx;
Kx7 = stif6 .* c44 .* mpml_dzz;

Kz1 = stif1 .* c44;
Kz2 = stif3 .* c44 + c13 .* stif4;
Kz3 = stif2 .* c33;
Kz4 = stif5 .* c44 .* mpml_dxx;
Kz5 = stif5 .* c13 .* mpml_dzz_pxz;
Kz6 = stif6 .* c44 .* mpml_dxx_pzx;
Kz7 = stif6 .* c33 .* mpml_dzz;


U1_now   = zeros(Node_num, 1); U2_now   = zeros(Node_num, 1); U3_now   = zeros(Node_num, 1);
U1_old   = zeros(Node_num, 1); U2_old   = zeros(Node_num, 1); U3_old   = zeros(Node_num, 1);
W1_now   = zeros(Node_num, 1); W2_now   = zeros(Node_num, 1); W3_now   = zeros(Node_num, 1);
W1_old   = zeros(Node_num, 1); W2_old   = zeros(Node_num, 1); W3_old   = zeros(Node_num, 1);
Lx1_now  = zeros(Node_num, 1); Lx2_now  = zeros(Node_num, 1); 
Lx3_now  = zeros(Node_num, 1); Lx4_now  = zeros(Node_num, 1);
Lz1_now  = zeros(Node_num, 1); Lz2_now  = zeros(Node_num, 1);
Lz3_now  = zeros(Node_num, 1); Lz4_now  = zeros(Node_num, 1);
U_now    = zeros(Node_num, 1);       
W_now    = zeros(Node_num, 1);
Source   = zeros(Node_num, 1);
U        = zeros(step, Node_num); 
W        = zeros(step, Node_num);
Energy   = zeros(step, 1);
a0 = 1 / dt^2;
a1 = 1 / (2*dt);
a2 = 2 * a0;
A1  =  a0 * M +  a1 * 2 * Cx;
A2  =  a0 * M +  a1 * (Cx + Cz);
A3  =  a0 * M +  a1 * 2 *  Cz;
AL  =  2  * a1 * M;
[Lower1 ,  Upper1]  = ilu(A1);
[Lower2 ,  Upper2]  = ilu(A2);
[Lower3 ,  Upper3]  = ilu(A3);
[LowerL ,  UpperL]  = ilu(AL);

U(1,:) = U1_old + U2_old + U3_old;
U(2,:) = U1_now + U2_now + U3_now;
W(1,:) = W1_old + W2_old + W3_old;
W(2,:) = W1_now + W2_now + W3_now;
Energy(1) = 0.0;
Energy(2) = 0.0;
tol = 10e-6;
for it = 3: 500
    time = it * dt;
    fprintf('Time:  %f s at Step: %d \n',time,it);
    Source(Source_node) = Seismic_Source( f0, t0, time);
    Source_x = Source * sin( Angle_Force * pi / 180.0 );
    Source_z = Source * cos( Angle_Force * pi / 180.0 );
 % CFD Scheme
    bx1 = Source_x + M * Lx1_now + M * Lx2_now -  Kx1 * U_now - (Cxx - a2 * M) * U1_now - (a0 * M - a1 * 2 * Cx ) * U1_old;
    bx2 =                                      -  Kx2 * W_now - (Cxz - a2 * M) * U2_now - (a0 * M - a1 * (Cx+Cz)) * U2_old;
    bx3 =            M * Lx3_now + M * Lx4_now -  Kx3 * U_now - (Czz - a2 * M) * U3_now - (a0 * M - a1 * 2 * Cz ) * U3_old;
    bx4 =   2 * a1 * M * Lx1_now               -  Kx4 * U_now -  Cx * Lx1_now;
    bx5 =   2 * a1 * M * Lx2_now               -  Kx5 * W_now -  Cz * Lx2_now;
    bx6 =   2 * a1 * M * Lx3_now               -  Kx6 * W_now -  Cx * Lx3_now;
    bx7 =   2 * a1 * M * Lx4_now               -  Kx7 * U_now -  Cz * Lx4_now;

    bz1 = Source_z + M * Lz1_now + M * Lz2_now -  Kz1 * W_now - (Cxx - a2 * M) * W1_now - (a0 * M - a1 * 2 * Cx ) * W1_old;
    bz2 =                                      -  Kz2 * U_now - (Cxz - a2 * M) * W2_now - (a0 * M - a1 * (Cx+Cz)) * W2_old;
    bz3 =            M * Lz3_now + M * Lz4_now -  Kz3 * W_now - (Czz - a2 * M) * W3_now - (a0 * M - a1 * 2 * Cz ) * W3_old;
    bz4 =   2 * a1 * M * Lz1_now               -  Kz4 * W_now -  Cx * Lz1_now;
    bz5 =   2 * a1 * M * Lz2_now               -  Kz5 * U_now -  Cz * Lz2_now;
    bz6 =   2 * a1 * M * Lz3_now               -  Kz6 * U_now -  Cx * Lz3_now;
    bz7 =   2 * a1 * M * Lz4_now               -  Kz7 * W_now -  Cz * Lz4_now;

    U1_old = U1_now;
    U2_old = U2_now;
    U3_old = U3_now;
    W1_old = W1_now;
    W2_old = W2_now;
    W3_old = W3_now;
       
    U1_now  = gmres(A1, bx1, 10, tol, 100, Lower1, Upper1);
    U2_now  = gmres(A2, bx2, 10, tol, 100, Lower2, Upper2);
    U3_now  = gmres(A3, bx3, 10, tol, 100, Lower3, Upper3);
    Lx1_now = gmres(AL, bx4, 10, tol, 100, LowerL, UpperL);
    Lx2_now = gmres(AL, bx5, 10, tol, 100, LowerL, UpperL);
    Lx3_now = gmres(AL, bx6, 10, tol, 100, LowerL, UpperL);
    Lx4_now = gmres(AL, bx7, 10, tol, 100, LowerL, UpperL);
    
    W1_now  = gmres(A1, bz1, 10, tol, 100, Lower1, Upper1);
    W2_now  = gmres(A2, bz2, 10, tol, 100, Lower2, Upper2);
    W3_now  = gmres(A3, bz3, 10, tol, 100, Lower3, Upper3);
    Lz1_now = gmres(AL, bz4, 10, tol, 100, LowerL, UpperL);
    Lz2_now = gmres(AL, bz5, 10, tol, 100, LowerL, UpperL);
    Lz3_now = gmres(AL, bz6, 10, tol, 100, LowerL, UpperL);
    Lz4_now = gmres(AL, bz7, 10, tol, 100, LowerL, UpperL);

    U_now  = U1_now + U2_now + U3_now;
    W_now  = W1_now + W2_now + W3_now;    
    U(it,:) = U_now;
    W(it,:) = W_now;
    
    % Energy Caculation
    Energy(it) = 0.5* U_now' * U_now + 0.5* W_now' * W_now;
    if (Energy(it) >= 10e8)
        fprintf('Warning: Time Evolution is Unstable, energy is more than 10e8 !!! \n');
        pause;
    end
   
end







