% Wave Equation Evolution
function [U, W, Energy] = Wave_Evolution_CFD(M, Cx, Cz, Cxx, Cxz, Czz, Kx1, Kx2, Kx3, Kx4, Kx5, Kx6, Kx7, ...
Kz1, Kz2, Kz3, Kz4, Kz5, Kz6, Kz7, Node_num, dt, step, f0, t0, Angle_Force, Source_node);
                
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
Source   = zeros(step+2,1);
Source_x   = zeros(Node_num, 1);
Source_z   = zeros(Node_num, 1);

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
for i= 1 : step + 2; 
    time = (i - 1) * dt;
    Source(i) = Seismic_Source( f0, t0, time);
end
Source = diff(Source) ;
Source = Source/max(Source);

for it = 3: step
    time = it * dt;
    fprintf('Time:  %f s at Step: %d \n',time,it);
    
    Source_x(Source_node) = Source(it) * sin( Angle_Force * pi / 180.0 );
    Source_z(Source_node) = Source(it) * cos( Angle_Force * pi / 180.0 );
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
end 
