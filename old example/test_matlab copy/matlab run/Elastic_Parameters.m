function [rho, c11, c13, c33, c44, Vp, Vs] = Elastic_Parameters(node_num, node, node_xy);

% Case 1 isotropic media by Komatisch and Tromp (2003)

Vp   = 2000;
Vs   = 1154.7;

% Vp   = 3000;
% Vs   = 1000;
rho  = 2200;
miu  = Vs^2*rho;
lamda= Vp^2*rho-2*miu;
c11  = lamda+2*miu;
c13  = lamda;
c33  = lamda+2*miu;
c44  = miu;

end