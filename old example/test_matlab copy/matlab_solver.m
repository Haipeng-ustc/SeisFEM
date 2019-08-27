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




