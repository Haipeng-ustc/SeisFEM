clear;
clc;
% mass_lump = load('test_lump.dat');
mass_p = load('stif_p.dat');
mass_j = load('stif_j.dat');
mass_x = load('stif_x.dat');
mass_p = mass_p + 1;
mass_j = mass_j + 1;

[np,junk ] = size(mass_p);
nrow = np - 1;

for i = nrow: -1 : 1
    k1 = mass_p(i+1) - 1;
    k2 = mass_p(i);
    for k = k1: -1:  k2
        mass_i(k) = i;
    end
end 
mass_i = mass_i';
mass = sparse(mass_i, mass_j, mass_x, nrow, nrow);
imagesc(mass);