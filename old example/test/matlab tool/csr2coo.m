clear;
clc;
% test_lump = load('test_lump.dat');
test_csr_p = load('mass_p.dat');
test_csr_j = load('mass_j.dat');
test_csr_x = load('mass_x.dat');
test_csr_p = test_csr_p + 1;
test_csr_j = test_csr_j + 1;

[np,junk ] = size(test_csr_p);
nrow = np - 1;

for i = nrow: -1 : 1
    k1 = test_csr_p(i+1) - 1;
    k2 = test_csr_p(i);
    for k = k1: -1:  k2
        test_csr_i(k) = i;
    end
end 
test_csr_i = test_csr_i';
test_csr = sparse(test_csr_i, test_csr_j, test_csr_x, nrow, nrow);

tol = 10e-6;
[Lower ,  Upper]  = ilu(test_csr);
rhs = [0:1:nrow-1]';
x = gmres(test_csr, rhs, 10, tol, 100, Lower, Upper);


