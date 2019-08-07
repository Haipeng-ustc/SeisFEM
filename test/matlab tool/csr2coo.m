clear;
clc;
% test_lump = load('test_lump.dat');
test_p = load('stif_p.dat');
test_j = load('stif_j.dat');
test_x = load('stif_x.dat');
test_p = test_p + 1;
test_j = test_j + 1;

[np,junk ] = size(test_p);
nrow = np - 1;

for i = nrow: -1 : 1
    k1 = test_p(i+1) - 1;
    k2 = test_p(i);
    for k = k1: -1:  k2
        test_i(k) = i;
    end
end 
test_i = test_i';
test = sparse(test_i, test_j, test_x, nrow, nrow);
imagesc(test);