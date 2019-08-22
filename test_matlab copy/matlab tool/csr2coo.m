function matrix = csr2coo(matrix_csr_p, matrix_csr_j, matrix_csr_x )

% convert to 1 based index
matrix_csr_p = matrix_csr_p + 1;
matrix_csr_j = matrix_csr_j + 1;

[np,junk ] = size(matrix_csr_p);
nrow = np - 1;

for i = nrow: -1 : 1
    k1 = matrix_csr_p(i+1) - 1;
    k2 = matrix_csr_p(i);
    for k = k1: -1:  k2
        matrix_csr_i(k) = i;
    end
end

matrix_csr_i = matrix_csr_i';

matrix = sparse(matrix_csr_i, matrix_csr_j, matrix_csr_x, nrow, nrow);

end 

