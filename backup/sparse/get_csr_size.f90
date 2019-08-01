
subroutine get_csr_size(nnz, Bp_size, Ai, Aj, Ax, get_size)
    ! Get csr size
    use sparse_matrix
    implicit none
    integer ( kind = 4 ) nnz, Bp_size, get_size
    integer ( kind = 4 ), dimension(nnz) :: Ai, Aj
    double precision    , dimension(nnz) :: Ax
    integer ( kind = 4 ), dimension(nnz) :: Bj
    double precision    , dimension(nnz) :: Bx
    integer ( kind = 4 ), dimension(Bp_size) :: Bp

    call coo2csr(Ai, Aj, Ax, Bp, Bj, Bx)
    call csr_sort_indices(Bp, Bj, Bx)
    call csr_sum_duplicates(Bp, Bj, Bx)
    get_size = Bp(Bp_size) - 1
    print *, "get csr size done"
    
end subroutine