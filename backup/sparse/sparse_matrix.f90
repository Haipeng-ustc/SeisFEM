module sparse_matrix
contains 

! subroutine for changing the matrices from coo format to csr format

subroutine csr_matvec(Ap_size, Ax_size, Ap, Aj, Ax, x_size, x, y)
    ! Compute y = A*x for CSR matrix A and dense vectors x, y
    ! input : 
    ! Ap_size: the size of Ap in csr
    ! Ax_size: the size of Aj and Ax in csr
    ! Ap, Aj, Ax: row offsets, column indices and values
    ! x_size : the vector size of x and y
    ! x : x vector
    !
    ! output:
    ! y:  y = A*x
    !
    integer( kind = 4 ) Ap_size, Ax_size, x_size, i
    integer( kind = 4 ), dimension(Ap_size) :: Ap
    integer( kind = 4 ), dimension(Ax_size) :: Aj
    double precision   , dimension(Ax_size) :: Ax
    double precision   , dimension(x_size ) :: x, y
    
    !$omp parallel default(none) shared(Ap, Aj, Ax, x, y) private(i)
    !$omp do
    do i = 1, size(Ap)-1
        y(i) = dot_product(Ax(Ap(i):Ap(i+1)-1), x(Aj(Ap(i):Ap(i+1)-1)))
    end do
    !$omp end do
    !$omp end parallel
    
end subroutine csr_matvec

subroutine coo2csr_canonical(nnz, Bp_size, csr_size, Ai, Aj, Ax, Bp, Bj, Bx)
    ! Get csr size and then 
    ! Converts from COO (Ai, Aj, Ax) into CSR (Bp, Bj, Bx)
    ! Row and column indices are *not* assumed to be ordered.
    ! Duplicate entries are summed up and the indices are ordered.
    implicit none
    integer ( kind = 4 ) nnz, Bp_size, csr_size
    integer ( kind = 4 ), dimension(nnz) :: Ai, Aj
    double precision    , dimension(nnz) :: Ax
    integer ( kind = 4 ), dimension(nnz) :: Bj
    double precision    , dimension(nnz) :: Bx
    integer ( kind = 4 ), dimension(Bp_size) :: Bp
    print *

    write(*,*) "coo2csr_canonical: Initialization."
    Ai = Ai + 1
    Aj = Aj + 1
    write(*,*) "coo2csr_canonical: convert Ai, Aj to 1-based index for Fortran use"
    call coo2csr(nnz, Bp_size, Ai, Aj, Ax, Bp, Bj, Bx)
    call csr_sort_indices(nnz, Bp_size, Bp, Bj, Bx)
    call csr_sum_duplicates(nnz, Bp_size, Bp, Bj, Bx)
    csr_size = Bp(Bp_size) - 1
    write(*,*) "coo2csr_canonical: Get csr size done, Bp and Bj size are:",Bp_size,csr_size
    write(*,*) "coo2csr_canonical: coo2csr -> csr_sort_indices -> csr_sum_duplicates"
    Bp = Bp - 1
    Bj = Bj - 1
    write(*,*) "coo2csr_canonical: convert Ai, Aj to 0-based index for C use"

    write(*,*) "coo2csr_canonical: Normal End"


    
end subroutine

!subroutine coo2csr_canonical(nnz, Bp_size, Bj_size, Ai, Aj, Ax, Bp, Bj, Bx)
!    ! Converts from COO (Ai, Aj, Ax) into CSR (Bp, Bj, Bx)
!    ! Row and column indices are *not* assumed to be ordered.
!    ! Duplicate entries are summed up and the indices are ordered.
!    implicit none
!    integer  :: nnz, Bp_size, Bj_size
!    integer ( kind = 4 ), dimension(  nnz  ) :: Ai, Aj
!    double precision    , dimension(  nnz  ) :: Ax, Bx_
!    integer ( kind = 4 ), dimension(  nnz  ) :: Bj_
!    integer ( kind = 4 ), dimension(Bp_size) :: Bp
!    integer ( kind = 4 ), dimension(Bj_size) :: Bj
!    double precision    , dimension(Bj_size) :: Bx
!    call coo2csr(nnz, Bp_size, Ai, Aj, Ax, Bp, Bj_, Bx_)
!    call csr_sort_indices(nnz, Bp_size, Bp, Bj_, Bx_)
!    call csr_sum_duplicates(nnz, Bp_size, Bp, Bj_, Bx_)
!    Bj  = Bj_(:Bj_size)
!    Bx  = Bx_(:Bj_size)
!    print *, "coo2csr -> csr_sort_indices -> csr_sum_duplicates -> csr done"
!end subroutine

subroutine coo2csr(nnz, Bp_size, Ai, Aj, Ax, Bp, Bj, Bx)
        ! Converts from COO (Ai, Aj, Ax) into CSR (Bp, Bj, Bx)
        ! Row and column indices are *not* assumed to be ordered.
        ! Duplicate entries are carried over to the CSR representation.
	implicit none 
	integer ( kind = 4 ) nnz, Bp_size
        integer ( kind = 4 ), dimension(  nnz  ) :: Ai, Aj
        double precision    , dimension(  nnz  ) :: Ax
        integer ( kind = 4 ), dimension(Bp_size) :: Bp
	integer ( kind = 4 ), dimension(  nnz  ) :: Bj
        double precision    , dimension(  nnz  ) :: Bx
	integer :: n, i, n_row, cumsum, temp, row, dest
        n_row = Bp_size - 1
        Bp = 0

       do n = 1,nnz
	  Bp(Ai(n)) = Bp(Ai(n)) + 1
       end do

        cumsum = 1
        do i = 1, n_row
            temp = Bp(i)
            Bp(i) = cumsum
            cumsum = cumsum + temp
        end do

        do n = 1, nnz
            row = Ai(n)
            dest = Bp(row)
            Bj(dest) = Aj(n)
            Bx(dest) = Ax(n)
            Bp(row) = Bp(row) + 1
        end do
        Bp(2:) = Bp(:n_row)
        Bp(1) = 1
end subroutine



subroutine csr_sum_duplicates(nnz, Bp_size, Bp, Bj, Bx)
        ! Sum together duplicate column entries in each row of CSR matrix A
        ! The column indicies within each row must be in sorted order.
        ! Explicit zeros are retained.
        ! Ap, Aj, and Ax will be modified *inplace*
   	integer ( kind = 4 ) nnz, Bp_size
        integer ( kind = 4 ), dimension(Bp_size) :: Bp
	integer ( kind = 4 ), dimension(  nnz  ) :: Bj
        double precision    , dimension(  nnz  ) :: Bx
        integer :: r1, r2, i, j, jj
        double precision :: x
        nnz = 1
        r2 = 1
        do i = 1, Bp_size - 1
            r1 = r2
            r2 = Bp(i+1)
            jj = r1
            do while (jj < r2)
                j = Bj(jj)
                x = Bx(jj)
                jj = jj + 1
                do while (jj < r2)
                    if (Bj(jj) == j) then
                        x = x + Bx(jj)
                        jj = jj + 1
                    else
                        exit
                    end if
                end do
                Bj(nnz) = j
                Bx(nnz) = x
                nnz = nnz + 1
            end do
            Bp(i+1) = nnz
        end do
end subroutine
        
subroutine csr_sort_indices(nnz, Bp_size, Bp, Bj, Bx)
        ! Sort CSR column indices inplace

	integer ( kind = 4 ) nnz, Bp_size
        integer ( kind = 4 ), dimension(Bp_size) :: Bp
	integer ( kind = 4 ), dimension(  nnz  ) :: Bj
        double precision    , dimension(  nnz  ) :: Bx
        integer :: i, r1, r2, l, idx(nnz)
        do i = 1, Bp_size-1
            r1 = Bp(i)
            r2 = Bp(i+1)-1
            l = r2-r1+1
            idx(:l) = iargsort_quicksort(Bj(r1:r2))
            Bj(r1:r2) = Bj(r1+idx(:l)-1)
            Bx(r1:r2) = Bx(r1+idx(:l)-1)
        end do

end subroutine

pure function iargsort_quicksort(vec_) result(map)
        integer, intent(in) :: vec_(:)
        integer :: map(size(vec_))
        integer, parameter :: levels = 300
        integer, parameter :: max_interchange_sort_size = 20
        integer :: i,left,right,l_bound(levels),u_bound(levels)
        integer :: pivot
        integer :: vec(size(vec_))

        vec = vec_

        forall(i=1:size(vec)) map(i) = i

        l_bound(1) = 1
        u_bound(1) = size(vec)
        i = 1
        do while(i >= 1)
            left = l_bound(i)
            right = u_bound(i)
            if (right - left < max_interchange_sort_size) then
                if (left < right) call interchange_sort_map_int(vec(left:right),map(left:right))
                i = i - 1
            else
                pivot = (vec(left) + vec(right)) / 2
                left = left - 1
                right = right + 1
                do
                do
                    left = left + 1
                    if (vec(left) >= pivot) exit
                end do
                do
                    right = right - 1
                    if (vec(right) <= pivot) exit
                end do
                if (left < right) then
                    call swap_int(vec(left),vec(right))
                    call swap_int(map(left),map(right))
                elseif(left == right) then
                    if (left == l_bound(i)) then
                        left = left + 1
                    else
                        right = right - 1
                    end if
                    exit
                else
                    exit
                end if
                end do
                u_bound(i + 1) = u_bound(i)
                l_bound(i + 1) = left
                u_bound(i) = right
                i = i + 1
            end if
        end do
end function

pure elemental subroutine swap_int(x,y)
            integer, intent(in out) :: x,y
            integer :: z
            z = x
            x = y
            y = z
end subroutine

pure subroutine interchange_sort_map_int(vec,map)
            integer, intent(in out) :: vec(:)
            integer, intent(in out) :: map(:)
            integer :: i,j
            do i = 1,size(vec) - 1
                j = minloc(vec(i:),1)
                if (j > 1) then
                    call swap_int(vec(i),vec(i + j - 1))
                    call swap_int(map(i),map(i + j - 1))
                end if
            end do
end subroutine

end module