module coo2csr_lib
contains 

! subroutine for changing the matrices from coo format to csr format

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
    integer ( kind = 4 ), dimension(nnz) :: Bj_temp
    double precision    , dimension(nnz) :: Bx_temp
    integer ( kind = 4 ), dimension(Bp_size) :: Bp_temp
    print *

    write(*,*) "coo2csr_canonical: Initialization."
    Bp = 0
    Bj = 0
    Bx = 0.0
    Ai = Ai + 1
    Aj = Aj + 1
    call coo2csr(nnz, Bp_size, Ai, Aj, Ax, Bp_temp, Bj_temp, Bx_temp)
    call csr_sort_indices(nnz, Bp_size, Bp_temp, Bj_temp, Bx_temp)
    call csr_sum_duplicates(nnz, Bp_size, Bp_temp, Bj_temp, Bx_temp)
    call csr_eliminate_zeros(nnz, Bp_size, Bp_temp, Bj_temp, Bx_temp, Bp, Bj, Bx)
    csr_size = Bp(Bp_size) - 1
    write(*,*) "coo2csr_canonical: csr_p_size, csr_j_size:", Bp_size, csr_size
    write(*,*) "coo2csr_canonical: coo2csr -> csr_sort_indices -> csr_sum_duplicates -> csr_eliminate_zeros"
    Bp = Bp - 1
    Bj = Bj - 1
    write(*,*) "coo2csr_canonical: convert to 0-based index for C use"
    write(*,*) "coo2csr_canonical: Normal End"
    
end subroutine

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

subroutine csr_sort_indices(nnz, Bp_size, Bp, Bj, Bx)
    ! Sort CSR column indices inplace
    implicit none
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

subroutine csr_sum_duplicates(nnz, Bp_size, Bp, Bj, Bx)
        ! Sum together duplicate column entries in each row of CSR matrix A
        ! The column indicies within each row must be in sorted order.
        ! Explicit zeros are retained.
        ! Ap, Aj, and Ax will be modified *inplace*
   	    integer ( kind = 4 ) nnz, Bp_size
        integer ( kind = 4 ), dimension(Bp_size) :: Bp
	      integer ( kind = 4 ), dimension(  nnz  ) :: Bj
        double precision    , dimension(  nnz  ) :: Bx
        integer :: r1, r2, i, j, jj, kk
        double precision :: x
        kk = 1
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
                Bj(kk) = j
                Bx(kk) = x
                kk = kk + 1
            end do
            Bp(i+1) = kk
        end do
end subroutine
        
subroutine csr_eliminate_zeros(nnz, Bp_size, Bp_temp, Bj_temp, Bx_temp, Bp, Bj, Bx)
    ! Sort CSR column indices inplace
    implicit none
    integer ( kind = 4 ) nnz, Bp_size, nnz_real, zero_count, k, i
    integer ( kind = 4 ), dimension(Bp_size) :: Bp_temp
    integer ( kind = 4 ), dimension(  nnz  ) :: Bj_temp
    double precision    , dimension(  nnz  ) :: Bx_temp
    integer ( kind = 4 ), dimension(Bp_size) :: Bp
    integer ( kind = 4 ), dimension(  nnz  ) :: Bj
    double precision    , dimension(  nnz  ) :: Bx
    Bp = Bp_temp
    nnz_real = 1
    zero_count = 0
    do i = 1, Bp_size - 1
        do k = Bp_temp(i), Bp_temp(i+1) - 1
            if(Bx_temp(k) .ne. 0.0 ) then 
                Bj( nnz_real ) = Bj_temp( k )
                Bx( nnz_real ) = Bx_temp( k )
                nnz_real = nnz_real + 1
            else 
                zero_count = zero_count + 1
            end if 
        end do
        Bp( i + 1 ) =  Bp( i + 1 ) - zero_count
    end do
end subroutine

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
    do i = 1, Ap_size - 1
        y(i) = dot_product(Ax(Ap(i):Ap(i+1)-1), x(Aj(Ap(i):Ap(i+1)-1)))
    end do
    !$omp end do
    !$omp end parallel
    
end subroutine csr_matvec

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


subroutine csrcoo ( nrow, job, nzmax, a, ja, ia, nnz, ao, ir, jc, ierr )

    !*****************************************************************************80
    !
    !! CSRCOO converts Compressed Sparse Row to Coordinate format.
    !
    !  Discussion:
    !
    !   This routine converts a matrix that is stored in row general sparse 
    !   A, JA, IA format into coordinate format AO, IR, JC. 
    !
    !  Modified:
    !
    !    07 January 2004
    !
    !  Author:
    !
    !    Youcef Saad
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) NROW, the row dimension of the matrix.
    !
    ! job   = integer ( kind = 4 ) serving as a job indicator.
    !         if job = 1 fill in only the array ir, ignore jc, and ao.
    !         if job = 2 fill in ir, and jc but not ao
    !         if job = 3 fill in everything.
    !         The reason why these options are provided is that on return
    !         ao and jc are the same as a, ja. So when job = 3, a and ja are
    !         simply copied into ao, jc.  When job=2, only jc and ir are
    !         returned. With job=1 only the array ir is returned. Moreover,
    !         the algorithm is in place:
    !           call csrcoo (nrow,1,nzmax,a,ja,ia,nnz,a,ia,ja,ierr)
    !         will write the output matrix in coordinate format on a, ja,ia.
    !         (Important: note the order in the output arrays a, ja, ia. )
    !         i.e., ao can be the same as a, ir can be the same as ia
    !         and jc can be the same as ja.
    !
    !    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
    !    Compressed Sparse Row format.
    !
    ! nzmax = length of space available in ao, ir, jc.
    !         the code will stop immediatly if the number of
    !         nonzero elements found in input matrix exceeds nzmax.
    !
    ! on return:
    !-
    ! ao, ir, jc = matrix in coordinate format.
    !
    ! nnz        = number of nonzero elements in matrix.
    !
    ! ierr       = integer ( kind = 4 ) error indicator.
    !         ierr == 0 means normal retur
    !         ierr == 1 means that the the code stopped
    !         because there was no space in ao, ir, jc
    !         (according to the value of  nzmax).
    !
      implicit none
    
      integer ( kind = 4 ) nrow
    
      real ( kind = 8 ) a(*)
      real ( kind = 8 ) ao(*)
      integer ( kind = 4 ) i
      integer ( kind = 4 ) ia(nrow+1)
      integer ( kind = 4 ) ierr
      integer ( kind = 4 ) ir(*)
      integer ( kind = 4 ) ja(*)
      integer ( kind = 4 ) jc(*)
      integer ( kind = 4 ) job
      integer ( kind = 4 ) k
      integer ( kind = 4 ) k1
      integer ( kind = 4 ) k2
      integer ( kind = 4 ) nnz
      integer ( kind = 4 ) nzmax
    
      ierr = 0
      nnz = ia(nrow+1)-1
    
      if ( nzmax < nnz ) then
        ierr = 1
        return
      end if
    
      if ( 3 <= job ) then
        ao(1:nnz) = a(1:nnz)
      end if
    
      if ( 2 <= job ) then
        jc(1:nnz) = ja(1:nnz)
      end if
    !
    !  Copy backward.
    !
      do i = nrow, 1, -1
        k1 = ia(i+1) - 1
        k2 = ia(i)
        do k = k1, k2, -1
          ir(k) = i
        end do
      end do
    
      return
    end


    subroutine csrcsc ( n, job, ipos, a, ja, ia, ao, jao, iao )
 
        !*****************************************************************************80
        !
        !! CSRCSC converts Compressed Sparse Row to Compressed Sparse Column.
        !
        !  Discussion:
        !
        !    This is essentially a transposition operation.  
        !
        !    It is NOT an in-place algorithm.
        !
        !    This routine transposes a matrix stored in a, ja, ia format.
        !
        !  Modified:
        !
        !    07 January 2004
        !
        !  Author:
        !
        !    Youcef Saad
        !
        !  Parameters:
        !
        !    Input, integer ( kind = 4 ) N, the order of the matrix.
        !
        !    Input, integer ( kind = 4 ) JOB, indicates whether or not to fill the values of the
        !    matrix AO or only the pattern (IA, and JA).  Enter 1 for yes.
        !
        ! ipos  = starting position in ao, jao of the transposed matrix.
        !         the iao array takes this into account (thus iao(1) is set to ipos.)
        !         Note: this may be useful if one needs to append the data structure
        !         of the transpose to that of A. In this case use
        !                call csrcsc (n,1,n+2,a,ja,ia,a,ja,ia(n+2))
        !        for any other normal usage, enter ipos=1.
        !
        !    Input, real A(*), integer ( kind = 4 ) JA(*), IA(N+1), the matrix in CSR
        !    Compressed Sparse Row format.
        !
        !    Output, real AO(*), JAO(*), IAO(N+1), the matrix in CSC
        !    Compressed Sparse Column format.
        !
          implicit none
        
          integer ( kind = 4 ) n
        
          real ( kind = 8 ) a(*)
          real ( kind = 8 ) ao(*)
          integer ( kind = 4 ) i
          integer ( kind = 4 ) ia(n+1)
          integer ( kind = 4 ) iao(n+1)
          integer ( kind = 4 ) ipos
          integer ( kind = 4 ) j
          integer ( kind = 4 ) ja(*)
          integer ( kind = 4 ) jao(*)
          integer ( kind = 4 ) job
          integer ( kind = 4 ) k
          integer ( kind = 4 ) next
        !
        !  Compute lengths of rows of A'.
        !
          iao(1:n+1) = 0
        
          do i = 1, n
            do k = ia(i), ia(i+1)-1
              j = ja(k) + 1
              iao(j) = iao(j) + 1
            end do
          end do
        !
        !  Compute pointers from lengths.
        !
          iao(1) = ipos
          do i = 1, n
            iao(i+1) = iao(i) + iao(i+1)
          end do
        !
        !  Do the actual copying.
        !
          do i = 1, n
            do k = ia(i), ia(i+1)-1
              j = ja(k)
              next = iao(j)
              if ( job == 1 ) then
                ao(next) = a(k)
              end if
              jao(next) = i
              iao(j) = next + 1
            end do
          end do
        !
        !  Reshift IAO and leave.
        !
          do i = n, 1, -1
            iao(i+1) = iao(i)
          end do
          iao(1) = ipos
        
          return
        end

end module
