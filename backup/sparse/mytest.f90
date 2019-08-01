subroutine mytest( n, a, b, c  )

	use myadd
	implicit none 
	integer n
	integer, dimension(n):: a, b, c
	call subadd(n,a,b,c)
	call submutil(n,a,b,c)
	!c = a + b;

end subroutine 

