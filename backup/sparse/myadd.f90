module myadd
contains

subroutine subadd( n, a, b, c  )
	implicit none 
	integer n
	integer, dimension(n):: a, b, c
	
	c = a + b;

end subroutine 

subroutine submutil( n, a, b, c  )
	implicit none 
	integer n
	integer, dimension(n):: a, b, c	
	c = a * b;

end subroutine 

end module 



