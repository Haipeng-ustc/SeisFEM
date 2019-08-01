module tesmodule 
        
        implicit none 
        contains
        
        subroutine tes1(n, a, b,c)
        integer :: n
        integer, dimension(n) ::  a,b,c

        c = a + b
        
        end subroutine tes1

end module tesmodule 
