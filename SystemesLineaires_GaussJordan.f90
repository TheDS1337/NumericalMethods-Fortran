program SystemesLineaires_GaussJordan

implicit none

integer, parameter:: dim = 3

real, dimension(1:2, dim, dim):: A
real, dimension(dim):: X
real, dimension(1:2, dim):: B
real, dimension(dim, dim):: P

real:: epsilon = 1.E-4
real:: norm

integer:: i, j, k, l, z

X(1) = 0.
X(2) = 0.
X(3) = 0.

do z = 1, 100, 1
	do i = 1, dim, 1
        if (i == 1) then
	        A(1, 1, 1) = -1.
        	A(1, 1, 2) = X(2) / 6.
        	A(1, 1, 3) = - exp(-X(3)) / 6.
            	B(1, 1) = X(1) - ((X(2)**2 + 2. * exp(-X(3))) / 12.)
            
        	A(1, 2, 1) = - 1. / 6.
        	A(1, 2, 2) = -1.
        	A(1, 2, 3) = cos(X(3)) / 6.
        	B(1, 2) = X(2) - ((1. - X(1) + sin(X(3))) / 6.)
            
        	A(1, 3, 1) = X(1) / 3.
        	A(1, 3, 2) = X(2) / 3.
        	A(1, 3, 3) = (X(3) / 3.) - 1.
        	B(1, 3) = X(3) - ((dot_product(X, X)) / 6.)
	    else 
	        do j = 1, dim, 1
    	    		B(1, j) = B(2, j)
    	       		do k = 1, dim, 1
    	            		A(1, j, k) = A(2, j, k)
    	    		end do  
	    	end do
	    end if
	    
		! Calcule de P(n)
		do j = 1, dim, 1
		  	do k = 1, dim, 1
		    	if (k == i) then
			      	if (j == k) then
			        	P(j, k) = 1. / A(1, j, k)
			      	else
			        	P(j, k) = - A(1, j, k) / A(1, k, k)
			      	end if
			    else
			      	if (j == k) then
			        	P(j, k) = 1.
			      	else
			        	P(j, k) = 0.
			      	end if
		    	end if
		 	end do 
		end do

	    ! Calcule de A(n) et B(n)
	    do j = 1, dim, 1
	    	B(2, j) = 0.
	       	do k = 1, dim, 1
	            A(2, j, k) = 0.
	            do l = 1, dim, 1
	              A(2, j, k) = A(2, j, k) + P(j, l) * A(1, l, k)
		    	end do
		    	B(2, j) = B(2, j) + P(j, k) * B(1, k)
	    	end do  
	 	end do
	end do
	
	do j = 1, dim, 1
        X(j) = X(j) + B(2, j)
    end do

	norm = sqrt(B(2, 1)**2 + B(2, 2)**2 + B(2, 3)**2)

	if (norm <= epsilon) then 
        exit
    end if

end do

write(*, *) z

do i = 1, dim, 1
	write(*, *) i, X(i)
end do

stop

end program SystemesLineaires_GaussJordan