program SystemesNonLineaires_ApproxSuccessive

implicit none

integer, parameter:: dim = 2, iterations = 100

real, dimension(0:iterations, dim):: X
real, dimension(dim):: diff

real:: epsilon, norm
integer:: i, j, l      

X(0, 1) = 1.
X(0, 2) = 1.

epsilon = 1.E-4

do i = 1, iterations, 1    
    l = i - 1

    X(i, 1) = sin(X(l, 1) + X(l, 2)) 
    X(i, 2) = cos(X(l, 1) + X(l, 2))

    diff = X(i, :) - X(l, :)
    norm = sqrt(dot_product(diff, diff))

    if (epsilon >= norm) then 
        exit
    end if
end do

do j = 1, dim, 1
    write(*, *) j, X(i, j)
end do

stop

end program SystemesNonLineaires_ApproxSuccessive