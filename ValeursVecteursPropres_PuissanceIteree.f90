program ValeursVecteursPropres_PuissanceIteree

implicit none

integer, parameter:: dim = 3
integer, parameter:: max_iterations = 100

double precision, dimension(dim, dim, dim):: A
double precision, dimension(dim, dim):: V, W
double precision, dimension(dim):: P, R, U, X
double precision, dimension(1:2):: D
double precision, dimension(dim):: Lambda

double precision:: epsilon, mdl

integer:: i, j, k, l, y, z

epsilon = 1D-400

data A(1, :, :) / -1., 3., -3., 3., -1., 3., -3., 3., 5. /
data X / 1., 2., 1. /

z = 0

do i = 1, dim, 1
    V(i, :) = X(:)
    y = 0
    
    do while (y < max_iterations)
        y = y + 1
        z = z + 1
        
        do j = 1, dim, 1
            P(j) = 0.
            
            do k = 1, dim, 1
                P(j) = P(j) + A(i, j, k) * V(i, k)
            end do
        end do

        l = 1
        
        do j = 1, dim, 1
            if (V(i, j) == 0.) then
                R(j) = 0.
                cycle
            end if

            R(j) = P(j) / V(i, j)
            D(1) = (Lambda(i) - R(j)) * (P(j) - V(i, j))
            D(2) = (Lambda(i) - R(j - 1)) * (P(j - 1) - V(i, j - 1))
            
            ! Le rapport le plus loin
            if ((j > 1) .AND. (abs(D(1)) > abs(D(2)))) then
                l = j
            end if 
        end do

        if (abs(Lambda(i) - R(l)) <= epsilon) then
            Lambda(i) = R(l)
            V(i, :) = P(:)
            exit
        end if 

        Lambda(i) = R(l)
        V(i, :) = P(:)
    end do
    
    if (i /= dim) then
        U(:) = V(i, :)
         
        do j = 1, dim, 1
            do k = 1, dim, 1
                W(j, k) = U(j) * U(k)
            end do
        end do
        
        A(i + 1, :, :) = A(i, :, :) - (Lambda(i) / dot_product(U, U)) * W(:, :)
    end if
end do

write(*, *) 'Iterations: ', z
write(*, *) '					'
write(*, *) '					'

do i = 1, dim, 1		
    write(*, *) 'Valeur Propre: ', Lambda(i)
    write(*, *) 'Vectur Propre: '

    mdl = 1 / sqrt(dot_product(V(i, :), V(i, :)))
    do j = 1, dim, 1
        write(*, *) '		', mdl * V(i, j)
    end do 

    write(*, *) '					'
end do

stop

end program ValeursVecteursPropres_PuissanceIteree