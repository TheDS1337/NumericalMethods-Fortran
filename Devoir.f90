program Devoir

implicit none

integer, parameter:: iterations_x = 4
integer, parameter:: iterations_t = 2
integer, parameter:: dim = iterations_x - 1
double precision, PARAMETER:: PI = 4.D0 * datan(1.0D0)

real, dimension(0:iterations_x, 0:iterations_t):: F, E
real, dimension(dim, dim):: A
real, dimension(dim):: B
real, dimension(dim, dim):: P

real:: epsilon = 1.E-4
real:: alpha = 0.16
real:: dx = 0.5
real:: dt = 0.5
real:: r
double precision:: lambda, prod

integer:: i, j, k, l, m

r = alpha * dt / (dx**2)

! Valeurs Analytiques
do k = 0, iterations_x, 1
    do l = 0, iterations_t, 1
        E(k, l) = 0.

        do m = 0, 1000000, 2
            lambda = (m + 0.5) * PI
            prod = 200 / (lambda**2)
            prod = prod * sin(lambda * k * dx) * exp(-((alpha * lambda)**2) * l * dt)
            E(k, l) = E(k, l) + prod
        end do
        
        do m = 1, 1000000, 2
            lambda = (m + 0.5) * PI
            prod = -200 / (lambda**2)
            prod = prod * sin(lambda * k * dx) * exp(-((alpha * lambda)**2) * l * dt)
            E(k, l) = E(k, l) + prod
        end do
        
    end do
end do

! Valeurs Numérique (Initialisation)
do k = 0, iterations_x, 1
    do l = 0, iterations_t, 1
        if (k == 0 .OR. k == iterations_x) then
            F(k, l) = 0.
        else if (l == 0) then
            if (k <= iterations_x / 2 ) then
                F(k, l) = 100 * dx * k
            else 
                F(k, l) = 100 * (2 - dx * k)
            end if                
        end if 
    end do
end do

do j = 0, iterations_t - 1, 1    
    do i = 1, dim, 1
        if (i == 1) then 
            A(1, 1) = 1. + 2 * r
            A(1, 2) = -r
            A(1, 3) = 0.
            A(2, 1) = -r
            A(2, 2) = 1. + 2 * r
            A(2, 3) = -r
            A(3, 1) = 0.
            A(3, 2) = -r
            A(3, 3) = 1. + 2 * r

            B(1) = F(1, j)
            B(2) = F(2, j)
            B(3) = F(3, j)     
        end if 

        ! Calcule de P(n)
        do k = 1, dim, 1
            do l = 1, dim, 1
                if (l == i) then
                    if (k == l) then
                        P(k, l) = 1. / A(k, l)
                    else
                        P(k, l) = - A(k, l) / A(l, l)
                    end if
                else
                    if (k == l) then
                        P(k, l) = 1.
                    else
                        P(k, l) = 0.
                    end if
                end if
            end do
        end do

        A(:, :) = matmul(P(:, :), A(:, :))
        B(:) = matmul(P(:, :), B(:))
    end do

    do m = 1, iterations_x - 1, 1 
        F(m, j + 1) = B(m)        
    end do     
end do
                                                                          
write(*, *) '       x           ', '      t          ', '   f(x, t) num   ', '        f(x, t)      ', '       erreur'
do i = 0, iterations_x, 1
    do j = 0, iterations_t, 1
        write(*, *) i * dx, '|', j * dt, '|', F(i, j), '|', E(i, j), '|', E(i, j) - F(i, j)
    end do
end do

stop

end program Devoir