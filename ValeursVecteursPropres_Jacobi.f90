program ValeursVecteursPropres_Jacobi

implicit none

integer, parameter:: dim = 4
integer, parameter:: max_iterations = 100
double precision, PARAMETER:: PI = 4.D0 * datan(1.0D0)

real, dimension(dim, dim):: A, T
real, dimension(1:max_iterations, dim, dim):: U

real:: epsilon, diff, theta, theta_cos, theta_sin

integer:: i, j, k, g, l, z

logical:: solution_found = .false.

epsilon = 1.E-4

data A /4.0, -30.0, 60.0, -35.0, -30.0, 300.0, -675.0, 420.0, 60.0, -675.0, 1620.0, -1050.0, -35.0, 420.0, -1050.0, 700.0/
    
do i = 1, dim, 1
    do j = 1, dim, 1
        if (A(i, j) /= A(j, i)) then
            write(*, *) 'La matrice n''est pas symétrique'	    
            stop
        end if
    end do  
end do

z = 1

do while ((solution_found .eqv. .false.) .AND. (z < max_iterations))
    do g = 1, dim - 1, 1
        do l = g + 1, dim, 1
            ! Des fois, on ne se trouve pas dans le cas A(g, l) = 0.0
            if (abs(A(g, l)) <= epsilon) then
                cycle
            endif

            diff = A(g, g) - A(l, l)

            if (diff == 0) then
                theta = PI / 4
            else
                theta = atan(2 * A(g, l) / diff) / 2
            endif

            theta_cos = cos(theta)
            theta_sin = sin(theta)

            ! Calcule de U(n)
            do i = 1, dim, 1
                do j = 1, dim, 1
                    if (j == i) then
                        U(z, i, j) = 1.
                    else
                        U(z, i, j) = 0.
                    end if	    			
                end do 
            end do

            U(z, g, g) = theta_cos
            U(z, g, l) = theta_sin
            U(z, l, g) = theta_sin
            U(z, l, l) = -theta_cos

            ! Calcule de U(n) * A(n) * U(n)        
            A(:, :) = matmul(U(z, :, :), matmul(A(:, :), U(z, :, :)))        

            ! Préservation de la symétrie... (des fois, les valeurs Aij et Aji sont différentes avec une difference de l'ordre de ~10^-7)
            do i = 1, dim - 1, 1
                do j = i + 1, dim, 1				
                    A(j, i) = A(i, j)
                end do 
            end do

            ! Calcule de T = T * U(n) (pour le produit U(1)...x...U(n))
            if (z == 1) then
                T(:, :) = U(z, :, :)
            else
                T(:, :) = matmul(T(:, :), U(z, :, :))  
            end if
        
            z = z + 1
        end do   
    end do

    solution_found = .true.

    do g = 1, dim - 1, 1
        do l = g + 1, dim, 1
            if (abs(A(g, l)) > epsilon) then
                solution_found = .false.
                exit
            end if
        end do

        if (solution_found .eqv. .false.) then             
            exit
        end if	
    end do
end do

if (solution_found .eqv. .true.) then
    z = z - 1

    write(*, *) 'Iterations: ', z
    write(*, *) '					'
    write(*, *) '					'

    do i = 1, dim, 1		
        write(*, *) 'Valeur Propre: ', A(i, i)
        write(*, *) 'Vectur Propre: '

        do j = 1, dim, 1							
            write(*, *) '		', T(j, i)
        end do 

        write(*, *) '					'
    end do
else
    write(*, *) 'La solution n''était pas trouvée'
end if

stop

end program ValeursVecteursPropres_Jacobi