program SystemesNonLineaires_ApproxSuccessive

real, dimension(0:100, 2):: X
real:: epsilon, norm

integer:: i, j, l      

X(0, 1) = 1.0
X(0, 2) = 1.0

epsilon = 1.E-4

i = 0

DO 1 i = 1, 100, 1
  norm = 0.

  l = i - 1

  X(i, 1) = SIN(X(l, 1) + X(l, 2)) 
  X(i, 2) = COS(X(l, 1) + X(l, 2)) 

  DO 2 j = 1, 2, 1
    norm = norm + ((X(i, j) - X(l, j)) ** 2)
 CONTINUE

  IF (epsilon .GE. SQRT(norm)) GOTO 3
CONTINUE

write(*, *) i, X(i, 1), X(i, 2)

stop

end program SystemesNonLineaires_ApproxSuccessive