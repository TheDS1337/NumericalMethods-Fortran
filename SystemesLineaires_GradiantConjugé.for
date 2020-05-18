      PROGRAM SystemesLineaires_GradiantConjuge

      REAL A(1:2, 1:2), B(1:2), X(0:100, 1:2)
      REAL Alpha(1:100), Beta(1:100), R(1:100, 1:2), P(1:100, 1:2)
      REAL epsilon, norm

!     Alpha = - G / F, Beta = - G / F
      REAL G, F

!     Le vector AP(n)
      REAL H(1:2), M(1:2)

      INTEGER i, j, k, l

!     Definition des matrices A, B et X(0)
      DATA A /2., -1., -1., 2./
      DATA B /1., 0./

      X(0, 1) = 0.
      X(0, 2) = 0.

      epsilon = 1.E-4

      i = 0
      j = 0
      k = 0

      DO 1 i = 1, 100, 1
        norm = 0.

        G = 0.
        F = 0.

        H(1) = 0.
        H(2) = 0.

        l = i - 1

!       Calcule de R(n) (= -P(n))
        DO 2 j = 1, 2, 1
          R(i, j) = 0.

          DO 3 k = 1, 2, 1
            R(i, j) = R(i, j) + A(j, k) * X(l, k)
 3        CONTINUE  

          R(i, j) = R(i, j) - B(j)      
          P(i, j) = -R(i, j)

          IF (l .GT. 0) P(i, j) = P(i, j) + Beta(l) * P(l, j)
 2      CONTINUE

!       Calcule de Alpha(n)
        DO 4 j = 1, 2, 1

          DO 5 k = 1, 2, 1
            H(j) = H(j) + A(j, k) * P(i, k)
 5        CONTINUE

        G = G + R(i, j) * P(i, j)
        F = F + H(j)    * P(i, j)
 4      CONTINUE

        Alpha(i) = - G / F

        G = 0.
        F = 0.

        H(1) = 0.
        H(2) = 0.

        M(1) = 0.
        M(2) = 0.

!       Calcule de Beta(n)
        DO 8 j = 1, 2, 1

          DO 9 k = 1, 2, 1
            H(j) = H(j) + A(j, k) * P(i, k)
 9        CONTINUE

          DO 10 k = 1, 2, 1
            M(j) = M(j) + A(j, k) * R(i, k)
 10       CONTINUE

        G = G + M(j) * P(i, j)
        F = F + H(j)    * P(i, j)
 8      CONTINUE

        Beta(i) = - G / F

        DO 6 j = 1, 2, 1
          X(i, j) = X(l, j) + Alpha(i) * P(i, j)
          norm = norm + ((X(i, j) - X(l, j)) ** 2)
 6      CONTINUE

        IF (epsilon .GE. SQRT(norm)) GOTO 7
 1    CONTINUE

 7    WRITE(*, *) i, X(i, 1), X(i, 2)

      STOP
      END