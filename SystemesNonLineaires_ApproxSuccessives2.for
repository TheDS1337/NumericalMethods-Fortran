      PROGRAM SystemesNonLineaires_ApproxSuccessive

      REAL Xanc(1:2), Xnouv(1:2)
      REAL epsilon, norm

      INTEGER i, j    

      Xanc(1) = 1.0
      Xanc(2) = 1.0

      epsilon = 1.E-4

      i = 0

      DO 1 i = 1, 100, 1
        norm = 0.

        Xnouv(1) = SIN(Xanc(1) + Xanc(2)) 
        Xnouv(2) = COS(Xanc(1) + Xanc(2)) 

        DO 2 j = 1, 2, 1
          norm = norm + ((Xnouv(j) - Xanc(j)) ** 2)
 2      CONTINUE

        Xanc(1) = Xnouv(1)
        Xanc(2) = Xnouv(2)

        IF (epsilon .GE. SQRT(norm)) GOTO 3
 1    CONTINUE

 3    WRITE(*, *) i, Xnouv(1), Xnouv(2)

      STOP
      END
