      SUBROUTINE GRDFCN(N,X,G,NPROB)                                    00000010
      INTEGER N,NPROB
      DOUBLE PRECISION X(N),G(N)
C     **********
C
C     SUBROUTINE GRDFCN
C
C     THIS SUBROUTINE DEFINES THE GRADIENT VECTORS OF EIGHTEEN
C     NONLINEAR UNCONSTRAINED MINIMIZATION PROBLEMS. THE PROBLEM
C     DIMENSIONS ARE AS DESCRIBED IN THE PROLOGUE COMMENTS OF OBJFCN.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE GRDFCN(N,X,G,NPROB)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       X IS AN INPUT ARRAY OF LENGTH N.
C
C       G IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE COMPONENTS
C         OF THE GRADIENT VECTOR OF THE NPROB OBJECTIVE FUNCTION
C         EVALUATED AT X.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... DABS,DATAN,DCOS,DEXP,DLOG,DSIGN,DSIN,
C                            DSQRT
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER I,IEV,IVAR,J
      DOUBLE PRECISION AP,ARG,BP,C2PDM6,CP0001,CP1,CP2,CP25,CP5,C1P5,
     1                 C2P25,C2P625,C3P5,C19P8,C20P2,C25,C29,C100,
     2                 C180,C200,C10000,C1PD6,D1,D2,EIGHT,FIFTY,FIVE,
     3                 FOUR,ONE,R,S1,S2,S3,T,T1,T2,T3,TEN,TH,THREE,
     4                 TPI,TWENTY,TWO,ZERO
      DOUBLE PRECISION FVEC(50),Y(15)
      DOUBLE PRECISION DFLOAT
      DATA ZERO,ONE,TWO,THREE,FOUR,FIVE,EIGHT,TEN,TWENTY,FIFTY
     1     /0.0D0,1.0D0,2.0D0,3.0D0,4.0D0,5.0D0,8.0D0,1.0D1,2.0D1,
     2      5.0D1/
      DATA C2PDM6,CP0001,CP1,CP2,CP25,CP5,C1P5,C2P25,C2P625,C3P5,
     1     C19P8,C20P2,C25,C29,C100,C180,C200,C10000,C1PD6
     2     /2.0D-6,1.0D-4,1.0D-1,2.0D-1,2.5D-1,5.0D-1,1.5D0,2.25D0,
     3      2.625D0,3.5D0,1.98D1,2.02D1,2.5D1,2.9D1,1.0D2,1.8D2,2.0D2,
     4      1.0D4,1.0D6/
      DATA AP,BP /1.0D-5,1.0D0/
      DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),Y(9),Y(10),Y(11),
     1     Y(12),Y(13),Y(14),Y(15)
     2     /9.0D-4,4.4D-3,1.75D-2,5.4D-2,1.295D-1,2.42D-1,3.521D-1,
     3      3.989D-1,3.521D-1,2.42D-1,1.295D-1,5.4D-2,1.75D-2,4.4D-3,
     4      9.0D-4/
      DFLOAT(IVAR) = IVAR
C
C     GRADIENT ROUTINE SELECTOR.
C
      GO TO (10,20,50,70,80,100,130,190,220,260,270,290,310,350,370,
     1       390,400,410), NPROB
C
C     HELICAL VALLEY FUNCTION.
C
   10 CONTINUE
      TPI = EIGHT*DATAN(ONE)
      TH = DSIGN(CP25,X(2))
      IF (X(1) .GT. ZERO) TH = DATAN(X(2)/X(1))/TPI
      IF (X(1) .LT. ZERO) TH = DATAN(X(2)/X(1))/TPI + CP5
      ARG = X(1)**2 + X(2)**2
      R = DSQRT(ARG)
      T = X(3) - TEN*TH
      S1 = TEN*T/(TPI*ARG)
      G(1) = C200*(X(1) - X(1)/R + X(2)*S1)
      G(2) = C200*(X(2) - X(2)/R - X(1)*S1)
      G(3) = TWO*(C100*T + X(3))
      GO TO 490
C
C     BIGGS EXP6 FUNCTION.
C
   20 CONTINUE
      DO 30 J = 1, N
         G(J) = ZERO
   30    CONTINUE
      DO 40 I = 1, 13
         D1 = DFLOAT(I)/TEN
         D2 = DEXP(-D1) - FIVE*DEXP(-TEN*D1) + THREE*DEXP(-FOUR*D1)
         S1 = DEXP(-D1*X(1))
         S2 = DEXP(-D1*X(2))
         S3 = DEXP(-D1*X(5))
         T = X(3)*S1 - X(4)*S2 + X(6)*S3 - D2
         TH = D1*T
         G(1) = G(1) - S1*TH
         G(2) = G(2) + S2*TH
         G(3) = G(3) + S1*T
         G(4) = G(4) - S2*T
         G(5) = G(5) - S3*TH
         G(6) = G(6) + S3*T
   40    CONTINUE
      G(1) = TWO*X(3)*G(1)
      G(2) = TWO*X(4)*G(2)
      G(3) = TWO*G(3)
      G(4) = TWO*G(4)
      G(5) = TWO*X(6)*G(5)
      G(6) = TWO*G(6)
      GO TO 490
C
C     GAUSSIAN FUNCTION.
C
   50 CONTINUE
      G(1) = ZERO
      G(2) = ZERO
      G(3) = ZERO
      DO 60 I = 1, 15
         D1 = CP5*DFLOAT(I-1)
         D2 = C3P5 - D1 - X(3)
         ARG = -CP5*X(2)*D2**2
         R = DEXP(ARG)
         T = X(1)*R - Y(I)
         S1 = R*T
         S2 = D2*S1
         G(1) = G(1) + S1
         G(2) = G(2) - D2*S2
         G(3) = G(3) + S2
   60    CONTINUE
      G(1) = TWO*G(1)
      G(2) = X(1)*G(2)
      G(3) = TWO*X(1)*X(2)*G(3)
      GO TO 490
C
C     POWELL BADLY SCALED FUNCTION.
C
   70 CONTINUE
      T1 = C10000*X(1)*X(2) - ONE
      S1 = DEXP(-X(1))
      S2 = DEXP(-X(2))
      T2 = S1 + S2 - ONE - CP0001
      G(1) = TWO*(C10000*X(2)*T1 - S1*T2)
      G(2) = TWO*(C10000*X(1)*T1 - S2*T2)
      GO TO 490
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
   80 CONTINUE
      G(1) = ZERO
      G(2) = ZERO
      G(3) = ZERO
      DO 90 I = 1, 10
         D1 = DFLOAT(I)
         D2 = D1/TEN
         S1 = DEXP(-D2*X(1))
         S2 = DEXP(-D2*X(2))
         S3 = DEXP(-D2) - DEXP(-D1)
         T = S1 - S2 - S3*X(3)
         TH = D2*T
         G(1) = G(1) - S1*TH
         G(2) = G(2) + S2*TH
         G(3) = G(3) - S3*T
   90    CONTINUE
      G(1) = TWO*G(1)
      G(2) = TWO*G(2)
      G(3) = TWO*G(3)
      GO TO 490
C
C     VARIABLY DIMENSIONED FUNCTION.
C
  100 CONTINUE
      T1 = ZERO
      DO 110 J = 1, N
         T1 = T1 + DFLOAT(J)*(X(J) - ONE)
  110    CONTINUE
      T = T1*(ONE + TWO*T1**2)
      DO 120 J = 1, N
         G(J) = TWO*(X(J) - ONE + DFLOAT(J)*T)
  120    CONTINUE
      GO TO 490
C
C     WATSON FUNCTION.
C
  130 CONTINUE
      DO 140 J = 1, N
         G(J) = ZERO
  140    CONTINUE
      DO 180 I = 1, 29
         D1 = DFLOAT(I)/C29
         S1 = ZERO
         D2 = ONE
         DO 150 J = 2, N
            S1 = S1 + DFLOAT(J-1)*D2*X(J)
            D2 = D1*D2
  150       CONTINUE
         S2 = ZERO
         D2 = ONE
         DO 160 J = 1, N
            S2 = S2 + D2*X(J)
            D2 = D1*D2
  160       CONTINUE
         T = S1 - S2**2 - ONE
         S3 = TWO*D1*S2
         D2 = TWO/D1
         DO 170 J = 1, N
            G(J) = G(J) + D2*(DFLOAT(J-1) - S3)*T
            D2 = D1*D2
  170       CONTINUE
  180    CONTINUE
      T1 = X(2) - X(1)**2 - ONE
      G(1) = G(1) + X(1)*(TWO - FOUR*T1)
      G(2) = G(2) + TWO*T1
      GO TO 490
C
C     PENALTY FUNCTION I.
C
  190 CONTINUE
      T1 = -CP25
      DO 200 J = 1, N
         T1 = T1 + X(J)**2
  200    CONTINUE
      D1 = TWO*AP
      TH = FOUR*BP*T1
      DO 210 J = 1, N
         G(J) = D1*(X(J) - ONE) + X(J)*TH
  210    CONTINUE
      GO TO 490
C
C     PENALTY FUNCTION II.
C
  220 CONTINUE
      T1 = -ONE
      DO 230 J = 1, N
         T1 = T1 + DFLOAT(N-J+1)*X(J)**2
  230    CONTINUE
      D1 = DEXP(CP1)
      D2 = ONE
      TH = FOUR*BP*T1
      DO 250 J = 1, N
         G(J) = DFLOAT(N-J+1)*X(J)*TH
         S1 = DEXP(X(J)/TEN)
         IF (J .EQ. 1) GO TO 240
         S3 = S1 + S2 - D2*(D1 + ONE)
         G(J) = G(J) + AP*S1*(S3 + S1 - ONE/D1)/FIVE
         G(J-1) = G(J-1) + AP*S2*S3/FIVE
  240    CONTINUE
         S2 = S1
         D2 = D1*D2
  250    CONTINUE
      G(1) = G(1) + TWO*BP*(X(1) - CP2)
      GO TO 490
C
C     BROWN BADLY SCALED FUNCTION.
C
  260 CONTINUE
      T1 = X(1) - C1PD6
      T2 = X(2) - C2PDM6
      T3 = X(1)*X(2) - TWO
      G(1) = TWO*(T1 + X(2)*T3)
      G(2) = TWO*(T2 + X(1)*T3)
      GO TO 490
C
C     BROWN AND DENNIS FUNCTION.
C
  270 CONTINUE
      G(1) = ZERO
      G(2) = ZERO
      G(3) = ZERO
      G(4) = ZERO
      DO 280 I = 1, 20
         D1 = DFLOAT(I)/FIVE
         D2 = DSIN(D1)
         T1 = X(1) + D1*X(2) - DEXP(D1)
         T2 = X(3) + D2*X(4) - DCOS(D1)
         T = T1**2 + T2**2
         S1 = T1*T
         S2 = T2*T
         G(1) = G(1) + S1
         G(2) = G(2) + D1*S1
         G(3) = G(3) + S2
         G(4) = G(4) + D2*S2
  280    CONTINUE
      G(1) = FOUR*G(1)
      G(2) = FOUR*G(2)
      G(3) = FOUR*G(3)
      G(4) = FOUR*G(4)
      GO TO 490
C
C     GULF RESEARCH AND DEVELOPMENT FUNCTION.
C
  290 CONTINUE
      G(1) = ZERO
      G(2) = ZERO
      G(3) = ZERO
      D1 = TWO/THREE
      DO 300 I = 1, 99
         ARG = DFLOAT(I)/C100
         R = DABS((-FIFTY*DLOG(ARG))**D1+C25-X(2))
         T1 = R**X(3)/X(1)
         T2 = DEXP(-T1)
         T = T2 - ARG
         S1 = T1*T2*T
         G(1) = G(1) + S1
         G(2) = G(2) + S1/R
         G(3) = G(3) - S1*DLOG(R)
  300    CONTINUE
      G(1) = TWO*G(1)/X(1)
      G(2) = TWO*X(3)*G(2)
      G(3) = TWO*G(3)
      GO TO 490
C
C     TRIGONOMETRIC FUNCTION.
C
  310 CONTINUE
      S1 = ZERO
      DO 320 J = 1, N
         G(J) = DCOS(X(J))
         S1 = S1 + G(J)
  320    CONTINUE
      S2 = ZERO
      DO 330 J = 1, N
         TH = DSIN(X(J))
         T = DFLOAT(N+J) - TH - S1 - DFLOAT(J)*G(J)
         S2 = S2 + T
         G(J) = (DFLOAT(J)*TH - G(J))*T
  330    CONTINUE
      DO 340 J = 1, N
         G(J) = TWO*(G(J) + DSIN(X(J))*S2)
  340    CONTINUE
      GO TO 490
C
C     EXTENDED ROSENBROCK FUNCTION.
C
  350 CONTINUE
      DO 360 J = 1, N, 2
         T1 = ONE - X(J)
         G(J+1) = C200*(X(J+1) - X(J)**2)
         G(J) = -TWO*(X(J)*G(J+1) + T1)
  360    CONTINUE
      GO TO 490
C
C     EXTENDED POWELL FUNCTION.
C
  370 CONTINUE
      DO 380 J = 1, N, 4
         T = X(J) + TEN*X(J+1)
         T1 = X(J+2) - X(J+3)
         S1 = FIVE*T1
         T2 = X(J+1) - TWO*X(J+2)
         S2 = FOUR*T2**3
         T3 = X(J) - X(J+3)
         S3 = TWENTY*T3**3
         G(J) = TWO*(T + S3)
         G(J+1) = TWENTY*T + S2
         G(J+2) = TWO*(S1 - S2)
         G(J+3) = -TWO*(S1 + S3)
  380    CONTINUE
      GO TO 490
C
C     BEALE FUNCTION.
C
  390 CONTINUE
      S1 = ONE - X(2)
      T1 = C1P5 - X(1)*S1
      S2 = ONE - X(2)**2
      T2 = C2P25 - X(1)*S2
      S3 = ONE - X(2)**3
      T3 = C2P625 - X(1)*S3
      G(1) = -TWO*(S1*T1 + S2*T2 + S3*T3)
      G(2) = TWO*X(1)*(T1 + X(2)*(TWO*T2 + THREE*X(2)*T3))
      GO TO 490
C
C     WOOD FUNCTION.
C
  400 CONTINUE
      S1 = X(2) - X(1)**2
      S2 = ONE - X(1)
      S3 = X(2) - ONE
      T1 = X(4) - X(3)**2
      T2 = ONE - X(3)
      T3 = X(4) - ONE
      G(1) = -TWO*(C200*X(1)*S1 + S2)
      G(2) = C200*S1 + C20P2*S3 + C19P8*T3
      G(3) = -TWO*(C180*X(3)*T1 + T2)
      G(4) = C180*T1 + C20P2*T3 + C19P8*S3
      GO TO 490
C
C     CHEBYQUAD FUNCTION.
C
  410 CONTINUE
      DO 420 I = 1, N
         FVEC(I) = ZERO
  420    CONTINUE
      DO 440 J = 1, N
         T1 = ONE
         T2 = TWO*X(J) - ONE
         T = TWO*T2
         DO 430 I = 1, N
            FVEC(I) = FVEC(I) + T2
            TH = T*T2 - T1
            T1 = T2
            T2 = TH
  430       CONTINUE
  440    CONTINUE
      D1 = ONE/DFLOAT(N)
      IEV = -1
      DO 450 I = 1, N
         FVEC(I) = D1*FVEC(I)
         IF (IEV .GT. 0) FVEC(I) = FVEC(I) + ONE/(DFLOAT(I)**2 - ONE)
         IEV = -IEV
  450    CONTINUE
      DO 470 J = 1, N
         G(J) = ZERO
         T1 = ONE
         T2 = TWO*X(J) - ONE
         T = TWO*T2
         S1 = ZERO
         S2 = TWO
         DO 460 I = 1, N
            G(J) = G(J) + FVEC(I)*S2
            TH = FOUR*T2 + T*S2 - S1
            S1 = S2
            S2 = TH
            TH = T*T2 - T1
            T1 = T2
            T2 = TH
  460       CONTINUE
  470    CONTINUE
      D2 = TWO*D1
      DO 480 J = 1, N
         G(J) = D2*G(J)
  480    CONTINUE
  490 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE GRDFCN.
C
      END
