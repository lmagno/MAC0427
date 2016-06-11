      SUBROUTINE OBJFCN(N,X,F,NPROB)                                    00000010
      INTEGER N,NPROB
      DOUBLE PRECISION F
      DOUBLE PRECISION X(N)
C     **********
C
C     SUBROUTINE OBJFCN
C
C     THIS SUBROUTINE DEFINES THE OBJECTIVE FUNCTIONS OF EIGHTEEN
C     NONLINEAR UNCONSTRAINED MINIMIZATION PROBLEMS. THE VALUES
C     OF N FOR FUNCTIONS 1,2,3,4,5,10,11,12,16 AND 17 ARE
C     3,6,3,2,3,2,4,3,2 AND 4, RESPECTIVELY.
C     FOR FUNCTION 7, N MAY BE 2 OR GREATER BUT IS USUALLY 6 OR 9.
C     FOR FUNCTIONS 6,8,9,13,14,15 AND 18 N MAY BE VARIABLE,
C     HOWEVER IT MUST BE EVEN FOR FUNCTION 14, A MULTIPLE OF 4 FOR
C     FUNCTION 15, AND NOT GREATER THAN 50 FOR FUNCTION 18.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE OBJFCN(N,X,F,NPROB)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       X IS AN INPUT ARRAY OF LENGTH N.
C
C       F IS AN OUTPUT VARIABLE WHICH CONTAINS THE VALUE OF
C         THE NPROB OBJECTIVE FUNCTION EVALUATED AT X.
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
     1                 C2P25,C2P625,C3P5,C25,C29,C90,C100,C10000,
     2                 C1PD6,D1,D2,EIGHT,FIFTY,FIVE,FOUR,ONE,R,S1,S2,
     3                 S3,T,T1,T2,T3,TEN,TH,THREE,TPI,TWO,ZERO
      DOUBLE PRECISION FVEC(50),Y(15)
      DOUBLE PRECISION DFLOAT
      DATA ZERO,ONE,TWO,THREE,FOUR,FIVE,EIGHT,TEN,FIFTY
     1     /0.0D0,1.0D0,2.0D0,3.0D0,4.0D0,5.0D0,8.0D0,1.0D1,5.0D1/
      DATA C2PDM6,CP0001,CP1,CP2,CP25,CP5,C1P5,C2P25,C2P625,C3P5,C25,
     1     C29,C90,C100,C10000,C1PD6
     2     /2.0D-6,1.0D-4,1.0D-1,2.0D-1,2.5D-1,5.0D-1,1.5D0,2.25D0,
     3      2.625D0,3.5D0,2.5D1,2.9D1,9.0D1,1.0D2,1.0D4,1.0D6/
      DATA AP,BP /1.0D-5,1.0D0/
      DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),Y(9),Y(10),Y(11),
     1     Y(12),Y(13),Y(14),Y(15)
     2     /9.0D-4,4.4D-3,1.75D-2,5.4D-2,1.295D-1,2.42D-1,3.521D-1,
     3      3.989D-1,3.521D-1,2.42D-1,1.295D-1,5.4D-2,1.75D-2,4.4D-3,
     4      9.0D-4/
      DFLOAT(IVAR) = IVAR
C
C     FUNCTION ROUTINE SELECTOR.
C
      GO TO (10,20,40,60,70,90,110,150,170,200,210,230,250,280,300,
     1       320,330,340), NPROB
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
      F = C100*(T**2 + (R - ONE)**2) + X(3)**2
      GO TO 390
C
C     BIGGS EXP6 FUNCTION.
C
   20 CONTINUE
      F = ZERO
      DO 30 I = 1, 13
         D1 = DFLOAT(I)/TEN
         D2 = DEXP(-D1) - FIVE*DEXP(-TEN*D1) + THREE*DEXP(-FOUR*D1)
         S1 = DEXP(-D1*X(1))
         S2 = DEXP(-D1*X(2))
         S3 = DEXP(-D1*X(5))
         T = X(3)*S1 - X(4)*S2 + X(6)*S3 - D2
         F = F + T**2
   30    CONTINUE
      GO TO 390
C
C     GAUSSIAN FUNCTION.
C
   40 CONTINUE
      F = ZERO
      DO 50 I = 1, 15
         D1 = CP5*DFLOAT(I-1)
         D2 = C3P5 - D1 - X(3)
         ARG = -CP5*X(2)*D2**2
         R = DEXP(ARG)
         T = X(1)*R - Y(I)
         F = F + T**2
   50    CONTINUE
      GO TO 390
C
C     POWELL BADLY SCALED FUNCTION.
C
   60 CONTINUE
      T1 = C10000*X(1)*X(2) - ONE
      S1 = DEXP(-X(1))
      S2 = DEXP(-X(2))
      T2 = S1 + S2 - ONE - CP0001
      F = T1**2 + T2**2
      GO TO 390
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
   70 CONTINUE
      F = ZERO
      DO 80 I = 1, 10
         D1 = DFLOAT(I)
         D2 = D1/TEN
         S1 = DEXP(-D2*X(1))
         S2 = DEXP(-D2*X(2))
         S3 = DEXP(-D2) - DEXP(-D1)
         T = S1 - S2 - S3*X(3)
         F = F + T**2
   80    CONTINUE
      GO TO 390
C
C     VARIABLY DIMENSIONED FUNCTION.
C
   90 CONTINUE
      T1 = ZERO
      T2 = ZERO
      DO 100 J = 1, N
         T1 = T1 + DFLOAT(J)*(X(J) - ONE)
         T2 = T2 + (X(J) - ONE)**2
  100    CONTINUE
      F = T2 + T1**2*(ONE + T1**2)
      GO TO 390
C
C     WATSON FUNCTION.
C
  110 CONTINUE
      F = ZERO
      DO 140 I = 1, 29
         D1 = DFLOAT(I)/C29
         S1 = ZERO
         D2 = ONE
         DO 120 J = 2, N
            S1 = S1 + DFLOAT(J-1)*D2*X(J)
            D2 = D1*D2
  120       CONTINUE
         S2 = ZERO
         D2 = ONE
         DO 130 J = 1, N
            S2 = S2 + D2*X(J)
            D2 = D1*D2
  130       CONTINUE
         T = S1 - S2**2 - ONE
         F = F + T**2
  140    CONTINUE
      T1 = X(2) - X(1)**2 - ONE
      F = F + X(1)**2 + T1**2
      GO TO 390
C
C     PENALTY FUNCTION I.
C
  150 CONTINUE
      T1 = -CP25
      T2 = ZERO
      DO 160 J = 1, N
         T1 = T1 + X(J)**2
         T2 = T2 + (X(J) - ONE)**2
  160    CONTINUE
      F = AP*T2 + BP*T1**2
      GO TO 390
C
C     PENALTY FUNCTION II.
C
  170 CONTINUE
      T1 = -ONE
      T2 = ZERO
      T3 = ZERO
      D1 = DEXP(CP1)
      D2 = ONE
      DO 190 J = 1, N
         T1 = T1 + DFLOAT(N-J+1)*X(J)**2
         S1 = DEXP(X(J)/TEN)
         IF (J .EQ. 1) GO TO 180
         S3 = S1 + S2 - D2*(D1 + ONE)
         T2 = T2 + S3**2
         T3 = T3 + (S1 - ONE/D1)**2
  180    CONTINUE
         S2 = S1
         D2 = D1*D2
  190    CONTINUE
      F = AP*(T2 + T3) + BP*(T1**2 + (X(1) - CP2)**2)
      GO TO 390
C
C     BROWN BADLY SCALED FUNCTION.
C
  200 CONTINUE
      T1 = X(1) - C1PD6
      T2 = X(2) - C2PDM6
      T3 = X(1)*X(2) - TWO
      F = T1**2 + T2**2 + T3**2
      GO TO 390
C
C     BROWN AND DENNIS FUNCTION.
C
  210 CONTINUE
      F = ZERO
      DO 220 I = 1, 20
         D1 = DFLOAT(I)/FIVE
         D2 = DSIN(D1)
         T1 = X(1) + D1*X(2) - DEXP(D1)
         T2 = X(3) + D2*X(4) - DCOS(D1)
         T = T1**2 + T2**2
         F = F + T**2
  220    CONTINUE
      GO TO 390
C
C     GULF RESEARCH AND DEVELOPMENT FUNCTION.
C
  230 CONTINUE
      F = ZERO
      D1 = TWO/THREE
      DO 240 I = 1, 99
         ARG = DFLOAT(I)/C100
         R = DABS((-FIFTY*DLOG(ARG))**D1+C25-X(2))
         T1 = R**X(3)/X(1)
         T2 = DEXP(-T1)
         T = T2 - ARG
         F = F + T**2
  240    CONTINUE
      GO TO 390
C
C     TRIGONOMETRIC FUNCTION.
C
  250 CONTINUE
      S1 = ZERO
      DO 260 J = 1, N
         S1 = S1 + DCOS(X(J))
  260    CONTINUE
      F = ZERO
      DO 270 J = 1, N
         T = DFLOAT(N+J) - DSIN(X(J)) - S1 - DFLOAT(J)*DCOS(X(J))
         F = F + T**2
  270    CONTINUE
      GO TO 390
C
C     EXTENDED ROSENBROCK FUNCTION.
C
  280 CONTINUE
      F = ZERO
      DO 290 J = 1, N, 2
         T1 = ONE - X(J)
         T2 = TEN*(X(J+1) - X(J)**2)
         F = F + T1**2 + T2**2
  290    CONTINUE
      GO TO 390
C
C     EXTENDED POWELL FUNCTION.
C
  300 CONTINUE
      F = ZERO
      DO 310 J = 1, N, 4
         T = X(J) + TEN*X(J+1)
         T1 = X(J+2) - X(J+3)
         S1 = FIVE*T1
         T2 = X(J+1) - TWO*X(J+2)
         S2 = T2**3
         T3 = X(J) - X(J+3)
         S3 = TEN*T3**3
         F = F + T**2 + S1*T1 + S2*T2 + S3*T3
  310    CONTINUE
      GO TO 390
C
C     BEALE FUNCTION.
C
  320 CONTINUE
      S1 = ONE - X(2)
      T1 = C1P5 - X(1)*S1
      S2 = ONE - X(2)**2
      T2 = C2P25 - X(1)*S2
      S3 = ONE - X(2)**3
      T3 = C2P625 - X(1)*S3
      F = T1**2 + T2**2 + T3**2
      GO TO 390
C
C     WOOD FUNCTION.
C
  330 CONTINUE
      S1 = X(2) - X(1)**2
      S2 = ONE - X(1)
      S3 = X(2) - ONE
      T1 = X(4) - X(3)**2
      T2 = ONE - X(3)
      T3 = X(4) - ONE
      F = C100*S1**2 + S2**2 + C90*T1**2 + T2**2 + TEN*(S3 + T3)**2
     1    + (S3 - T3)**2/TEN
      GO TO 390
C
C     CHEBYQUAD FUNCTION.
C
  340 CONTINUE
      DO 350 I = 1, N
         FVEC(I) = ZERO
  350    CONTINUE
      DO 370 J = 1, N
         T1 = ONE
         T2 = TWO*X(J) - ONE
         T = TWO*T2
         DO 360 I = 1, N
            FVEC(I) = FVEC(I) + T2
            TH = T*T2 - T1
            T1 = T2
            T2 = TH
  360       CONTINUE
  370    CONTINUE
      F = ZERO
      D1 = ONE/DFLOAT(N)
      IEV = -1
      DO 380 I = 1, N
         T = D1*FVEC(I)
         IF (IEV .GT. 0) T = T + ONE/(DFLOAT(I)**2 - ONE)
         F = F + T**2
         IEV = -IEV
  380    CONTINUE
  390 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE OBJFCN.
C
      END
