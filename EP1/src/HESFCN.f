      SUBROUTINE HESFCN (N,X,HESD,HESL,NPROB)
      INTEGER N,NPROB
      DOUBLE PRECISION X(N),HESD(N),HESL(*)
C     **********
C
C     SUBROUTINE HESFCN
C
C     THIS SUBROUTINE DEFINES THE HESSIAN MATRICES OF 18
C     NONLINEAR UNCONSTRAINED MINIMIZATION PROBLEMS.  THE PROBLEM
C     DIMENSIONS ARE AS DESCRIBED IN OBJFCN.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE HESFCN (N, X, HESD, HESL, NPROB)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       X IS AN INPUT ARRAY OF LENGTH N.
C
C       HESD IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE
C         DIAGONAL COMPONENTS OF THE HESSIAN MATRIX OF THE NPROB
C         OBJECTIVE FUNCTION EVALUATED AT X.
C
C       HESL IS AN OUTPUT ARRAY OF LENGTH N*(N-1)/2 WHICH CONTAINS
C         THE LOWER TRIANGULAR PART OF THE HESSIAN MATRIX OF THE
C         NPROB OBJECTIVE FUNCTION EVALUATED AT X.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM.  NPROB MUST NOT EXCEED 18.
C
C     SUBPROGRAMS CALLED
C
C       FORTRAN-SUPPLIED ... ABS, ATAN, COS, EXP, LOG, SIGN, SIN,
C                            SQRT
C
C       INTEGER INLINE FUNCTION IX GIVES THE LOCATION OF A HESSIAN
C       ELEMENT (I,J), I>J, IN HESL
C
C     VICTORIA Z. AVERBUKH, SAMUEL A. FIGUEROA, AND
C     TAMAR SCHLICK, 1993.
C     **********
      INTEGER I, J, K, M, II, JJ, IX, IVAR
      DOUBLE PRECISION ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX, EIGHT,
     1      NINE, TEN, FIFTY, CP0001, CP1, CP2, CP25, CP5, C1P5, C2P25,
     2      C2P625, C3P5, C12, C19P8, C25, C29, C50, C90, C100, C120,
     3      C180, C200, C200P2, C202, C220P2, C360, C400, C1000, C1080,
     4      C1200, C2000, C20000, C2E8, C4E8, AP, BP, PI
      PARAMETER (ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0,
     1      FOUR=4.0D0, FIVE=5.0D0, SIX=6.0D0, EIGHT=8.0D0, NINE=9.0D0,
     2      TEN=1.0D1, FIFTY=5.0D1, CP0001=1.0D-4, CP1=1.0D-1,
     3      CP2=2.0D-1, CP25=2.5D-1, CP5=5.0D-1, C1P5=1.5D0,
     4      C2P25=2.25D0, C2P625=2.625D0, C3P5=3.5D0, C12=1.2D1,
     5      C19P8=1.98D1, C25=2.5D1, C29=2.9D1, C50=5.0D1, C90=9.0D1,
     6      C100=1.0D2, C120=1.2D2, C180=1.8D2, C200=2.0D2,
     7      C200P2=2.002D2, C202=2.02D2, C220P2=2.202D2, C360=3.6D2,
     8      C400=4.0D2, C1000=1.0D3, C1080=1.08D3, C1200=1.2D3,
     9      C2000=2.0D3, C20000=2.0D4, C2E8=2.0D8, C4E8=4.0D8,
     1      AP=1.0D-5, BP=ONE, PI=3.141592653589793D0)
      DOUBLE PRECISION ARG, D1, D2, D3, LOGR, P1, P2, PIARG, PIARG2,
     1      R, R3INV, S1, S2, S3, S1S2, S1S3, S2S3, SS1, SS2,
     2      T, T1, T2, T3, TH, TT, TT1, TT2, TTH
      DOUBLE PRECISION FVEC(50), GVEC(50), Y(15)
      LOGICAL IEV
      DOUBLE PRECISION DFLOAT
      IX(II,JJ)=(II-1)*(II-2)/2+JJ

      DFLOAT(IVAR) = IVAR
      DATA Y /9.0D-4, 4.4D-3, 1.75D-2, 5.4D-2, 1.295D-1, 2.42D-1,
     1      3.521D-1, 3.989D-1, 3.521D-1, 2.42D-1, 1.295D-1, 5.4D-2,
     2      1.75D-2, 4.4D-3, 9.0D-4/

C
C     HESSIAN ROUTINE SELECTOR.
C
      GO TO (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
     1      1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800), NPROB
C
C     HELICAL VALLEY FUNCTION.
C
  100 CONTINUE
C
      IF (X(1) .EQ. ZERO) THEN
         TH = SIGN(CP25,X(2))
      ELSE
         TH = ATAN(X(2)/X(1)) / (TWO*PI)
         IF (X(1) .LT. ZERO) TH = TH + CP5
      END IF
      ARG = X(1)**2 + X(2)**2
      PIARG = PI * ARG
      PIARG2 = PIARG * ARG
      R3INV = ONE / SQRT(ARG)**3
      T = X(3) - TEN*TH
      S1 = FIVE*T / PIARG
      P1 = C2000*X(1)*X(2)*T / PIARG2
      P2 = (FIVE/PIARG)**2
      HESD(1) = C200 - C200*(R3INV-P2)*X(2)**2 - P1
      HESD(2) = C200 - C200*(R3INV-P2)*X(1)**2 + P1
      HESD(3) = C202
      HESL(1) = C200*X(1)*X(2)*R3INV +
     1      C1000/PIARG2 * ( T*(X(1)**2-X(2)**2) - FIVE*X(1)*X(2)/PI )
      HESL(2) =  C1000*X(2) / PIARG
      HESL(3) = -C1000*X(1) / PIARG
      RETURN
C
C     BIGGS EXP6 FUNCTION.
C
  200 CONTINUE
      DO 210 I = 1, 6
         HESD(I) = ZERO
  210 CONTINUE
      DO 220 I = 1, 15
         HESL(I) = ZERO
  220 CONTINUE
      DO 230 I = 1, 13
         D1 = DFLOAT(I)/TEN
         D2 = EXP(-D1) - FIVE*EXP(-TEN*D1) + THREE*EXP(-FOUR*D1)
         S1 = EXP(-D1*X(1))
         S2 = EXP(-D1*X(2))
         S3 = EXP(-D1*X(5))
         T = X(3)*S1 - X(4)*S2 + X(6)*S3 - D2
         D2 = D1**2
         S1S2 = S1 * S2
         S1S3 = S1 * S3
         S2S3 = S2 * S3
         HESD(1) = HESD(1) + D2*S1*(T+X(3)*S1)
         HESD(2) = HESD(2) - D2*S2*(T-X(4)*S2)
         HESD(3) = HESD(3) + S1**2
         HESD(4) = HESD(4) + S2**2
         HESD(5) = HESD(5) + D2*S3*(T+X(6)*S3)
         HESD(6) = HESD(6) + S3**2
         HESL(1) = HESL(1) - D2*S1S2
         HESL(2) = HESL(2) - D1*S1*(T+X(3)*S1)
         HESL(3) = HESL(3) + D1*S1S2
         HESL(4) = HESL(4) + D1*S1S2
         HESL(5) = HESL(5) + D1*S2*(T-X(4)*S2)
         HESL(6) = HESL(6) - S1S2
         HESL(7) = HESL(7) + D2*S1S3
         HESL(8) = HESL(8) - D2*S2S3
         HESL(9) = HESL(9) - D1*S1S3
         HESL(10) = HESL(10) + D1*S2S3
         HESL(11) = HESL(11) - D1*S1S3
         HESL(12) = HESL(12) + D1*S2S3
         HESL(13) = HESL(13) + S1S3
         HESL(14) = HESL(14) - S2S3
         HESL(15) = HESL(15) - D1*S3*(T+X(6)*S3)
  230 CONTINUE
      HESD(1) = X(3)*HESD(1)
      HESD(2) = X(4)*HESD(2)
      HESD(5) = X(6)*HESD(5)
      HESL(1) = X(3)*X(4)*HESL(1)
      HESL(3) = X(4)*HESL(3)
      HESL(4) = X(3)*HESL(4)
      HESL(7) = X(3)*X(6)*HESL(7)
      HESL(8) = X(4)*X(6)*HESL(8)
      HESL(9) = X(6)*HESL(9)
      HESL(10) = X(6)*HESL(10)
      HESL(11) = X(3)*HESL(11)
      HESL(12) = X(4)*HESL(12)
      DO 240 I = 1, 6
         HESD(I) = TWO*HESD(I)
  240 CONTINUE
      DO 250 I = 1, 15
         HESL(I) = TWO*HESL(I)
  250 CONTINUE
      RETURN
C
C     GAUSSIAN FUNCTION.
C
  300 CONTINUE
      HESD(1) = ZERO
      HESD(2) = ZERO
      HESD(3) = ZERO
      HESL(1) = ZERO
      HESL(2) = ZERO
      HESL(3) = ZERO
      DO 310 I = 1, 15
         D1 = CP5*DFLOAT(I-1)
         D2 = C3P5 - D1 - X(3)
         ARG = CP5*X(2)*D2**2
         R = EXP(-ARG)
         T = X(1)*R - Y(I)
         T1 = TWO*X(1)*R - Y(I)
         HESD(1) = HESD(1) + R**2
         HESD(2) = HESD(2) + R*T1*D2**4
         HESD(3) = HESD(3) + R*(X(2)*T1*D2**2-T)
         HESL(1) = HESL(1) - R*T1*D2**2
         HESL(2) = HESL(2) + D2*R*T1
         HESL(3) = HESL(3) + D2*R*(T-ARG*T1)
  310 CONTINUE
      HESD(1) = TWO*HESD(1)
      HESD(2) = CP5*X(1)*HESD(2)
      HESD(3) = TWO*X(1)*X(2)*HESD(3)
      HESL(2) = TWO*X(2)*HESL(2)
      HESL(3) = TWO*X(1)*HESL(3)
      RETURN
C
C     POWELL BADLY SCALED FUNCTION.
C
  400 CONTINUE
      S1 = EXP(-X(1))
      S2 = EXP(-X(2))
      T2 = S1 + S2 - ONE - CP0001
      HESD(1) = C2E8*X(2)**2 + TWO*S1*(S1+T2)
      HESD(2) = C2E8*X(1)**2 + TWO*S2*(S2+T2)
      HESL(1) = C4E8*X(1)*X(2) + TWO*S1*S2 - C20000
      RETURN
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
  500 CONTINUE
      HESD(1) = ZERO
      HESD(2) = ZERO
      HESD(3) = ZERO
      HESL(1) = ZERO
      HESL(2) = ZERO
      HESL(3) = ZERO
      DO 510 I = 1, 10
         D1 = DFLOAT(I)
         D2 = D1/TEN
         S1 = EXP(-D2*X(1))
         S2 = EXP(-D2*X(2))
         S3 = EXP(-D2) - EXP(-D1)
         T = S1 - S2 - S3*X(3)
         TH = T*D2**2
         HESD(1) = HESD(1) + TH*S1 + (D2*S1)**2
         HESD(2) = HESD(2) - TH*S2 + (D2*S2)**2
         HESD(3) = HESD(3) + S3**2
         HESL(1) = HESL(1) - S1*S2*D2**2
         HESL(2) = HESL(2) + D2*S1*S3
         HESL(3) = HESL(3) - D2*S2*S3
  510 CONTINUE
      HESD(1) = TWO*HESD(1)
      HESD(2) = TWO*HESD(2)
      HESD(3) = TWO*HESD(3)
      HESL(1) = TWO*HESL(1)
      HESL(2) = TWO*HESL(2)
      HESL(3) = TWO*HESL(3)
      RETURN
C
C     VARIABLY DIMENSIONED FUNCTION.
C
  600 CONTINUE
      T1 = ZERO
      DO 610 J = 1, N
         T1 = T1 + DFLOAT(J)*(X(J)-ONE)
  610 CONTINUE
      T = ONE + SIX*T1**2
      M = 0
      DO 630 J = 1, N
         HESD(J) = TWO + TWO*T*DFLOAT(J)**2
         DO 620 K = 1, J-1
            M = M + 1
            HESL(M) = TWO*T*DFLOAT(J*K)
  620    CONTINUE
  630 CONTINUE
      RETURN
C
C     WATSON FUNCTION.
C
  700 CONTINUE
      DO 710 J = 1, N
         HESD(J) = ZERO
  710 CONTINUE
      DO 720 J = 1, N*(N-1)/2
         HESL(J) = ZERO
  720 CONTINUE
      DO 760 I = 1, 29
         D1 = DFLOAT(I)/C29
         D2 = ONE
         S1 = ZERO
         S2 = X(1)
         DO 730 J = 2, N
            S1 = S1 + DFLOAT(J-1)*D2*X(J)
            D2 = D1*D2
            S2 = S2 + D2*X(J)
  730    CONTINUE
         T = TWO * (S1-S2**2-ONE) * D1**2
         S3 = TWO*D1*S2
         D2 = ONE/D1
         M = 0
         DO 750 J = 1, N
            T1 = DFLOAT(J-1) - S3
            HESD(J) = HESD(J) + (T1**2-T)*D2**2
            D3 = ONE/D1
            DO 740 K = 1, J-1
               M = M + 1
               HESL(M) = HESL(M) + (T1*(DFLOAT(K-1)-S3) - T) * D2*D3
               D3 = D1*D3
  740       CONTINUE
            D2 = D1*D2
  750    CONTINUE
  760 CONTINUE
      T3 = X(2) - X(1)**2 - ONE
      HESD(1) = HESD(1) + ONE - TWO*(T3-TWO*X(1)**2)
      HESD(2) = HESD(2) + ONE
      HESL(1) = HESL(1) - TWO*X(1)
      DO 770 J = 1, N
         HESD(J) = TWO * HESD(J)
  770 CONTINUE
      DO 780 J = 1, N*(N-1)/2
         HESL(J) = TWO * HESL(J)
  780 CONTINUE
      RETURN
C
C     PENALTY FUNCTION I.
C
  800 CONTINUE
      T1 = -CP25
      DO 810 J = 1, N
         T1 = T1 + X(J)**2
  810 CONTINUE
      D1 = TWO*AP
      TH = FOUR*BP*T1
      M = 0
      DO 830 J = 1, N
         HESD(J) = D1 + TH + EIGHT*X(J)**2
         DO 820 K = 1, J-1
            M = M + 1
            HESL(M) = EIGHT*X(J)*X(K)
  820    CONTINUE
  830 CONTINUE
      RETURN
C
C     PENALTY FUNCTION II.
C
  900 CONTINUE
      T1 = -ONE
      DO 910 J = 1, N
         T1 = T1 + DFLOAT(N-J+1)*X(J)**2
  910 CONTINUE
      D1 = EXP(CP1)
      D2 = ONE
      TH = FOUR*BP*T1
      M = 0
      DO 930 J = 1, N
         HESD(J) = EIGHT*BP*(DFLOAT(N-J+1)*X(J))**2 + DFLOAT(N-J+1)*TH
         S1 = EXP(X(J)/TEN)
         IF (J .GT. 1) THEN
            S3 = S1 + S2 - D2*(D1 + ONE)
            HESD(J) = HESD(J) + AP*S1*(S3 + S1 - ONE/D1 + TWO*S1)/C50
            HESD(J-1) = HESD(J-1) + AP*S2*(S2+S3)/C50
            DO 920 K = 1, J-1
               M = M + 1
               T1 = EXP(DFLOAT(K)/TEN)
               HESL(M) = EIGHT*DFLOAT(N-J+1)*DFLOAT(N-K+1)*X(J)*X(K)
  920       CONTINUE
            HESL(M) = HESL(M) + AP*S1*S2/C50
         END IF
         S2 = S1
         D2 = D1*D2
  930 CONTINUE
      HESD(1) = HESD(1) + TWO*BP
      RETURN
C
C     BROWN BADLY SCALED FUNCTION.
C
 1000 CONTINUE
      HESD(1) = TWO + TWO*X(2)**2
      HESD(2) = TWO + TWO*X(1)**2
      HESL(1) = FOUR*X(1)*X(2) - FOUR
      RETURN
C
C     BROWN AND DENNIS FUNCTION.
C
 1100 CONTINUE
      DO 1110 I = 1, 4
         HESD(I) = ZERO
 1110 CONTINUE
      DO 1120 I = 1, 6
         HESL(I) = ZERO
 1120 CONTINUE
      DO 1130 I = 1, 20
         D1 = DFLOAT(I)/FIVE
         D2 = SIN(D1)
         T1 = X(1) + D1*X(2) - EXP(D1)
         T2 = X(3) + D2*X(4) - COS(D1)
         T = EIGHT * T1 * T2
         S1 = C12*T1**2 + FOUR*T2**2
         S2 = C12*T2**2 + FOUR*T1**2
         HESD(1) = HESD(1) + S1
         HESD(2) = HESD(2) + S1*D1**2
         HESD(3) = HESD(3) + S2
         HESD(4) = HESD(4) + S2*D2**2
         HESL(1) = HESL(1) + S1*D1
         HESL(2) = HESL(2) + T
         HESL(4) = HESL(4) + T*D2
         HESL(3) = HESL(3) + T*D1
         HESL(5) = HESL(5) + T*D1*D2
         HESL(6) = HESL(6) + S2*D2
 1130 CONTINUE
      RETURN
C
C     GULF RESEARCH AND DEVELOPMENT FUNCTION.
C
 1200 CONTINUE
      DO 1210 I = 1, 3
         HESD(I) = ZERO
         HESL(I) = ZERO
 1210 CONTINUE
      D1 = TWO/THREE
      DO 1220 I = 1, 99
         ARG = DFLOAT(I)/C100
         R = (-FIFTY*LOG(ARG))**D1+C25-X(2)
         T1 = ABS(R)**X(3)/X(1)
         T2 = EXP(-T1)
         T3 = T1 * T2 * (T1*T2+(T1-ONE)*(T2-ARG))
         T = T1 * T2 * (T2-ARG)
         LOGR = LOG(ABS(R))
         HESD(1) = HESD(1) + T3 - T
         HESD(2) = HESD(2) + (T+X(3)*T3)/R**2
         HESD(3) = HESD(3) + T3*LOGR**2
         HESL(1) = HESL(1) + T3/R
         HESL(2) = HESL(2) - T3*LOGR
         HESL(3) = HESL(3) + (T-X(3)*T3*LOGR)/R
 1220 CONTINUE
      HESD(1) = HESD(1) / X(1)**2
      HESD(2) = HESD(2) * X(3)
      HESL(1) = HESL(1) * X(3)/X(1)
      HESL(2) = HESL(2) / X(1)
      DO 1230 I = 1, 3
         HESD(I) = TWO * HESD(I)
         HESL(I) = TWO * HESL(I)
 1230 CONTINUE
      RETURN
C
C     TRIGONOMETRIC FUNCTION.
C
 1300 CONTINUE
      S1 = ZERO
      DO 1310 J = 1, N
         HESD(J) = SIN(X(J))
         S1 = S1 + COS(X(J))
 1310 CONTINUE
      S2 = ZERO
      M = 0
      DO 1330 J = 1, N
         TH = COS(X(J))
         T = DFLOAT(N+J) - HESD(J) - S1 - DFLOAT(J)*TH
         S2 = S2 + T
         DO 1320 K = 1, J-1
            M = M + 1
            HESL(M) = SIN(X(K))*(DFLOAT(N+J+K)*HESD(J)-TH) -
     *            HESD(J)*COS(X(K))
            HESL(M) = TWO*HESL(M)
 1320    CONTINUE
         HESD(J) = DFLOAT(J*(J+2)+N)*HESD(J)**2 +
     *         TH*(TH-DFLOAT(2*J+2)*HESD(J)) + T*(DFLOAT(J)*TH+HESD(J))
 1330 CONTINUE
      DO 1340 J = 1, N
         HESD(J) = TWO*(HESD(J) + COS(X(J))*S2)
 1340 CONTINUE
      RETURN
C
C     EXTENDED ROSENBROCK FUNCTION.
C
 1400 CONTINUE
      DO 1410 J = 1, N*(N-1)/2
         HESL(J) = ZERO
 1410 CONTINUE
      DO 1420 J = 1, N, 2
         HESD(J+1) = C200
         HESD(J) = C1200*X(J)**2 - C400*X(J+1) + TWO
         HESL(IX(J+1,J)) = -C400*X(J)
 1420 CONTINUE
      RETURN
C
C     EXTENDED POWELL FUNCTION.
C
 1500 CONTINUE
      DO 1510 J = 1, N*(N-1)/2
         HESL(J) = ZERO
 1510 CONTINUE
      DO 1520 J = 1, N, 4
         T2 = X(J+1) - TWO*X(J+2)
         T3 = X(J) - X(J+3)
         S1 = C12 * T2**2
         S2 = C120 * T3**2
         HESD(J) = TWO + S2
         HESD(J+1) = C200 + S1
         HESD(J+2) = TEN + FOUR*S1
         HESD(J+3) = TEN + S2
         HESL(IX(J+1,J)) = TWO*TEN
         HESL(IX(J+2,J)) = ZERO
         HESL(IX(J+2,J+1)) = -TWO*S1
         HESL(IX(J+3,J)) = -S2
         HESL(IX(J+3,J+1)) = ZERO
         HESL(IX(J+3,J+2)) = -TEN
 1520 CONTINUE
      RETURN
C
C     BEALE FUNCTION.
C
 1600 CONTINUE
      S1 = ONE - X(2)
      T1 = C1P5 - X(1)*S1
      S2 = ONE - X(2)**2
      T2 = C2P25 - X(1)*S2
      S3 = ONE - X(2)**3
      T3 = C2P625 - X(1)*S3
      HESD(1) = TWO * (S1**2 + S2**2 + S3**2)
      HESD(2) = TWO*X(1) * (X(1) + TWO*T2 + FOUR*X(1)*X(2)**2 +
     1      SIX*X(2)*T3 + NINE*X(1)*X(2)**4)
      HESL(1) = TWO*(T1-X(1)*S1) + FOUR*X(2)*(T2-X(1)*S2) +
     2      SIX*(T3-X(1)*S3)*X(2)**2
      RETURN
C
C     WOOD FUNCTION.
C
 1700 CONTINUE
      HESD(1) = C1200*X(1)**2 - C400*X(2) + TWO
      HESD(2) = C220P2
      HESD(3) = C1080*X(3)**2 - C360*X(4) + TWO
      HESD(4) = C200P2
      HESL(1) = -C400*X(1)
      HESL(2) = ZERO
      HESL(3) = ZERO
      HESL(4) = ZERO
      HESL(5) = C19P8
      HESL(6) = -C360*X(3)
      RETURN
C
C     CHEBYQUAD FUNCTION.
C
 1800 CONTINUE
      DO 1810 I = 1, N
         FVEC(I) = ZERO
 1810 CONTINUE
      DO 1830 J = 1, N
         T1 = ONE
         T2 = TWO*X(J) - ONE
         T = TWO*T2
         DO 1820 I = 1, N
            FVEC(I) = FVEC(I) + T2
            TH = T*T2 - T1
            T1 = T2
            T2 = TH
 1820    CONTINUE
 1830 CONTINUE
      D1 = ONE/FLOAT(N)
      IEV = .FALSE.
      DO 1840 I = 1, N
         FVEC(I) = D1*FVEC(I)
         IF (IEV) FVEC(I) = FVEC(I) + ONE/(DFLOAT(I)**2 - ONE)
         IEV = .NOT. IEV
 1840 CONTINUE
      D2 = TWO*D1
      M = 0
      DO 1880 J = 1, N
         HESD(J) = FOUR*D1
         T1 = ONE
         T2 = TWO*X(J) - ONE
         T = TWO*T2
         S1 = ZERO
         S2 = TWO
         P1 = ZERO
         P2 = ZERO
         GVEC(1) = S2
         DO 1850 I = 2, N
            TH = FOUR*T2 + T*S2 - S1
            S1 = S2
            S2 = TH
            TH = T*T2 - T1
            T1 = T2
            T2 = TH
            TH = EIGHT*S1 + T*P2 - P1
            P1 = P2
            P2 = TH
            GVEC(I) = S2
            HESD(J) = HESD(J) + FVEC(I)*TH + D1*S2**2
 1850    CONTINUE
         HESD(J) = D2*HESD(J)
         DO 1870 K = 1, J-1
            M = M + 1
            HESL(M) = ZERO
            TT1 = ONE
            TT2 = TWO*X(K) - ONE
            TT = TWO*TT2
            SS1 = ZERO
            SS2 = TWO
            DO 1860 I = 1, N
               HESL(M) = HESL(M) + SS2*GVEC(I)
               TTH = FOUR*TT2 + TT*SS2 - SS1
               SS1 = SS2
               SS2 = TTH
               TTH = TT*TT2 - TT1
               TT1 = TT2
               TT2 = TTH
 1860       CONTINUE
            HESL(M) = D2*D1*HESL(M)
 1870    CONTINUE
 1880 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE HESFCN.
C
      END
