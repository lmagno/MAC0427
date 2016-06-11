      SUBROUTINE INITPT(N,X,NPROB,FACTOR)                               00000010
      INTEGER N,NPROB
      DOUBLE PRECISION FACTOR
      DOUBLE PRECISION X(N)
C     **********
C
C     SUBROUTINE INITPT
C
C     THIS SUBROUTINE SPECIFIES THE STANDARD STARTING POINTS FOR THE
C     FUNCTIONS DEFINED BY SUBROUTINE OBJFCN. THE SUBROUTINE RETURNS
C     IN X A MULTIPLE (FACTOR) OF THE STANDARD STARTING POINT. FOR
C     THE SEVENTH FUNCTION THE STANDARD STARTING POINT IS ZERO, SO IN
C     THIS CASE, IF FACTOR IS NOT UNITY, THEN THE SUBROUTINE RETURNS
C     THE VECTOR  X(J) = FACTOR, J=1,...,N.
C
C     THE SUBROUTINE STATEMENT IS
C
C       SUBROUTINE INITPT(N,X,NPROB,FACTOR)
C
C     WHERE
C
C       N IS A POSITIVE INTEGER INPUT VARIABLE.
C
C       X IS AN OUTPUT ARRAY OF LENGTH N WHICH CONTAINS THE STANDARD
C         STARTING POINT FOR PROBLEM NPROB MULTIPLIED BY FACTOR.
C
C       NPROB IS A POSITIVE INTEGER INPUT VARIABLE WHICH DEFINES THE
C         NUMBER OF THE PROBLEM. NPROB MUST NOT EXCEED 18.
C
C       FACTOR IS AN INPUT VARIABLE WHICH SPECIFIES THE MULTIPLE OF
C         THE STANDARD STARTING POINT. IF FACTOR IS UNITY, NO
C         MULTIPLICATION IS PERFORMED.
C
C     MINPACK. VERSION OF JULY 1978.
C     BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE
C
C     **********
      INTEGER IVAR,J
      DOUBLE PRECISION C1,C2,C3,C4,FIVE,H,HALF,ONE,TEN,THREE,TWENTY,
     1                 TWNTF,TWO,ZERO
      DOUBLE PRECISION DFLOAT
      DATA ZERO,HALF,ONE,TWO,THREE,FIVE,TEN,TWENTY,TWNTF
     1     /0.0D0,0.5D0,1.0D0,2.0D0,3.0D0,5.0D0,1.0D1,2.0D1,2.5D1/
      DATA C1,C2,C3,C4 /4.0D-1,2.5D0,1.5D-1,1.2D0/
      DFLOAT(IVAR) = IVAR
C
C     SELECTION OF INITIAL POINT.
C
      GO TO (10,20,30,40,50,60,80,100,120,140,150,160,170,190,210,230,
     1       240,250), NPROB
C
C     HELICAL VALLEY FUNCTION.
C
   10 CONTINUE
      X(1) = -ONE
      X(2) = ZERO
      X(3) = ZERO
      GO TO 270
C
C     BIGGS EXP6 FUNCTION.
C
   20 CONTINUE
      X(1) = ONE
      X(2) = TWO
      X(3) = ONE
      X(4) = ONE
      X(5) = ONE
      X(6) = ONE
      GO TO 270
C
C     GAUSSIAN FUNCTION.
C
   30 CONTINUE
      X(1) = C1
      X(2) = ONE
      X(3) = ZERO
      GO TO 270
C
C     POWELL BADLY SCALED FUNCTION.
C
   40 CONTINUE
      X(1) = ZERO
      X(2) = ONE
      GO TO 270
C
C     BOX 3-DIMENSIONAL FUNCTION.
C
   50 CONTINUE
      X(1) = ZERO
      X(2) = TEN
      X(3) = TWENTY
      GO TO 270
C
C     VARIABLY DIMENSIONED FUNCTION.
C
   60 CONTINUE
      H = ONE/DFLOAT(N)
      DO 70 J = 1, N
         X(J) = ONE - DFLOAT(J)*H
   70    CONTINUE
      GO TO 270
C
C     WATSON FUNCTION.
C
   80 CONTINUE
      DO 90 J = 1, N
         X(J) = ZERO
   90    CONTINUE
      GO TO 270
C
C     PENALTY FUNCTION I.
C
  100 CONTINUE
      DO 110 J = 1, N
         X(J) = DFLOAT(J)
  110    CONTINUE
      GO TO 270
C
C     PENALTY FUNCTION II.
C
  120 CONTINUE
      DO 130 J = 1, N
         X(J) = HALF
  130    CONTINUE
      GO TO 270
C
C     BROWN BADLY SCALED FUNCTION.
C
  140 CONTINUE
      X(1) = ONE
      X(2) = ONE
      GO TO 270
C
C     BROWN AND DENNIS FUNCTION.
C
  150 CONTINUE
      X(1) = TWNTF
      X(2) = FIVE
      X(3) = -FIVE
      X(4) = -ONE
      GO TO 270
C
C     GULF RESEARCH AND DEVELOPMENT FUNCTION.
C
  160 CONTINUE
      X(1) = FIVE
      X(2) = C2
      X(3) = C3
      GO TO 270
C
C     TRIGONOMETRIC FUNCTION.
C
  170 CONTINUE
      H = ONE/DFLOAT(N)
      DO 180 J = 1, N
         X(J) = H
  180    CONTINUE
      GO TO 270
C
C     EXTENDED ROSENBROCK FUNCTION.
C
  190 CONTINUE
      DO 200 J = 1, N, 2
         X(J) = -C4
         X(J+1) = ONE
  200    CONTINUE
      GO TO 270
C
C     EXTENDED POWELL SINGULAR FUNCTION.
C
  210 CONTINUE
      DO 220 J = 1, N, 4
         X(J) = THREE
         X(J+1) = -ONE
         X(J+2) = ZERO
         X(J+3) = ONE
  220    CONTINUE
      GO TO 270
C
C     BEALE FUNCTION.
C
  230 CONTINUE
      X(1) = ONE
      X(2) = ONE
      GO TO 270
C
C     WOOD FUNCTION.
C
  240 CONTINUE
      X(1) = -THREE
      X(2) = -ONE
      X(3) = -THREE
      X(4) = -ONE
      GO TO 270
C
C     CHEBYQUAD FUNCTION.
C
  250 CONTINUE
      H = ONE/DFLOAT(N+1)
      DO 260 J = 1, N
         X(J) = DFLOAT(J)*H
  260    CONTINUE
  270 CONTINUE
C
C     COMPUTE MULTIPLE OF INITIAL POINT.
C
      IF (FACTOR .EQ. ONE) GO TO 320
      IF (NPROB .EQ. 7) GO TO 290
         DO 280 J = 1, N
            X(J) = FACTOR*X(J)
  280       CONTINUE
         GO TO 310
  290 CONTINUE
         DO 300 J = 1, N
            X(J) = FACTOR
  300       CONTINUE
  310 CONTINUE
  320 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE INITPT.
C
      END