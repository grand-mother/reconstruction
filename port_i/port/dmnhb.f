      SUBROUTINE  DMNHB(N, D, X, B, CALCF, CALCGH, IV, LIV, LV, V,
     1                  UIPARM, URPARM, UFPARM)
C
C  ***  MINIMIZE GENERAL SIMPLY BOUNDED OBJECTIVE FUNCTION USING   ***
C  ***  (ANALYTIC) GRADIENT AND HESSIAN PROVIDED BY THE CALLER.    ***
C
      INTEGER LIV, LV, N
C/6S
C     INTEGER IV(LIV), UIPARM(1)
C     DOUBLE PRECISION B(2,N), D(N), X(N), V(LV), URPARM(1)
C/7S
      INTEGER IV(LIV), UIPARM(*)
      DOUBLE PRECISION B(2,N), D(N), X(N), V(LV), URPARM(*)
C/
C     DIMENSION IV(59 + 3*N), V(78 + N*(N+15))
      EXTERNAL CALCF, CALCGH, UFPARM
C
C------------------------------  DISCUSSION  ---------------------------
C
C        THIS ROUTINE IS LIKE  DMNGB, EXCEPT THAT THE SUBROUTINE PARA-
C     METER CALCG OF  DMNGB (WHICH COMPUTES THE GRADIENT OF THE OBJEC-
C     TIVE FUNCTION) IS REPLACED BY THE SUBROUTINE PARAMETER CALCGH,
C     WHICH COMPUTES BOTH THE GRADIENT AND (LOWER TRIANGLE OF THE)
C     HESSIAN OF THE OBJECTIVE FUNCTION.  THE CALLING SEQUENCE IS...
C             CALL CALCGH(N, X, NF, G, H, UIPARM, URPARM, UFPARM)
C     PARAMETERS N, X, NF, G, UIPARM, URPARM, AND UFPARM ARE THE SAME
C     AS FOR  DMNGB, WHILE H IS AN ARRAY OF LENGTH N*(N+1)/2 IN WHICH
C     CALCGH MUST STORE THE LOWER TRIANGLE OF THE HESSIAN AT X.  START-
C     ING AT H(1), CALCGH MUST STORE THE HESSIAN ENTRIES IN THE ORDER
C     (1,1), (2,1), (2,2), (3,1), (3,2), (3,3), ...
C        THE VALUE PRINTED (BY DITSUM) IN THE COLUMN LABELLED STPPAR
C     IS THE LEVENBERG-MARQUARDT USED IN COMPUTING THE CURRENT STEP.
C     ZERO MEANS A FULL NEWTON STEP.  IF THE SPECIAL CASE DESCRIBED IN
C     REF. 1 IS DETECTED, THEN STPPAR IS NEGATED.  THE VALUE PRINTED
C     IN THE COLUMN LABELLED NPRELDF IS ZERO IF THE CURRENT HESSIAN
C     IS NOT POSITIVE DEFINITE.
C        IT SOMETIMES PROVES WORTHWHILE TO LET D BE DETERMINED FROM THE
C     DIAGONAL OF THE HESSIAN MATRIX BY SETTING IV(DTYPE) = 1 AND
C     V(DINIT) = 0.  THE FOLLOWING IV AND V COMPONENTS ARE RELEVANT...
C
C IV(DTOL)..... IV(59) GIVES THE STARTING SUBSCRIPT IN V OF THE DTOL
C             ARRAY USED WHEN D IS UPDATED.  (IV(DTOL) CAN BE
C             INITIALIZED BY CALLING  DMNHB WITH IV(1) = 13.)
C IV(DTYPE).... IV(16) TELLS HOW THE SCALE VECTOR D SHOULD BE CHOSEN.
C             IV(DTYPE) .LE. 0 MEANS THAT D SHOULD NOT BE UPDATED, AND
C             IV(DTYPE) .GE. 1 MEANS THAT D SHOULD BE UPDATED AS
C             DESCRIBED BELOW WITH V(DFAC).  DEFAULT = 0.
C V(DFAC)..... V(41) AND THE DTOL AND D0 ARRAYS (SEE V(DTINIT) AND
C             V(D0INIT)) ARE USED IN UPDATING THE SCALE VECTOR D WHEN
C             IV(DTYPE) .GT. 0.  (D IS INITIALIZED ACCORDING TO
C             V(DINIT), DESCRIBED IN  DMNG.)  LET
C                  D1(I) = MAX(SQRT(ABS(H(I,I))), V(DFAC)*D(I)),
C             WHERE H(I,I) IS THE I-TH DIAGONAL ELEMENT OF THE CURRENT
C             HESSIAN.  IF IV(DTYPE) = 1, THEN D(I) IS SET TO D1(I)
C             UNLESS D1(I) .LT. DTOL(I), IN WHICH CASE D(I) IS SET TO
C                  MAX(D0(I), DTOL(I)).
C             IF IV(DTYPE) .GE. 2, THEN D IS UPDATED DURING THE FIRST
C             ITERATION AS FOR IV(DTYPE) = 1 (AFTER ANY INITIALIZATION
C             DUE TO V(DINIT)) AND IS LEFT UNCHANGED THEREAFTER.
C             DEFAULT = 0.6.
C V(DTINIT)... V(39), IF POSITIVE, IS THE VALUE TO WHICH ALL COMPONENTS
C             OF THE DTOL ARRAY (SEE V(DFAC)) ARE INITIALIZED.  IF
C             V(DTINIT) = 0, THEN IT IS ASSUMED THAT THE CALLER HAS
C             STORED DTOL IN V STARTING AT V(IV(DTOL)).
C             DEFAULT = 10**-6.
C V(D0INIT)... V(40), IF POSITIVE, IS THE VALUE TO WHICH ALL COMPONENTS
C             OF THE D0 VECTOR (SEE V(DFAC)) ARE INITIALIZED.  IF
C             V(DFAC) = 0, THEN IT IS ASSUMED THAT THE CALLER HAS
C             STORED D0 IN V STARTING AT V(IV(DTOL)+N).  DEFAULT = 1.0.
C
C  ***  REFERENCE  ***
C
C 1. GAY, D.M. (1981), COMPUTING OPTIMAL LOCALLY CONSTRAINED STEPS,
C         SIAM J. SCI. STATIST. COMPUT. 2, PP. 186-197.
C.
C  ***  GENERAL  ***
C
C     CODED BY DAVID M. GAY (WINTER, SPRING 1983).
C
C----------------------------  DECLARATIONS  ---------------------------
C
      EXTERNAL DIVSET, DRMNHB
C
C DIVSET.... PROVIDES DEFAULT INPUT VALUES FOR IV AND V.
C DRMNHB... REVERSE-COMMUNICATION ROUTINE THAT DOES  DMNHB ALGORITHM.
C
      INTEGER G1, H1, IV1, LH, NF
      DOUBLE PRECISION F
C
C  ***  SUBSCRIPTS FOR IV   ***
C
      INTEGER G, H, NEXTV, NFCALL, NFGCAL, TOOBIG, VNEED
C
C/6
C     DATA NEXTV/47/, NFCALL/6/, NFGCAL/7/, G/28/, H/56/, TOOBIG/2/,
C    1     VNEED/4/
C/7
      PARAMETER (NEXTV=47, NFCALL=6, NFGCAL=7, G=28, H=56, TOOBIG=2,
     1           VNEED=4)
C/
C
C+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
C
      LH = N * (N + 1) / 2
      IF (IV(1) .EQ. 0) CALL DIVSET(2, IV, LIV, LV, V)
      IV1 = IV(1)
      IF (IV1 .EQ. 14) GO TO 10
      IF (IV1 .GT. 2 .AND. IV1 .LT. 12) GO TO 10
      IF (IV1 .EQ. 12) IV(1) = 13
      IF (IV(1) .EQ. 13) IV(VNEED) = IV(VNEED) + N*(N+3)/2
      CALL DRMNHB(B, D, F, V, V, IV, LH, LIV, LV, N, V, X)
      IF (IV(1) .NE. 14) GO TO 999
C
C  ***  STORAGE ALLOCATION
C
      IV(G) = IV(NEXTV)
      IV(H) = IV(G) + N
      IV(NEXTV) = IV(H) + N*(N+1)/2
      IF (IV1 .EQ. 13) GO TO 999
C
 10   G1 = IV(G)
      H1 = IV(H)
C
 20   CALL DRMNHB(B, D, F, V(G1), V(H1), IV, LH, LIV, LV, N, V, X)
      IF (IV(1) - 2) 30, 40, 999
C
 30   NF = IV(NFCALL)
      CALL CALCF(N, X, NF, F, UIPARM, URPARM, UFPARM)
      IF (NF .LE. 0) IV(TOOBIG) = 1
      GO TO 20
C
 40   NF = IV(NFGCAL)
      CALL CALCGH(N, X, NF, V(G1), V(H1), UIPARM, URPARM, UFPARM)
      IF (NF .LE. 0) IV(TOOBIG) = 1
      GO TO 20
C
 999  RETURN
C  ***  LAST CARD OF  DMNHB FOLLOWS  ***
      END
