      INTEGER FUNCTION I8SAVE(ISW,IVALUE,SET)
C
C  IF (ISW = 1) I8SAVE RETURNS THE CURRENT ERROR NUMBER AND
C               SETS IT TO IVALUE IF SET = .TRUE. .
C
C  IF (ISW = 2) I8SAVE RETURNS THE CURRENT RECOVERY SWITCH AND
C               SETS IT TO IVALUE IF SET = .TRUE. .
C
      LOGICAL SET
C
      INTEGER IPARAM(2)
      EQUIVALENCE (IPARAM(1),LERROR) , (IPARAM(2),LRECOV)
C
C  START EXECUTION ERROR FREE AND WITH RECOVERY TURNED OFF.
C
      DATA LERROR/0/ , LRECOV/2/
C
      I8SAVE=IPARAM(ISW)
      IF (SET) IPARAM(ISW)=IVALUE
C
      RETURN
C
      END
