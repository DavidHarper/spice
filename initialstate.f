      PROGRAM ISTATE
      IMPLICIT NONE
      INTEGER STRLEN
      PARAMETER (STRLEN=32)
      CHARACTER*(STRLEN) UTC0
      CHARACTER*(STRLEN) ABCORR
      CHARACTER*(STRLEN) FRAME
      CHARACTER*7 BODIES(0:4)
      DOUBLE PRECISION GM,STATE(6),ET0,JD0,LT
      INTEGER I,N
      DATA BODIES/'SUN','VENUS','EARTH','JUPITER','CASSINI'/
C
      CALL FURNSH('cassini.meta')
C
      CALL PROMPT('Enter the beginning UTC time: ', UTC0)
      CALL STR2ET(UTC0, ET0)
C
      OPEN(UNIT=8,FILE='initialstate.out',STATUS='UNKNOWN')
C
      JD0=2451545.0D0+ET0/86400.0D0
      WRITE(8,1000) JD0
 1000 FORMAT(F15.6)
C
      DO I=0,3
        CALL BODVRD(BODIES(I) , 'GM', 1, N, GM )
        WRITE(8,1001) GM
      END DO
 1001 FORMAT(D20.12)
C
      ABCORR='NONE'
      FRAME='J2000'
C
      DO I=1,4
        CALL SPKEZR(BODIES(I), ET0, FRAME, ABCORR, 'SUN', STATE, LT)
        WRITE(8,1002) STATE
      END DO
 1002 FORMAT(6(1X,D20.12))
C
      CLOSE(UNIT=8,STATUS='KEEP')
C
      CALL BYEBYE( 'SUCCESS' )
      END
