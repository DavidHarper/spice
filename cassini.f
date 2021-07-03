      PROGRAM CASSINI
      IMPLICIT NONE
      INTEGER     STRLEN
      PARAMETER   ( STRLEN =  32 )
      CHARACTER*(STRLEN)    UTCBEG
      CHARACTER*(STRLEN)    UTCEND
      CHARACTER*(STRLEN)    ABCORR
      CHARACTER*(STRLEN)    FRAME
      CHARACTER*(STRLEN)    UTC
      CHARACTER*(4)         FORMAT
      DOUBLE PRECISION      DELTA
      DOUBLE PRECISION      ET
      DOUBLE PRECISION      ETBEG
      DOUBLE PRECISION      ETEND
      DOUBLE PRECISION      LT
      DOUBLE PRECISION      STATE  ( 6 )
      DOUBLE PRECISION DEARTH,DVENUS,DJUPTR
      DOUBLE PRECISION VHELIO, GMSUN
      DOUBLE PRECISION VNORM,TYEAR,TSTEP
      DOUBLE PRECISION ELMNTS(8)
      INTEGER I, N
      INTEGER PREC

      CALL FURNSH ( 'cassini.meta' )

      CALL PROMPT ( 'Enter the beginning UTC time: ', UTCBEG )

      CALL PROMPT ( 'Enter the ending UTC time: ', UTCEND )

      WRITE(6,*) 'Enter the time step in days: '
      READ(5,*) TSTEP

      CALL STR2ET ( UTCBEG, ETBEG )
      CALL STR2ET ( UTCEND, ETEND )

      ABCORR='NONE'
      FRAME='J2000'
      PREC=0
      FORMAT = 'ISOC'
      ET=ETBEG

      CALL BODVRD ( 'SUN', 'GM', 1, N, GMSUN )

      OPEN(UNIT=8,FILE='cassini.out',STATUS='UNKNOWN')

      DO WHILE (ET.LE.ETEND)
        CALL SPKEZR ( 'CASSINI', ET, FRAME, ABCORR, 'EARTH', STATE, LT )
        DEARTH=VNORM(STATE)
        CALL SPKEZR ( 'CASSINI', ET, FRAME, ABCORR, 'VENUS', STATE, LT )
        DVENUS=VNORM(STATE)
        CALL SPKEZR ( 'CASSINI', ET, FRAME, ABCORR, 'JUPITER', STATE,LT)
        DJUPTR=VNORM(STATE)
        CALL SPKEZR ( 'CASSINI', ET, FRAME, ABCORR, 'SUN', STATE, LT)
        VHELIO=VNORM(STATE(4))
        CALL OSCELT(STATE, ET, GMSUN, ELMNTS)

        TYEAR=2000.0D0+ET/(365.25D0*86400.0D0)
        CALL ET2UTC ( ET, FORMAT, PREC, UTC )
        WRITE(8,1000) TYEAR,DEARTH,DVENUS,DJUPTR,VHELIO,ELMNTS(2),UTC
        ET=ET+86400.0D0*TSTEP
      END DO

 1000 FORMAT(F11.6,3(1X,F13.3),1X,F13.6,1X,F10.6,1X,A32)

      CLOSE(UNIT=8,STATUS='KEEP')

      CALL BYEBYE( 'SUCCESS' )
      END
