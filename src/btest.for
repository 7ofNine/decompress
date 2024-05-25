      FUNCTION BTEST(IVAL,IBIT)
C***************************************************************************
C
C_TITLE BTEST - test for bit turned on or off
C
C_ARGS   
      LOGICAL   BTEST
      INTEGER*4 IVAL
      INTEGER*4 IBIT
C
C_DESCR This function emulates the VAX/VMS system function for extracting
C       bit information from an integer long word. The function tests
C       the bit value at location 'bit' in the longword 'ival'. 
C         btest = .true.  if the bit is turned on
C         btest = .false. if the bit is turned off
C
C       The routine should be used only by non VAX/VMS users. The least
C       significant bit is bit 0, the most significant bit is bit 31.
C       The routine uses the Fortran system modulus routine 'mod'.
C
C_HIST 01-Mar-88 Eric Eliason USGS, Flagstaff, Original verison
C*****************************************************************************
       INTEGER*4 IWORK
       BTEST = .FALSE.
       IF (IBIT.GT.31.OR.IBIT.LT.0) RETURN
       IWORK = MOD(IVAL/(2**IBIT),2)
       IF (IWORK.EQ.1) BTEST = .TRUE.
       RETURN
       END
