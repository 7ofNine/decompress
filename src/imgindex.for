C  INDXTEST:  Test Fortran Formatted Read of Voyager Image Index Table

      INTEGER*4    UNIT /10/
      CHARACTER*9  SCNAME
      CHARACTER*17 MSNPHSNM
      CHARACTER*8  TARGETNM
      CHARACTER*10 IMAGEID
      REAL*4       IMAGENUM
      CHARACTER*20 IMAGETM,ETHRCDTM
      CHARACTER*19 INSTRNM
      CHARACTER*7  SCNMODID,SHTMODID,GNMODID,EDTMODID,FILTERNM
      INTEGER*4    FLTRNUM
      REAL*4       EXPOSRDU
      CHARACTER*80 NOTE
      CHARACTER*8  SMPBTMSK
      CHARACTER*6  DATANMTP
      CHARACTER*8  IMGVOLNM,BRSVOLNM
      CHARACTER*31 IMGFILNM
      CHARACTER*38 BRSFILNM

      OPEN (10,FILE='IMGINDEX.TAB',STATUS='OLD')
      IREC=0

 10   READ (10,100,END=80) SCNAME,MSNPHSNM,TARGETNM,
     1 IMAGEID,IMAGENUM,IMAGETM,ETHRCDTM,INSTRNM,
     2 SCNMODID,SHTMODID,GNMODID,EDTMODID,FILTERNM,
     3 FLTRNUM,EXPOSRDU,NOTE,SMPBTMSK,DATANMTP,
     4 IMGVOLNM,IMGFILNM,BRSVOLNM,BRSFILNM
      IREC=IREC+1

 100  FORMAT (1X,A9,3X,A17,3X,A8,
     1 3X,A10,2X,F8.2,2X,A20,3X,A20,3X,A19,
     2 3X,A7,3X,A7,3X,A7,3X,A7,3X,A7,
     3 2X,I4,1X,F7.4,2X,A80,3X,A8,3X,A6,
     4 3X,A8,3X,A31,3X,A8,3X,A38)

      WRITE (*,200)IREC,SCNAME,MSNPHSNM,TARGETNM,
     1 IMAGEID,IMAGENUM,IMAGETM,ETHRCDTM,INSTRNM,
     2 SCNMODID,SHTMODID,GNMODID,EDTMODID,FILTERNM,
     3 FLTRNUM,EXPOSRDU,NOTE,SMPBTMSK,DATANMTP,
     4 IMGVOLNM,IMGFILNM,BRSVOLNM,BRSFILNM

 200  FORMAT ('0RECORD',I5,1X,A9,3X,A16,3X,A8/
     1 3X,A10,2X,F8.2,2X,A20,3X,A20/
     2 3X,A19,3X,A7,3X,A7,3X,A7,3X,A7,3X,A7/
     3 2X,I4,1X,F7.4/
     4 2X,A80/
     5 3X,A8,3X,A6/
     6 3X,A8,3X,A31/
     6 3X,A8,3X,A38)

      GO TO 10

 80   CLOSE (10)
      STOP
      END
