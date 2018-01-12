PROGRAM CHK_DATA

!-------------------------------------------------------------------------------
! Simple program to dump the first few elements of the data array for each
! record of an ARL packed meteorological file. Used for diagnostic testing.
! Created: 23 Nov 1999 (RRD)
!          14 Dec 2000 (RRD) - fortran90 upgrade
!          18 Oct 2001 (RRD) - expanded grid domain
!          03 Jun 2008 (RRD) - embedded blanks
!-------------------------------------------------------------------------------
!!	 modified to dump the elements over the Auger South Array


  REAL,          ALLOCATABLE :: RDATA(:,:)
  CHARACTER(1),  ALLOCATABLE :: CPACK(:)

  CHARACTER(4)               :: KVAR, MODEL
  CHARACTER(50)              :: LABEL
  CHARACTER(80)              :: FDIR, FILE
  CHARACTER(3072)            :: HEADER
  LOGICAL                    :: FTEST


  character(3)	:: month_char
  character(2)	:: year, month_digi
  character(1)	:: week
  integer	:: i

!-------------------------------------------------------------------------------
  INTERFACE
  SUBROUTINE UNPACK(CPACK,RDATA,NX,NY,NEXP,VAR1)
  CHARACTER(1),INTENT(IN)  :: CPACK(:)
  REAL,        INTENT(OUT) :: RDATA(:,:)
  INTEGER,     INTENT(IN)  :: NX,NY,NEXP
  REAL,        INTENT(IN)  :: VAR1
  END SUBROUTINE
  END INTERFACE
!-------------------------------------------------------------------------------

  call getarg( 1, FILE )
  FDIR = ADJUSTL(FILE)

  KLEN=LEN_TRIM(FILE)
  INQUIRE(FILE=FILE(1:KLEN), EXIST=FTEST)
  IF(.NOT.FTEST)THEN
     WRITE(*,*)'Unable to find file: ',FILE(1:KLEN)
     STOP
  END IF

  ! open file to decode the standard label (50) plus the
  ! fixed portion (108) of the extended header
  OPEN(unit=11,FILE=FILE(1:KLEN),RECL=158,ACCESS='DIRECT',FORM='UNFORMATTED')

  ! decode the standard portion of the index record
  READ(11,REC=1)LABEL,HEADER(1:108)
  READ(LABEL,'(5I2,4X,A4)')IYR,IMO,IDA,IHR,IFC,KVAR
  WRITE(*,'(A,4I5)')'Opened file       : ',IYR,IMO,IDA,IHR

  IF(KVAR.NE.'INDX')THEN
     WRITE(*,*)'WARNING Old format meteo data grid'
     WRITE(*,*)LABEL
     WRITE(*,*)HEADER(1:108)
     STOP
  END IF

  ! decode extended portion of the header
  READ(HEADER(1:108),'(A4,I3,I2,12F7.0,3I3,I2,I4)',ERR=900)                    &
         MODEL,    ICX,       MN,                                              &
         POLE_LAT, POLE_LON,  REF_LAT,                                         &
         REF_LON,  SIZE,      ORIENT,                                          &
         TANG_LAT, SYNC_XP,   SYNC_YP,                                         &
         SYNC_LAT, SYNC_LON,  DUMMY,                                           &
         NX,       NY,        NZ,                                              &
         K_FLAG,   LENH

  ! close file and reopen with proper length
  CLOSE (11)
  NXY = NX*NY
  LEN = NXY+50
  OPEN(unit=11,FILE=FILE(1:KLEN),RECL=LEN,ACCESS='DIRECT',FORM='UNFORMATTED')

  ! print file diagnostic
  WRITE(*,'(A,4I5)')'Grid size and lrec: ',NX,NY,NXY,LEN
  WRITE(*,'(A,I5)') 'Header record size: ',LENH

  ! allocate array space
  ALLOCATE (RDATA(NX,NY), STAT=KRET)
  ALLOCATE (CPACK(NXY),   STAT=KRET)

  ! read entire file and print headers
  KREC=1
  100 READ(11,REC=KREC,ERR=800)LABEL,(CPACK(K),K=1,NXY)
      READ(LABEL,'(6I2,2X,A4,I4,2E14.7)',ERR=900) IY,IM,ID,IH,IF,KL,  &
                                             KVAR,NEXP,PREC,VAR1

  open(unit=12,file=FILE(1:KLEN)//'_magic1')
  ! start at level 1 (=1000hPa): year month day hour level
  IF(KVAR.EQ.'HGTS'.AND.KL.GE.1) WRITE(12,'(5I6$)')IY,IM,ID,IH,KL
  IF(KL.GE.1) THEN
    CALL UNPACK(CPACK,RDATA,NX,NY,NEXP,VAR1)
    IF(KVAR.EQ.'RELH') WRITE(12,'(A)') " "
    IF(KVAR.EQ.'VWND'.AND.KL.GE.22) THEN
      WRITE(12,'(F10.2)',advance="no") 0.
      WRITE(12,'(F10.2)') 0.
    END IF
  END IF
  IF(KVAR.EQ.'PRSS'.AND.KL.EQ.0) WRITE(12,'(5I6$)')IY,IM,ID,IH, KL
  IF((KVAR.NE.'INDX').AND.KL.EQ.0) THEN
    CALL UNPACK(CPACK,RDATA,NX,NY,NEXP,VAR1)
    IF(KVAR.EQ.'HCLD') WRITE(12,'(A)') " "
  END IF

  KREC=KREC+1
  GO TO 100

  900 WRITE(*,*)'ERROR: decoding header'
      WRITE(*,*)LABEL
      WRITE(*,*)HEADER(1:108)
      CLOSE (11)

  800 close(11)


END PROGRAM chk_data

!-------------------------------------------------------------------------------

SUBROUTINE UNPACK(CPACK,RDATA,NX,NY,NEXP,VAR1)

  CHARACTER(1),INTENT(IN)  :: CPACK(:)
  REAL,        INTENT(OUT) :: RDATA(:,:)
  INTEGER,     INTENT(IN)  :: NX,NY,NEXP
  REAL,        INTENT(IN)  :: VAR1

! only required when dealing with F95 compilers
! replace ICHAR below with internally defined JCHAR function
! CHARACTER MYCHR*1
! JCHAR(MYCHR)=IAND(ICHAR(MYCHR),255)

  SCALE=2.0**(7-NEXP)
  VOLD=VAR1
  INDX=0
  DO J=1,NY
    DO I=1,NX
      INDX=INDX+1
      RDATA(I,J)=(ICHAR(CPACK(INDX))-127.)/SCALE+VOLD
      VOLD=RDATA(I,J)
      !!	grid point north of La Palma at 29.00/-18.00
      IF(I.EQ.343.AND.J.EQ.120) WRITE(12,'(F10.2$)')RDATA(I,J)
    END DO
    VOLD=RDATA(1,J)
  END DO

END SUBROUTINE unpack
