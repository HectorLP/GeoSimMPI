 program comfile
  
      INTEGER*4               :: ncx, ncy, ncz
      INTEGER*4               :: ind, ix, iy, iz  
      INTEGER*4               :: inum, ijk, lnx, lny, lnz, N         
      REAL*4                  :: perm          
      CHARACTER (LEN=4)       :: intform
      CHARACTER (LEN=64)      :: filename, filenamea, filein, fileina, fileout, fileouta, fileout2
      DATA intform / '(I )' /
       
      write(*,*)'enter the number of compound bars'
      Read(*,*)  inum
      !inum =  40
      !filein = 'rpsults.out.'
      filein = 'pixdat.out.' !'pixasc.out.'
      fileina = 'pixdata.out.' 
      fileout = 'indicator.out'
      fileout2 = 'indicatorXYZ.out'
    
      OPEN(4, file=fileout, FORM='FORMATTED', STATUS='UNKNOWN')
      OPEN(5, file=fileout2, FORM='UNFORMATTED', STATUS='UNKNOWN')
!vlf
      OPEN(11, file='indicatorXYZ.asc', FORM='FORMATTED', STATUS='UNKNOWN')
!vlf
     j=0
     DO ijk = 1, inum
       ii=ijk-1
       write(*,*)
       WRITE(intform(3:3), '(I1)')icount_digits(ii)
       Length=LEN_TRIM (filein)
       filename(1:Length)=filein
       WRITE(filename(Length+1:), intform)ijk-1
       Lengtha=LEN_TRIM (fileina)
       filenamea(1:Lengtha)=fileina
       WRITE(filenamea(Lengtha+1:), intform)ijk-1
!vlf   WRITE(100,*)ii,icount_digits(ii),ijk
       WRITE(*,*)filename, filenamea        
       OPEN(2, file=filename, FORM='UNFORMATTED', STATUS='old')
       OPEN(3, file=filenamea, FORM='FORMATTED', STATUS='old')
!      OPEN(2, file=filename, FORM='FORMATTED', STATUS='old')
       READ(3,*)lnx, lny, lnz
       N = lnx*lny*lnz
!      print *, ijk,':', lnx,lny,lnz,n
!      DO i = 1, N
!        READ(2,*) ind !, ix, iy, iz                 
!      ENDDO
!      N = lnx*lny*lnz               
       DO i = 1, N
         READ(2) ind, ix, iy, iz        
         j = j+1
         WRITE(4) ind
         WRITE(5) ind , ix, iy, iz       
         WRITE(11,*) ind , ix, iy, iz       
       ENDDO
!      print *, 'ijk: ', j
       CLOSE(2)
       CLOSE(3)
       WRITE(*,*)'read file', ijk
      ENDDO
      CLOSE(4)
         
      WRITE(*,*)'FINISHED!!!'
      STOP
      END program comfile


      
!***********************************************************************
!  Function icount_digits counts the number of digits in an integer.
!***********************************************************************

      FUNCTION icount_digits(ic)

      INTEGER :: icount_digits,ic


      icount_digits = 0
      DO
        icount_digits = icount_digits+1
        ic = ic/10
        IF( ic <= 0 ) EXIT
      END DO

      END FUNCTION iCOUNT_DIGITS


