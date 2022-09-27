 program comfile
  
      INTEGER*4, ALLOCATABLE  :: ind(:,:,:)
      INTEGER*4, ALLOCATABLE :: ind1(:,:,:)
      INTEGER*4               :: nx, ny, nz      
      INTEGER*4               :: maxx, maxy, maxz
      INTEGER*4               :: minx, miny, minz
      INTEGER*4               :: ik1, ik2, ik3
      INTEGER*4               :: inx, iny, inz
      INTEGER*4               :: ijk, i, j, k,ii
      INTEGER*4               :: iflag
      REAL*4                  :: i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24,i25,i26,i27
      REAL*4                  :: am, bm, cm
      REAL*4                  :: x, y, z
      REAL*4                  :: dx, dy, dz, N
      CHARACTER (LEN=4)       :: intform
      CHARACTER (LEN=64)      :: filename, filenamea, filein, fileina, fileout, fileouta, filepercent
      DATA intform / '(I )' /
       

      Write(*,*)'Enter input file name'
      Read(*,*) filein
      !  filein = 'indicator.out' !'pixdat.out.0' !'pixdatS.out'
      !Write(*,*)'Enter input file name 1'
      !Read(*,*) fileina

      Write(*,*)'Enter output file name'
      Read(*,*) fileout
      !  fileout = 'm.out'
      !Write(*,*)'Enter output file name 1'
      !Read(*,*) fileouta

      Write(*,*)'Enter nx'
      Read(*,*) nx
      !  nx = 100 !200
      Write(*,*)'Enter ny'
      Read(*,*) ny
      !  ny = 100
      Write(*,*)'Enter nz'
      Read(*,*) nz
      !  nz = 120
     ! Write(*,*)'Do you need binary file? Enter (1=yes, 0=no)'
     ! Read(*,*)iflag
      iflag = 1
      filepercent = 'percent.out'

      N = nx*ny*nz
      !N = nx*ny*100
      ALLOCATE(ind(nx, ny, nz), STAT=IERR)
      ALLOCATE(ind1(nx, ny, nz), STAT=IERR)
         
      filename=filein 
      !filenamea=fileina
 
      IF(iflag==1) THEN
        OPEN(2, file=filename, FORM='UNFORMATTED', STATUS='old')
      ELSE
        OPEN(2, file=filename, FORM='FORMATTED', STATUS='old')
      ENDIF       
      !OPEN(3, file=filenamea, FORM='FORMATTED', STATUS='old')
      !READ(3,*)nx, ny, nz

      i0=0
      i1=0
      i2=0
      i3=0
      i4=0
      i5=0
      i6=0
      i7=0
      i8=0
      i9=0
      i10=0
      i11=0
      i12=0
      i13=0
      i14=0
      i15=0
      i16=0
      i17=0
      i18=0
      i19=0
      i20=0
      i21=0
      i22=0
      i23=0
      i24=0
      i25=0
      i26=0
      i27=0
 
      !READ(2)ind(1:nx,1:ny,1:nz)
      
      DO j = 1,ny
        DO i = 1, nx  
          DO k =1, nz
            IF(iflag==1) THEN
              READ(2)ind(i,j,k) !,x,y,z
!             ind1(i,j,k) = ind(i,j,k)
              if(ind(i,j,k) == 0)then
                i0=i0+1
              else if(ind(i,j,k) == 1)then
                i1=i1+1
              else if(ind(i,j,k) == 2)then
                i2=i2+1
              else if(ind(i,j,k) == 3)then
                i3=i3+1
              else if(ind(i,j,k) == 4)then
                i4=i4+1
              else if(ind(i,j,k) == 5)then
                i5=i5+1
              else if(ind(i,j,k) == 6)then
                i6=i6+1
              else if(ind(i,j,k) == 7)then
                i7=i7+1
              else if(ind(i,j,k) == 8)then
                i8=i8+1
              else if(ind(i,j,k) == 9)then
                i9=i9+1
              else if(ind(i,j,k) == 10)then
                i10=i10+1           
              else if(ind(i,j,k) == 11)then
                i11=i11+1
              else if(ind(i,j,k) == 12)then
                i12=i12+1
              else if(ind(i,j,k) == 13)then
                i13=i13+1
              else if(ind(i,j,k) == 14)then
                i14=i14+1
              else if(ind(i,j,k) == 15)then
                i15=i15+1
              else if(ind(i,j,k) == 16)then
                i16=i16+1
              else if(ind(i,j,k) == 17)then
                i17=i17+1
              else if(ind(i,j,k) == 18)then
                i18=i18+1
              else if(ind(i,j,k) == 19)then
                i19=i19+1
              else if(ind(i,j,k) == 20)then
                i20=i20+1
              else if(ind(i,j,k) == 22)then
                i22=i22+1
              else if(ind(i,j,k) == 23)then
                i23=i23+1
              else if(ind(i,j,k) == 24)then
                i24=i24+1
              else if(ind(i,j,k) == 25)then
                i25=i25+1
              else if(ind(i,j,k) == 26)then
                i26=i26+1
              else if(ind(i,j,k) == 27)then
                i27=i27+1           
              endif
            ELSE
              !READ(2,*)x,y,z,ind(i,j,k)
              READ(2,*)ind(i,j,k)
            ENDIF
           ENDDO
         ENDDO
       ENDDO

       CLOSE(2)
       CLOSE(3)

!      iflag = 0
       IF(iflag==1) THEN
         OPEN(4, file=fileout, FORM='UNFORMATTED', STATUS='UNKNOWN')
       ELSE
         OPEN(4, file=fileout, FORM='FORMATTED', STATUS='UNKNOWN')
       ENDIF
       OPEN(5, file=filepercent, FORM='FORMATTED', STATUS='UNKNOWN')
       !OPEN(5, file=fileouta, FORM='FORMATTED', STATUS='UNKNOWN')
!      WRITE(4,*)ind(1:nx,1:ny,1:nz)
       WRITE(4)ind(1:nx,1:ny,1:nz)
       !WRITE(4)ind(1:nx,1:ny,21:120)
       !WRITE(44,*)ind(1:nx,1:ny,1:nz)
        WRITE(5,*) 100*i0/N,100*i1/N,100*i2/N,100*i3/N,100*i4/N,100*i5/N,100*i6/N,100*i7/N,100*i8/N,100*i9/N,100*i10/N,100*i11/N,100*i12/N,100*i13/N,100*i14/N,100*i15/N,100*i16/N,100*i17/N,100*i18/N,100*i19/N,100*i20/N,100*i21/N,100*i22/N,100*i23/N,100*i24/N,100*i25/N,100*i26/N,100*i27/N, 'open gravel:', 100*(i4+i7+i10+i13)/N 
       !write(*,*)am,bm,cm
       WRITE(*,*)'FINISHED!!!'
       STOP
       END PROGRAM comfile

      
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
