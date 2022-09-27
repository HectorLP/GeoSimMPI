!  This is the main calling program for the scalable utilities !!linegen and tplane
!Ramya Ramanathan, Pacific Northwest Natioanl Lab/Wright State ! !University. Last modified 2010/07/20.
!**********************************************************

      MODULE prcmpi

        INTEGER*4                                         :: nproc
        INTEGER*4                                         :: iproc
        INTEGER*4                                         :: ncbar_partition
        INTEGER*4						                  :: nc_determine
  
        INTEGER*4                                         :: ibin
        INTEGER*4                                         :: itest
        INTEGER*4                                         :: ichk
        CHARACTER (LEN=4)                                 :: intform
  
        INTEGER*4, ALLOCATABLE                          :: ncbar_process(:)
  
        DATA intform / '(I )' /
  
        END MODULE prcmpi
  
  !**********************************************************
  
        MODULE bar
  
        INTEGER*4                                         :: nocb
        INTEGER*4                                         :: ncbar
  
        INTEGER*4                                         :: noub
        INTEGER*4                                         :: nocs
        INTEGER*4                                         :: nl
  
        INTEGER*4                                         :: jcm
        INTEGER*4, ALLOCATABLE                            :: jcmarray(:)
        INTEGER*4, ALLOCATABLE, DIMENSION(:)            :: nncoset(:) !Naum
        INTEGER*4, ALLOCATABLE, DIMENSION(:,:)          :: nline(:,:) !Naum
  
        INTEGER*2, ALLOCATABLE                            :: iub(:)
  
        INTEGER*4                                         :: is1
        INTEGER*4                                         :: is2
        INTEGER*4                                         :: is3
      
      INTEGER*4, ALLOCATABLE						    :: noub_array(:) !Array to save noub number for each component bar when ncbar > nproc
  
        REAL*4, ALLOCATABLE, DIMENSION(:,:,:)             :: xmin
        REAL*4, ALLOCATABLE, DIMENSION(:,:,:)             :: xmax
        REAL*4, ALLOCATABLE, DIMENSION(:,:,:)             :: delh
        REAL*4, ALLOCATABLE, DIMENSION(:,:,:)             :: rllen
  
        REAL*4, ALLOCATABLE, DIMENSION(:)                 :: wcb
        REAL*4, ALLOCATABLE, DIMENSION(:)                 :: lcb
        REAL*4, ALLOCATABLE, DIMENSION(:)                 :: hcb
  
        REAL*4, ALLOCATABLE, DIMENSION(:)                 :: x
        REAL*4, ALLOCATABLE, DIMENSION(:)                 :: y
        REAL*4, ALLOCATABLE, DIMENSION(:)                 :: z
  
        REAL*4, ALLOCATABLE, DIMENSION(:,:)               :: width(:)
        REAL*4, ALLOCATABLE, DIMENSION(:,:)               :: length(:)
        REAL*4, ALLOCATABLE                               :: zhmin(:)
        REAL*4, ALLOCATABLE                               :: zhmax(:)
        REAL*4, ALLOCATABLE, DIMENSION(:,:)               :: hu(:)
        REAL*4, ALLOCATABLE, DIMENSION(:,:)               :: zi(:)
        REAL*4, ALLOCATABLE, DIMENSION(:,:,:)             :: xinfl(:,:)
        REAL*4, ALLOCATABLE, DIMENSION(:,:,:)             :: yinfl(:,:)
        REAL*4, ALLOCATABLE, DIMENSION(:,:,:)             :: a(:,:)
        REAL*4, ALLOCATABLE, DIMENSION(:,:,:)             :: b(:,:)
  
        REAL*4                                            :: mxlub
        REAL*4                                            :: wxub
      REAL*4                                            :: llenm
      REAL*4                                            :: twidthm
      REAL*4                                            :: delhm
  
        END MODULE bar
  !**********************************************************
  
        PROGRAM MAIN
  
        USE prcmpi
        USE bar
  
        INCLUDE 'mpif.h'
  
        INTEGER*4                                         :: ierr
        INTEGER*4                                         :: ans
        INTEGER*4                                         :: i_partition
        INTEGER*4                                         :: i
        

        INTEGER*4						                              :: ncbar_partition1
        INTEGER*4                                         ::i_cbs 
        
        CHARACTER (LEN=15)                                :: adum
  
      1 FORMAT(1X,a)
     11 FORMAT(a12)
    101 FORMAT(1X,/)
  
  ! Flag for writing binary or ascii files - Need ascii for error checking; itest
  ! keeps random numbers the same on each processor; ichk = 1 will write debug files.
  
        ibin = 1
        itest = 0
        !ichk = 1 !Naum
        ichk = 0 !Naum
  
  !  Initialize MPI
  
        CALL MPI_INIT(ierr)
        adum = 'mpiInit'
        CALL CHECK(adum,ierr)
  
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
        adum = 'mpiSize'
        CALL CHECK(adum,ierr)
  
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
        adum = 'mpiRank'
        CALL CHECK(adum,ierr)
  
  !  Get Number of compound bars from maxv.out
  
        IF (iproc == 0 ) THEN
          OPEN(55,FILE='maxv.out',STATUS='old')
          READ(55,*)nocb
      ! nc_determine = MOD(nocb, nproc)
      ! IF (nc_determine == 0) THEN
      	  ncbar_partition = nocb/nproc
      ! ELSE
      !     	ncbar_partition = nocb / nproc + 1
      ! ENDIF
          !WRITE(*, *) 'The number for each cores is ', ncbar_partition
          CLOSE(55)
        ENDIF
        ! CLOSE(55)
  
        CALL MPI_Bcast(nocb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        adum = 'nocb-mpi'
        CALL CHECK(adum,ierr)
        
        ncbar = 1!nocb/nproc
        CALL MPI_Bcast(ncbar_partition,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        adum = 'ncbar_partition'
        CALL CHECK(adum, ierr)
  
        ALLOCATE(noub_array(nocb),STAT=IERR)
        ALLOCATE(jcmarray(nocb),STAT=IERR)
        ALLOCATE(ncbar_process(nproc),STAT=IERR)
        ! Establish a cartesian topology communicator
        nc_determine = MOD(nocb, nproc)
        ! ncbar_partition1 = nocb/nproc
        ! ncbar_process(iproc) = ncbar_partition
        ! IF (nc_determine > 0) THEN
        ncbar_process(iproc+1) = ncbar_partition
        IF (iproc < nc_determine) THEN
            ncbar_process(iproc+1) = ncbar_partition + 1
        ENDIF
        WRITE(*,*) 'The component bars in each process is ',iproc,ncbar_process(iproc+1), ncbar_partition
        ! DO i_cbs = 0, nproc-1
        !   IF (i_cbs==iproc) THEN
        !     ncbar_partition=ncbar_process(i_cbs+1)
        !   ENDIF
        ! ENDDO
        
        
        !WRITE(*,*) 'The process is ', iproc
        !DO i = 1,ncbar_partition
        !  WRITE(*,*) 'The index of files for each core is', iproc*ncbar_partition + i
        !END DO
        !IF (nproc /= nocb) THEN
        !  WRITE(*,*) 'Number of compound bars '
        !  WRITE(*,*) 'must be equal to the number of processors'
        !  WRITE(*,*) 'Number of compound bars: ',nocb
        !  WRITE(*,*) 'Number of processors: ',nproc
        !  CALL MPI_Finalize(ierr)
        !  adum = 'mpiFinalize'
        !  CALL CHECK(adum,ierr)
        !  STOP
        !ENDIF
  
  !  OK...This will only work for one compound bar per processor...
  
        !ncbar = nocb/nproc
       
        !ll = iproc
        WRITE(intform(3:3),'(I1)') icount_digits(ll)
        DO i = 0, ncbar_partition-1
            ll = iproc * ncbar_partition + i
            IF(ll < nocb) THEN
                !WRITE(*,*) 'Read the parameters from index ',ll
                WRITE(intform(3:3),'(I1)') icount_digits(ll)
            ENDIF
        END DO
  
  
  !  Set initial values of random number generator seeds, and
  !  print to output:
  
        IF (iproc == 0) THEN
          !WRITE(*,101)
          !WRITE(*,1) 'Seed values specified (1) or random (0)? '
          !READ(*,*) ans !Naum
  !vlf
          !ans=0 !Naum
          ans=3 !Naum
          is1 = 10
          is2 = 20
          is3 = 30
  !vlf
          IF (ans == 1) THEN
            WRITE(*,1) '   Enter first integer seed: '
            READ(*,*) is1
            WRITE(*,1) '   Enter second integer seed: '
            READ(*,*) is2
            WRITE(*,1) '   Enter third integer seed: '
            READ(*,*) is3
          ELSE
            CALL tseed(is1,is2,is3)
          END IF
          WRITE(*,101)
          WRITE(*,1) 'Random seed values used are: '
          WRITE(*,101)
          WRITE(*,*) '     IS1: ',is1
          WRITE(*,*) '     IS2: ',is2
          WRITE(*,*) '     IS3: ',is3
          WRITE(*,101)
          WRITE(*,101)
        ENDIF
  
        CALL MPI_Bcast(is1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        adum = 'is1-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(is2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        adum = 'is2-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(is3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        adum = 'is3-mpi'
        CALL CHECK(adum,ierr)
  
  
        CALL cbch     !1 start 216
        CALL MPI_Barrier(MPI_COMM_WORLD,i_error)
        CALL arcgen   !2 start 1176 
        CALL MPI_Barrier(MPI_COMM_WORLD,i_error)
        CALL trgena   !3 start 1729
        CALL MPI_Barrier(MPI_COMM_WORLD,i_error)
        CALL ubch     !4 start 2003
        CALL MPI_Barrier(MPI_COMM_WORLD,i_error)
        CALL trgenao  !5 start 3084 
        CALL MPI_Barrier(MPI_COMM_WORLD,i_error)
        CALL cbcch    !6 start 3418
        CALL MPI_Barrier(MPI_COMM_WORLD,i_error)
        CALL linegen  !7 start 3623
        CALL MPI_Barrier(MPI_COMM_WORLD,i_error)
        CALL tplane   !8 start 5365
  
        CALL MPI_FINALIZE(ierr)
        adum = 'mpiFinalize'
        CALL CHECK(adum,ierr)
  
        END !main end
  
  !***********************************************************************************************************************************
  ! Code converted using TO_F90 by Alan Miller
  ! Date: 2007-07-27  Time: 16:56:30
  
  !     Program CBCH generates coefficient of planes for each compound bar.
  !     Orginally created by, Tim Scheibe, Modified by Guin,17 May  ,06
  !     Modified Feb11,07, modified(march 28,07), April30,07
  !     n-number of CBs...random width, length and height
  !     add bottom plane..(May 17,07)
  !     Number of planes: 30
  
  !     Input :   cbpts.out  (output file from crtpts.f)
  !               cblwh.out  (output file from crtpts.f)
  
  !     Output:  Plane equation of CB (used as input file in pixdat)
  
  !      6/9/07...added MCF
  !***********************************************************************************************************************************
        SUBROUTINE cbch
  
        USE prcmpi
        USE bar
  
        INCLUDE 'mpif.h'
  
        REAL*4                                           :: thi1
        REAL*4                                           :: thi2
        REAL*4                                           :: thi3
        REAL*4                                           :: thi4
        REAL*4                                           :: thi5
        REAL*4                                           :: tript(3,3)
  
        REAL*4,ALLOCATABLE                               :: aaa(:)
        REAL*4,ALLOCATABLE                               :: bbb(:)
        REAL*4,ALLOCATABLE                               :: ccc(:)
        REAL*4,ALLOCATABLE                               :: ddd(:)
  
        REAL*4                                           :: l1
        REAL*4                                           :: l2
        REAL*4                                           :: l3
        REAL*4                                           :: capl
        REAL*4                                           :: w
        REAL*4                                           :: wa
        REAL*4                                           :: wb
        REAL*4                                           :: wc
        REAL*4                                           :: cjj
        REAL*4                                           :: cthi
        REAL*4                                           :: thia
        REAL*4                                           :: thib
        REAL*4                                           :: thic
        REAL*4                                           :: thid
        REAL*4                                           :: cy
        REAL*4                                           :: cx
        REAL*4                                           :: cl
        REAL*4                                           :: h0
        REAL*4                                           :: h1
        REAL*4                                           :: h2
        REAL*4                                           :: h3
        REAL*4                                           :: h
        REAL*4                                           :: xi
        REAL*4                                           :: yi
        REAL*4                                           :: th
        REAL*4                                           :: ys
        REAL*4                                           :: wml
        REAL*4                                           :: hm
        REAL*4                                           :: h4
        REAL*4                                           :: det
  
        REAL*4,ALLOCATABLE                               :: ba(:)
        REAL*4,ALLOCATABLE                               :: bb(:)
        REAL*4,ALLOCATABLE                               :: bc(:)
        REAL*4,ALLOCATABLE                               :: bd(:)
  
        REAL*4,ALLOCATABLE                               :: aml(:)
        REAL*4,ALLOCATABLE                               :: bml(:)
        REAL*4,ALLOCATABLE                               :: cml(:)
        REAL*4,ALLOCATABLE                               :: dml(:)
  
        REAL*4,ALLOCATABLE                               :: amr(:)
        REAL*4,ALLOCATABLE                               :: bmr(:)
        REAL*4,ALLOCATABLE                               :: cmr(:)
        REAL*4,ALLOCATABLE                               :: dmr(:)
  
        REAL*4, ALLOCATABLE                              :: cb(:)
        REAL*4, ALLOCATABLE                              :: alpha(:)
        INTEGER*2, ALLOCATABLE                           :: l(:) !Naum was 4
        INTEGER*2, ALLOCATABLE                           :: ind(:) !Naum was 4
  
        INTEGER*4                                        :: jj
        INTEGER*4                                        :: j
        INTEGER*2                                        :: m !Naum was 4
        INTEGER*4                                        :: kk
        INTEGER*4                                        :: ii
        INTEGER*4                                        :: ml
        INTEGER*4                                        :: mr
        
        INTEGER*4                                        :: i_partition
        
  
        CHARACTER (LEN=64)                               :: filename
        CHARACTER (LEN=15)                                :: adum
  
      1 FORMAT(1X,a)
     10 FORMAT(a12)
     15 FORMAT(1X,6(f8.6,2X))
     16 FORMAT(1X,7(f8.6,2X))
     20 FORMAT(1X,3(f11.8,3X))
     25 FORMAT(1X,4(f25.12,3X))
     30 FORMAT(1X,a,i1,a,f15.7)
     31 FORMAT(2(i4,2X),2X,3(f13.2),2X)
     32 FORMAT(2(i4,2X),4(f13.2,2X),i4)
     33 FORMAT(3(f13.2),2X)
     35 FORMAT(1X,9(f15.4,2X))
     36 FORMAT(12(f13.2,1X))
     37 FORMAT(i4,2X,5(f13.2,1X),i4)
    101 FORMAT(1X,/)
  
        IF (iproc == 0 ) THEN
          WRITE(*,101)
          WRITE(*,101)
          WRITE(*,1) 'Set Characteristic Geometry: cbch'
          WRITE(*,1) '****************************'
          WRITE(*,101)
          WRITE(*,1)'Thickness of bottom bed'
          !READ(*,*)th !Naum
          th=0.5 !Naum
        ENDIF
  
        CALL MPI_Bcast(th,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'th-mpi'
        CALL CHECK(adum,ierr)
  
  ! Input files
  
        !OPEN(3,FILE='cbpts.out',FORM='unformatted',STATUS='old') !Naum      
        !OPEN(77,FILE='cblwh.out',FORM='unformatted',STATUS='old') !Naum
  
  ! Read input data
        
        ALLOCATE (l(nocb), STAT=IERR)
        adum = 'l'
        CALL CHECK(adum,ierr)
  
        ALLOCATE (x(nocb), STAT=IERR)
        adum = 'x'
        CALL CHECK(adum,ierr)
  
        ALLOCATE (y(nocb), STAT=IERR)
        adum = 'y'
        CALL CHECK(adum,ierr)
  
        ALLOCATE (z(nocb), STAT=IERR)
        adum = 'z'
        CALL CHECK(adum,ierr)
  
        ALLOCATE (wcb(nocb), STAT=IERR)
        adum = 'wcb'
        CALL CHECK(adum,ierr)
  
        ALLOCATE (lcb(nocb), STAT=IERR)
        adum = 'lcb'
        CALL CHECK(adum,ierr)
  
        ALLOCATE (hcb(nocb), STAT=IERR)
        adum = 'hcb'
        CALL CHECK(adum,ierr)
  
        ALLOCATE (alpha(nocb), STAT=IERR)
        adum = 'alpha'
        CALL CHECK(adum,ierr)
  
        ALLOCATE (ind(nocb), STAT=IERR)
        adum = 'ind'
        CALL CHECK(adum,ierr)
  
        ALLOCATE(aaa(15),STAT=IERR)
        adum = 'aaa'
        CALL CHECK(adum,ierr)
  
        ALLOCATE(bbb(15),STAT=IERR)
        adum = 'bbb'
        CALL CHECK(adum,ierr)
  
        ALLOCATE(ccc(15),STAT=IERR)
        adum = 'ccc'
        CALL CHECK(adum,ierr)
  
        ALLOCATE(ddd(15),STAT=IERR)
        adum = 'ddd'
        CALL CHECK(adum,ierr)
  
  
        ALLOCATE(ba(15),STAT=IERR)
        adum = 'ba'
        CALL CHECK(adum,ierr)
  
        ALLOCATE(bb(15),STAT=IERR)
        adum = 'bb'
        CALL CHECK(adum,ierr)
  
        ALLOCATE(bc(15),STAT=IERR)
        adum = 'bc'
        CALL CHECK(adum,ierr)
  
        ALLOCATE(bd(15),STAT=IERR)
        adum = 'bd'
        CALL CHECK(adum,ierr)
  
        ALLOCATE(aml(3),STAT=IERR)
        adum = 'aml'
        CALL CHECK(adum,ierr)
  
        ALLOCATE(bml(3),STAT=IERR)
        adum = 'bml'
        CALL CHECK(adum,ierr)
  
        ALLOCATE(cml(3),STAT=IERR)
        adum = 'cml'
        CALL CHECK(adum,ierr)
  
        ALLOCATE(dml(3),STAT=IERR)
        adum = 'dml'
        CALL CHECK(adum,ierr)
  
        ALLOCATE(amr(3),STAT=IERR)
        adum = 'amr'
        CALL CHECK(adum,ierr)
  
        ALLOCATE(bmr(3),STAT=IERR)
        adum = 'bmr'
        CALL CHECK(adum,ierr)
  
        ALLOCATE(cmr(3),STAT=IERR)
        adum = 'cmr'
        CALL CHECK(adum,ierr)
  
        ALLOCATE(dmr(3),STAT=IERR)
        adum = 'dmr'
        CALL CHECK(adum,ierr)
        WRITE(*,*) 'It is the process ',iproc
        WRITE(*,*) 'The partition size is ',ncbar_partition
        !nn = iproc+1
        !DO  k = 1, nn !Naum
        !    WRITE(*,*) 'The index of parameters is ',k
        !    READ(77) m,l(k),wcb(k),lcb(k),hcb(k),alpha(k),ind(k) !Naum         
        !    READ(3) m,l(k),x(k),y(k),z(k) !Naum 
        !END DO
        !Change the code for processing the more number of component bars 
        ! OPEN(3,FILE='cbpts.out',FORM='unformatted',STATUS='old')
        ! OPEN(77,FILE='cblwh.out',FORM='unformatted',STATUS='old')
        DO i_partition = 1,  ncbar_process(iproc+1)
            ! nn = iproc * ncbar_partition + i_partition
            IF (iproc<nc_determine) THEN
                nn=iproc*(ncbar_partition+1)+i_partition
              ELSE
                ! kk1=myY*ncbar_partition+i_nb-1
                nn= nc_determine*(ncbar_partition+1)+&
                      (iproc-nc_determine)*ncbar_partition+i_partition
              END IF
            !WRITE(*,*) 'The index of the parameters file is ',nn
            !WRITE(*,*) 'The number of the component bars is ',nocb
            IF (nn <= nocb) THEN
              OPEN(3,FILE='cbpts.out',FORM='unformatted',STATUS='old')
              OPEN(77,FILE='cblwh.out',FORM='unformatted',STATUS='old')
                DO  k = 1, nn !Naum
        !        IF (k < nocb) THEN
                    !WRITE(*,*) 'The index of the parameter files ',k
                    READ(77) m,l(k),wcb(k),lcb(k),hcb(k),alpha(k),ind(k) !Naum         
                    READ(3) m,l(k),x(k),y(k),z(k) !Naum 
                    WRITE(*,*) 'The parameters from k: ', k, m, iproc
                    WRITE(*,*) l(k), wcb(k), lcb(k), hcb(k), alpha(k), ind(k)
                    WRITE(*,*) l(k), x(k), y(k), z(k)
                !ENDIF
                END DO
                !
                iiii = nn-1
                WRITE(intform(3:3),'(I1)') icount_digits(iiii)
                filename(1:11)='points.out.'
                WRITE(filename(12:),intform) nn-1
                IF (ibin==1) THEN
                          OPEN(14,FILE=filename,FORM='unformatted',STATUS='unknown')
                          WRITE(14)x(nn),y(nn),wcb(nn),lcb(nn),alpha(nn),ind(nn)
                ELSE
                          OPEN(14,FILE=filename,FORM='unformatted',STATUS='unknown')
                          WRITE(14,*)x(nn),y(nn),wcb(nn),lcb(nn),alpha(nn),ind(nn)
                ENDIF
                CLOSE(14)
        !     ENDIF
        ! END DO
          CLOSE(3)
          CLOSE(77)
  
        ! DO i_partition = 1, ncbar_process(iproc+1)
        !     ! kk = iproc * ncbar_partition + i_partition
        !     IF (iproc<nc_determine) THEN
        !         kk=iproc*(ncbar_partition+1)+i_partition
        !     ELSE
        !         ! kk1=myY*ncbar_partition+i_nb-1
        !         kk= nc_determine*(ncbar_partition+1)+&
        !               (iproc-nc_determine)*ncbar_partition+i_partition
        !     END IF
        !   IF (kk <= nocb) THEN
          kk = nn
          !WRITE(*,*) 'Points for the component bar with the index',kk
          l1=lcb(kk)/3.0
          l2=lcb(kk)/3.0
          l3=lcb(kk)/3.0
          wa=(2*wcb(kk))/5.0                                     !WSU 4/4/08
          wb=wcb(kk)
          wc=(2*wcb(kk))/3.0
          h0=z(kk)-hcb(kk)
          h1=hcb(kk)
          h2=hcb(kk)
          h3=hcb(kk)
          h4=hcb(kk)
          wml = wcb(kk)/12+10                                     !WSU 5/19/08 
    ! Plane 1:
          WRITE(*,*) 'For the Plane 1 from the component bar ',kk
          tript(1,1) = wb/2.0-wa/2.0
          tript(1,2) = 0.0
          tript(1,3) = h0
          tript(2,1) = tript(1,1)+ wa/12.0
          tript(2,2) = 0.0
          tript(2,3) = h0+h1/1.33
          tript(3,1) = 0.0
          tript(3,2) = l1
          tript(3,3) = h0
          jj = 1
          CALL det3d(tript,det)
          ddd(jj) = -det
          CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
          
    ! Plane 2:
          !WRITE(*,*) 'For the Plane 2 from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1) + wa/12.0
          tript(2,2) = 0.0
          tript(2,3) = h0+h1
          tript(3,1) = tript(3,1)+wb/12.0
          tript(3,2) = l1
          tript(3,3) = h0+h2/1.33
          jj = 2
          CALL det3d(tript,det)
          ddd(jj) = -det
          CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
    
    ! Plane 3:
          !WRITE(*,*) 'For the Plane 3 from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1) + (2*wa)/3.0
          tript(2,2) = 0.0
          tript(2,3) = h0+h1
          tript(3,1) = tript(3,1)+wb/12.0
          tript(3,2) = l1
          tript(3,3) = h0+h2
          jj = 3
          CALL det3d(tript,det)
          ddd(jj) = -det
          CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
    ! Plane 4
          !WRITE(*,*) 'For the Plane 4 from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1) + wa/12.0
          tript(2,2) = 0.0
          tript(2,3) = h0+h1/1.33
          tript(3,1) = tript(3,1)+(2*wb)/3.0
          tript(3,2) = l1
          tript(3,3) = h0+h2
          jj = 4
          CALL det3d(tript,det)
          ddd(jj) = -det
          CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
  
    ! Plane 5
          !WRITE(*,*) 'For the Plane 5 from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1) + wa/12.0
          tript(2,2) = 0.0
          tript(2,3) = h0
          tript(3,1) = tript(3,1)+wb/6.0
          tript(3,2) = l1
          tript(3,3) = h0
          jj = 5
          CALL det3d(tript,det)
          ddd(jj) = -det
          CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
  
    ! Plane 6
          !WRITE(*,*) 'For the Plane 6 from the component bar ',kk
          tript(1,1) = 0.0
          tript(1,2) = l1
          tript(1,3) = h0
          tript(2,1) = tript(1,1) + wb/12.0
          tript(2,2) = l1
          tript(2,3) = h0+h2/1.33
          tript(3,1) = wb/2.0-wc/2.0
          tript(3,2) = l1+l2
          tript(3,3) = h0
          jj = 6
          CALL det3d(tript,det)
          ddd(jj) = -det
          CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
  
    ! Plane 7
          !WRITE(*,*) 'For the Plane 7 from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1) + wb/12.0
          tript(2,2) = l1
          tript(2,3) = h0+h2
          tript(3,1) = tript(3,1)+wc/12.0
          tript(3,2) = l1+l2
          tript(3,3) = h0+h3/1.33
          jj = 7
          CALL det3d(tript,det)
          ddd(jj) = -det
          CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
  
    ! Plane 8
          !WRITE(*,*) 'For the Plane 8 from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1) + (2*wb)/3.0
          tript(2,2) = l1
          tript(2,3) = h0 +h2
          tript(3,1) = tript(3,1)+wc/12.0
          tript(3,2) = l1+l2
          tript(3,3) = h0+h3
          jj = 8
          CALL det3d(tript,det)
          ddd(jj) = -det
          CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
  
    ! Plane 9
          !WRITE(*,*) 'For the Plane 9 from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1) + wb/12.0
          tript(2,2) = l1
          tript(2,3) = h0 +h2/1.33
          tript(3,1) = tript(3,1)+(2*wc)/3.0
          tript(3,2) = l1+l2
          tript(3,3) = h0+h3
          jj = 9
          CALL det3d(tript,det)
          ddd(jj) = -det
          CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
  
    !  Plane 10
          !WRITE(*,*) 'For the Plane 10 from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1) + wb/12.0
          tript(2,2) = l1
          tript(2,3) = h0
          tript(3,1) = tript(3,1)+wc/6.0
          tript(3,2) = l1+l2
          tript(3,3) = h0
          jj = 10
          CALL det3d(tript,det)
          ddd(jj) = -det
          CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
  
    ! Plane 11
          !WRITE(*,*) 'For the Plane 11 from the component bar ',kk
          tript(1,1) = wb/2.0-wc/2.0
          tript(1,2) = l1+l2
          tript(1,3) = h0
          tript(2,1) = tript(1,1) + wc/12.0
          tript(2,2) = l1+l2
          tript(2,3) = h0+h3/1.33
          tript(3,1) = wb/2.0
          tript(3,2) = l1+l2+l3
          tript(3,3) = h0
          jj = 11
          CALL det3d(tript,det)
          ddd(jj) = -det
          CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
  
    ! Plane 12
          !WRITE(*,*) 'For the Plane 12 from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1)+wc/12.0
          tript(2,2) = l1+l2
          tript(2,3) = h0+h3
          jj = 12
          CALL det3d(tript,det)
          ddd(jj) = -det
          CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
  
    ! Plane 13
          !WRITE(*,*) 'For the Plane 13 from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1)+(2*wc)/3.0
          tript(2,2) = l1+l2
          tript(2,3) = h0+h3
          jj = 13
          CALL det3d(tript,det)
          ddd(jj) = -det
          CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
  
    ! Plane 14
          !WRITE(*,*) 'For the Plane 14 from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1)+wc/12.0
          tript(2,2) = l1+l2
          tript(2,3) = h0+h3/1.33
          jj = 14
          CALL det3d(tript,det)
          ddd(jj) = -det
          CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
  
    ! Plane 15
          !WRITE(*,*) 'For the Plane 15 from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1)+wc/12.0
          tript(2,2) = l1+l2
          tript(2,3) = h0
          jj = 15
          CALL det3d(tript,det)
          ddd(jj) = -det
          CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
  
    ! Bottom planes
  
  
    ! Plane 1a:
          !WRITE(*,*) 'For the Plane 1a - bottom planes from the component bar ',kk
          tript(1,1) = wb/2.0-wa/2.0
          tript(1,2) = 0.0
          tript(1,3) = h0-th
          tript(2,1) = tript(1,1)+ wa/12.0
          tript(2,2) = 0.0
          tript(2,3) = (h0+h1/1.33)-th
          tript(3,1) = 0.0
          tript(3,2) = l1
          tript(3,3) = h0-th
          ii = 1
          CALL det3d(tript,det)
          bd(ii) = -det
          CALL findabc(tript,ba(ii),bb(ii),bc(ii))
  
    ! Plane 2a:
          !WRITE(*,*) 'For the Plane 2a - bottom planes from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1) + wa/12.0
          tript(2,2) = 0.0
          tript(2,3) = (h0+h1)-th
          tript(3,1) = tript(3,1)+wb/12.0
          tript(3,2) = l1
          tript(3,3) = (h0+h2/1.33)-th
          ii = 2
          CALL det3d(tript,det)
          bd(ii) = -det
          CALL findabc(tript,ba(ii),bb(ii),bc(ii))
  
    ! Plane 3a:
          !WRITE(*,*) 'For the Plane 3a - bottom planes from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1) + (2*wa)/3.0
          tript(2,2) = 0.0
          tript(2,3) = (h0+h1)-th
          tript(3,1) = tript(3,1)+wb/12.0
          tript(3,2) = l1
          tript(3,3) = (h0+h2)-th
          ii = 3
          CALL det3d(tript,det)
          bd(ii) = -det
          CALL findabc(tript,ba(ii),bb(ii),bc(ii))
  
    ! Plane 4a
          !WRITE(*,*) 'For the Plane 4a - bottom planes from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1) + wa/12.0
          tript(2,2) = 0.0
          tript(2,3) = (h0+h1/1.33)-th
          tript(3,1) = tript(3,1)+(2*wb)/3.0
          tript(3,2) = l1
          tript(3,3) = (h0+h2)-th
          ii = 4
          CALL det3d(tript,det)
          bd(ii) = -det
          CALL findabc(tript,ba(ii),bb(ii),bc(ii))
  
    ! Plane 5a
          !WRITE(*,*) 'For the Plane 5a - bottom planes from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1) + wa/12.0
          tript(2,2) = 0.0
          tript(2,3) = h0-th
          tript(3,1) = tript(3,1)+wb/6.0
          tript(3,2) = l1
          tript(3,3) = h0-th
          ii = 5
          CALL det3d(tript,det)
          bd(ii) = -det
          CALL findabc(tript,ba(ii),bb(ii),bc(ii))
  
    ! Plane 6a
          !WRITE(*,*) 'For the Plane 6a - bottom planes from the component bar ',kk
          tript(1,1) = 0.0
          tript(1,2) = l1
          tript(1,3) = h0-th
          tript(2,1) = tript(1,1) + wb/12.0
          tript(2,2) = l1
          tript(2,3) = (h0+h2/1.33)-th
          tript(3,1) = wb/2.0-wc/2.0
          tript(3,2) = l1+l2
          tript(3,3) = h0-th
          ii = 6
          CALL det3d(tript,det)
          bd(ii) = -det
          CALL findabc(tript,ba(ii),bb(ii),bc(ii))
  
    ! Plane 7a
          !WRITE(*,*) 'For the Plane 7a - bottom planes from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1) + wb/12.0
          tript(2,2) = l1
          tript(2,3) = (h0+h2)-th
          tript(3,1) = tript(3,1)+wc/12.0
          tript(3,2) = l1+l2
          tript(3,3) = (h0+h3/1.33)-th
          ii = 7
          CALL det3d(tript,det)
          bd(ii) = -det
          CALL findabc(tript,ba(ii),bb(ii),bc(ii))
  
    ! Plane 8a
          !WRITE(*,*) 'For the Plane 8a - bottom planes from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1) + (2*wb)/3.0
          tript(2,2) = l1
          tript(2,3) = h0 +h2-th
          tript(3,1) = tript(3,1)+wc/12.0
          tript(3,2) = l1+l2
          tript(3,3) = h0+h3-th
          ii = 8
          CALL det3d(tript,det)
          bd(ii) = -det
          CALL findabc(tript,ba(ii),bb(ii),bc(ii))
  
    ! Plane 9a
          !WRITE(*,*) 'For the Plane 9a - bottom planes from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1) + wb/12.0
          tript(2,2) = l1
          tript(2,3) = (h0 +h2/1.33)-th
          tript(3,1) = tript(3,1)+(2*wc)/3.0
          tript(3,2) = l1+l2
          tript(3,3) = (h0+h3)-th
          ii = 9
          CALL det3d(tript,det)
          bd(ii) = -det
          CALL findabc(tript,ba(ii),bb(ii),bc(ii))
  
    ! Plane 10a
          !WRITE(*,*) 'For the Plane 10a - bottom planes from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1) + wb/12.0
          tript(2,2) = l1
          tript(2,3) = h0-th
          tript(3,1) = tript(3,1)+wc/6.0
          tript(3,2) = l1+l2
          tript(3,3) = h0-th
          ii = 10
          CALL det3d(tript,det)
          bd(ii) = -det
          CALL findabc(tript,ba(ii),bb(ii),bc(ii))
  
    ! Plane 11a
          !WRITE(*,*) 'For the Plane 11a - bottom planes from the component bar ',kk
          tript(1,1) = wb/2.0-wc/2.0
          tript(1,2) = l1+l2
          tript(1,3) = h0-th
          tript(2,1) = tript(1,1) + wc/12.0
          tript(2,2) = l1+l2
          tript(2,3) = (h0+h3/1.33)-th
          tript(3,1) = wb/2.0
          tript(3,2) = l1+l2+l3
          tript(3,3) = h0-th
          ii = 11
          CALL det3d(tript,det)
          bd(ii) = -det
          CALL findabc(tript,ba(ii),bb(ii),bc(ii))
  
    ! Plane 12a
          !WRITE(*,*) 'For the Plane 12a - bottom planes from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1)+wc/12.0
          tript(2,2) = l1+l2
          tript(2,3) = h0+h3-th
          ii = 12
          CALL det3d(tript,det)
          bd(ii) = -det
          CALL findabc(tript,ba(ii),bb(ii),bc(ii))
  
    ! Plane 13a
          !WRITE(*,*) 'For the Plane 13a - bottom planes from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1)+(2*wc)/3.0
          tript(2,2) = l1+l2
          tript(2,3) = h0+h3-th
          ii = 13
          CALL det3d(tript,det)
          bd(ii) = -det
          CALL findabc(tript,ba(ii),bb(ii),bc(ii))
  
    ! Plane 14a
          !WRITE(*,*) 'For the Plane 14a - bottom planes from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1)+wc/12.0
          tript(2,2) = l1+l2
          tript(2,3) = (h0+h3/1.33)-th
          ii = 14
          CALL det3d(tript,det)
          bd(ii) = -det
          CALL findabc(tript,ba(ii),bb(ii),bc(ii))
  
    ! Plane 15a
          !WRITE(*,*) 'For the Plane 15a - bottom planes from the component bar ',kk
          DO  k = 1, 3
      tript(1,k) = tript(2,k)
          END DO
          tript(2,1) = tript(1,1)+wc/12.0
          tript(2,2) = l1+l2
          tript(2,3) = h0-th
          ii = 15
          CALL det3d(tript,det)
          bd(ii) = -det
          CALL findabc(tript,ba(ii),bb(ii),bc(ii))
    ! Mcfl1
          !WRITE(*,*) 'For the Mcf1 from the component bar ',kk
          tript(1,1) = wb/2.0-wa/2.0+wa/12.0
          tript(1,2) = 0.0
          tript(1,3) = h0+h1/1.33
          tript(2,1) = tript(1,1)+ (wml-wa/12.0)
          tript(2,2) = 0.0
          tript(2,3) = h0
          tript(3,1) = wml
          tript(3,2) = l1
          tript(3,3) = h0
          ml = 1
          CALL det3d(tript,det)
          dml(ml) = -det
          CALL findabc(tript,aml(ml),bml(ml),cml(ml))
  
    ! Mcfl2
          !WRITE(*,*) 'For the Mcf12 from the component bar ',kk
          tript(1,1) = wb/12.0
          tript(1,2) = l1
          tript(1,3) = h0+h2/1.33
          tript(2,1) = tript(1,1)+ (wml-wb/12.0)
          tript(2,2) = l1
          tript(2,3) = h0
          tript(3,1) = wb/2.0-wc/2.0+wml
          tript(3,2) = l1+l2
          tript(3,3) = h0
          ml = 2
          CALL det3d(tript,det)
          dml(ml) = -det
          CALL findabc(tript,aml(ml),bml(ml),cml(ml))
  
    ! Mcfl3
          !WRITE(*,*) 'For the Mcf13 from the component bar ',kk
          tript(1,1) = wb/2.0-wc/2.0+wc/12.0
          tript(1,2) = l1+l2
          tript(1,3) = h0+h3/1.33
          tript(2,1) = tript(1,1)+ (wml-wc/12.0)
          tript(2,2) = l1+l2
          tript(2,3) = h0
          tript(3,1) = wb/2.0 + wml                                           !WSU 4/4/08
          tript(3,2) = l1+l2+l3
          tript(3,3) = h0
          ml = 3
          CALL det3d(tript,det)
          dml(ml) = -det
          CALL findabc(tript,aml(ml),bml(ml),cml(ml))
  
    ! Mcfr1
          !WRITE(*,*) 'For the Mcfr1 from the component bar ',kk
          tript(1,1) = (wb/2.0+wa/2.0)-wml
          tript(1,2) = 0.0
          tript(1,3) = h0
          tript(2,1) = tript(1,1)+ (wml-wa/12.0)
          tript(2,2) = 0.0
          tript(2,3) = h0+h1/1.33
          tript(3,1) = wb-wml
          tript(3,2) = l1
          tript(3,3) = h0
          mr = 1
          CALL det3d(tript,det)
          dmr(mr) = -det
          CALL findabc(tript,amr(mr),bmr(mr),cmr(mr))
  
    ! Mcfr2
          !WRITE(*,*) 'For the Mcfr2 from the component bar ',kk
          tript(1,1) = wb-wml
          tript(1,2) = l1
          tript(1,3) = h0
          tript(2,1) = tript(1,1)+ wml-wb/12.0
          tript(2,2) = l1
          tript(2,3) = h0+h2/1.33
          tript(3,1) = (wb/2.0+wc/2.0)-wml
          tript(3,2) = l1+l2
          tript(3,3) = h0
          mr = 2
          CALL det3d(tript,det)
          dmr(mr) = -det
          CALL findabc(tript,amr(mr),bmr(mr),cmr(mr))
  
    ! Mcfr3
          !WRITE(*,*) 'For the Mcfr3 from the component bar ',kk
          tript(1,1) = (wb/2.0+wc/2.0)-wml
          tript(1,2) = l1+l2
          tript(1,3) = h0
          tript(2,1) = tript(1,1)+ wml-wc/12.0
          tript(2,2) = l1+l2
          tript(2,3) = h0+h3/1.33
          tript(3,1) = wb/2.0-wml                                                     !WSU 4/4/08       
          tript(3,2) = l1+l2+l3
          tript(3,3) = h0
          mr = 3
          CALL det3d(tript,det)
          dmr(mr) = -det
          CALL findabc(tript,amr(mr),bmr(mr),cmr(mr))
          
          iiii = kk - 1
          WRITE(intform(3:3),'(I1)') icount_digits(iiii)
          filename(1:9) = 'cbch.out.'
          WRITE(filename(10:),intform) kk-1 !iproc
          OPEN(81,FILE=filename,FORM='unformatted',STATUS='unknown')
  
          WRITE(81)aaa      
          WRITE(81)bbb
          WRITE(81)ccc
          WRITE(81)ddd
  
          WRITE(81)ba
          WRITE(81)bb
          WRITE(81)bc
          WRITE(81)bd
  
          WRITE(81)aml
          WRITE(81)bml
          WRITE(81)cml
          WRITE(81)dml
  
          WRITE(81)amr
          WRITE(81)bmr
          WRITE(81)cmr
          WRITE(81)dmr
          CLOSE(81)
      ENDIF
      END DO
  
        DEALLOCATE (aaa,STAT=IERR)
        adum = 'aaa-de'
        CALL CHECK(adum,ierr)
  
        DEALLOCATE (bbb,STAT=IERR)
        adum = 'bbb-de'
        CALL CHECK(adum,ierr)
  
        DEALLOCATE (ccc,STAT=IERR)
        adum = 'ccc-de'
        CALL CHECK(adum,ierr)
  
        DEALLOCATE (ddd,STAT=IERR)
        adum = 'ddd-de'
        CALL CHECK(adum,ierr)
  
  
        DEALLOCATE (ba,STAT=IERR)
        adum = 'ba-de'
        CALL CHECK(adum,ierr)
  
        DEALLOCATE (bb,STAT=IERR)
        adum = 'bb-de'
        CALL CHECK(adum,ierr)
  
        DEALLOCATE (bc,STAT=IERR)
        adum = 'bc-de'
        CALL CHECK(adum,ierr)
  
        DEALLOCATE (bd,STAT=IERR)
        adum = 'bd-de'
        CALL CHECK(adum,ierr)
  
  
        DEALLOCATE (aml,STAT=IERR)
        adum = 'aml-de'
        CALL CHECK(adum,ierr)
  
        DEALLOCATE (bml,STAT=IERR)
        adum = 'bml-de'
        CALL CHECK(adum,ierr)
  
        DEALLOCATE (cml,STAT=IERR)
        adum = 'cml-de'
        CALL CHECK(adum,ierr)
  
        DEALLOCATE (dml,STAT=IERR)
        adum = 'dml-de'
        CALL CHECK(adum,ierr)
  
  
        DEALLOCATE (amr,STAT=IERR)
        adum = 'amr-de'
        CALL CHECK(adum,ierr)
  
        DEALLOCATE (bmr,STAT=IERR)
        adum = 'bmr-de'
        CALL CHECK(adum,ierr)
  
        DEALLOCATE (cmr,STAT=IERR)
        adum = 'cmr-de'
        CALL CHECK(adum,ierr)
  
        DEALLOCATE (dmr,STAT=IERR)
        adum = 'dmr-de'
        CALL CHECK(adum,ierr)
  
  
        !DEALLOCATE (l,STAT=IERR)
        !adum = 'l-de'
        !CALL CHECK(adum,ierr)
  
        DEALLOCATE (alpha,STAT=IERR)
        adum = 'alpha-de'
        CALL CHECK(adum,ierr)
  
        
        !CLOSE(3) 
        !CLOSE(77) 
        !CLOSE(81)
  
        RETURN
        END SUBROUTINE cbch
   
  !********************************************************************************
   
  ! Code converted using TO_F90 by Alan Miller
  ! Date: 2007-07-27  Time: 15:56:07
   
  ! arcgen.for, the code generates input files for TRGENA. This code will generate
  ! inflection points of the center arc of n-number of unit bars randomly.
  ! The code creates middle first then left and right until xmax and xmin reached.
  ! Created on November 3rd,06(Guin)..almost full portion of UBs are preserved in 
  ! downstream part of CB...latest version (April 19,07),
  ! iub = 1 (middle); =2(left); =3(right)
  !
  ! Latest (April 25,07)..creates UBs in different layers along Z.
  ! Read input files from crtpts.out
  ! Output files: 1- Inflection points
  !                  2- UBZi.out
  !                  3- NUB.out
  
  !*********************************************************************************
        SUBROUTINE arcgen
  
          USE prcmpi
          USE bar
      
          INCLUDE 'mpif.h'
      
          REAL*4                                            :: xi
          REAL*4                                            :: yi
          REAL*4                                            :: xx
          REAL*4                                            :: yy
          REAL*4                                            :: l
          REAL*4                                            :: wx
          REAL*4                                            :: w
      !      REAL*4                                            :: aa
          REAL*4                                            :: lx
          REAL*4                                            :: rx
          REAL*4                                            :: yii
          REAL*4                                            :: h
          REAL*4                                            :: xc
          REAL*4                                            :: yc
          REAL*4                                            :: zc
          REAL*4                                            :: zx
          REAL*4                                            :: wa                       !WSU 4/4/08                         
          REAL*4                                            :: wb                       !WSU 4/4/08 
          REAL*4                                            :: wc                       !WSU 4/4/08 
          REAL*4                                            :: l1                       !WSU 4/4/08 
          REAL*4                                            :: l2                       !WSU 4/4/08
          REAL*4                                            :: l3                        !WSU 4/4/08
          REAL*4                                            :: ya                        !WSU 4/4/08
          REAL*4                                            :: xa                       !WSU 4/4/08
          REAL*4                                            :: ang                       !WSU 5/13/08
          REAL*4                                            :: angcb1                    !WSU 06-19-08
          REAL*4                                            :: angcb2                    !WSU 06-19-08
          REAL*4                                            :: delangs                    !WSU 06-19-08
          REAL*4                                            :: delango                    !WSU 06-19-08
          REAL*4                                            :: atii                    !WSU 06-19-08
          REAL*4                                            :: ango                    !WSU 06-19-08
          REAL*4                                            :: angs                    !WSU 06-19-08
          REAL*4                                            :: humn
          REAL*4                                            :: huvar
          REAL*4                                            :: wumn
          REAL*4                                            :: wuvar
          REAL*4                                            :: lumn
          REAL*4                                            :: luvar
      !      REAL*4                                            :: xxmn
      !      REAL*4                                            :: xxvar
          REAL*4                                            :: yymn
          REAL*4                                            :: yyvar
      !      REAL*4                                            :: aamn
      !      REAL*4                                            :: aavar
      
          REAL*4, ALLOCATABLE                               :: xxx(:)
          REAL*4, ALLOCATABLE                               :: yyy(:)
      
          REAL*4, ALLOCATABLE                               :: ap(:)
          REAL*4, ALLOCATABLE                               :: bp(:)
          REAL*4, ALLOCATABLE                               :: x1(:)
      
          REAL*4, ALLOCATABLE                               :: x2(:)
      !      REAL*4, ALLOCATABLE                               :: x3(:)
          REAL*4, ALLOCATABLE                               :: y1(:)
          REAL*4, ALLOCATABLE                               :: y2(:)
      !      REAL*4, ALLOCATABLE                               :: y3(:)
          REAL*4, ALLOCATABLE                               :: zzi(:)
      
          REAL*4, ALLOCATABLE                               :: wu(:)
          REAL*4, ALLOCATABLE                               :: lu(:)
          REAL*4, ALLOCATABLE                               :: lm(:)
          REAL*4, ALLOCATABLE                               :: hhu(:)
          REAL*4, ALLOCATABLE                               :: wm(:)
      
          REAL*4, ALLOCATABLE				   :: xinfl_1(:,:)
          REAL*4, ALLOCATABLE                :: yinfl_1(:,:)
          REAL*4, ALLOCATABLE                :: a_1(:,:)
          REAL*4, ALLOCATABLE                :: b_1(:,:)
          REAL*4, ALLOCATABLE                :: width_1(:)
          INTEGER*2, ALLOCATABLE             :: iub_1(:)
          REAL*4, ALLOCATABLE				   :: length_1(:)
          REAL*4, ALLOCATABLE                :: hu_1(:)
          REAL*4, ALLOCATABLE                :: zi_1(:)
          
          INTEGER*4							tmp_noub
            REAL*4											:: wm_max
            REAL*4											:: lm_max
      
          INTEGER*4                                         :: m
          INTEGER*4                                         :: idir                           !WSU 05-13-08
          INTEGER*4                                         :: il                           !WSU 06-19-08
          INTEGER*4                                         :: ir                           !WSU 06-19-08
          INTEGER*4                                         :: nninfl
          INTEGER*4                                         :: k
          INTEGER*4                                         :: kk
      
          INTEGER*4, ALLOCATABLE                            :: iiub(:)
          !INTEGER*4, ALLOCATABLE                            :: mm(:) !Naum never used
          INTEGER*4                                         :: mcb
          INTEGER*4                                         :: maxUB
      
          INTEGER*4                                         :: jtest
      
          CHARACTER (LEN=64)                                :: filename
          CHARACTER (LEN=15)                                :: adum

          REAL*4                                            :: tmpmxlub
          REAL*4                                            :: tmpwxub
          INTEGER*4                                         :: tmpnumub

          INTEGER*4                                         :: iiii
      
        1 FORMAT(1X,a)
       11 FORMAT(a12)
      101 FORMAT(1X,/)
       31 FORMAT(4(i4,2X),3(f15.5,2X,f15.5,2X),2(i4,2X),4(f15.5,2X),  &
                 i2,2X,2(f13.2,2X))
       33 FORMAT(2(i2,2X),f15.5)
       41 FORMAT(i4,2X,2(f15.5,2X))
       32 FORMAT(2(i4,2X),3(f13.2,2X))
       34 FORMAT(2(i4,2X))
       30 FORMAT(2(i4,2X),2X,3(f13.2),2X)
      
      !     FILENAME READS AND OPEN STATEMENTS
      
          IF (iproc == 0) THEN
            WRITE(*,101)
            WRITE(*,101)
            WRITE(*,101)
            WRITE(*,1) 'Program arcgen'
            WRITE(*,1) '**************'
            WRITE(*,101)
          ENDIF
      
          OPEN(14,FILE='UBLOC.dat',STATUS='old')
          OPEN(55,FILE='maxv.out',ACCESS='APPEND',STATUS='old') 
      
      !     filename(1:10) = 'arcpt.out.'
      !     WRITE(filename(11:),intform) iproc
      !     OPEN(15,FILE=filename,FORM='formatted',STATUS='unknown')
      
          IF (iproc == 0) THEN
      !        READ(14,*)ati
            READ(14,*)humn, huvar
            READ(14,*)wumn, wuvar
            READ(14,*)lumn, luvar
            READ(14,*)ang                                            !WSU 06-19-08
            READ(14,*)yymn, yyvar
      !        READ(14,*)aamn, aavar                                   !WSU 6/19/08
            
          ENDIF
          CLOSE(14)
      
      
          CALL MPI_Bcast(humn,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
          adum = 'humn-mpi'
          CALL CHECK(adum,ierr)
      
          CALL MPI_Bcast(huvar,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
          adum = 'huvar-mpi'
          CALL CHECK(adum,ierr)
      
          CALL MPI_Bcast(wumn,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
          adum = 'wumn-mpi'
          CALL CHECK(adum,ierr)
      
          CALL MPI_Bcast(wuvar,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
          adum = 'wuvar-mpi'
          CALL CHECK(adum,ierr)
      
          CALL MPI_Bcast(lumn,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
          adum = 'lumn-mpi'
          CALL CHECK(adum,ierr)
      
          CALL MPI_Bcast(luvar,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
          adum = 'luvar-mpi'
          CALL CHECK(adum,ierr)
      
      !      CALL MPI_Bcast(xxmn,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      !      adum = 'xxmn-mpi'
      !      CALL CHECK(adum,ierr)
      
      !      CALL MPI_Bcast(xxvar,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      !      adum = 'xxvar-mpi'
      !      CALL CHECK(adum,ierr)
      
          CALL MPI_Bcast(yymn,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
          adum = 'xxmn-mpi'
          CALL CHECK(adum,ierr)
      
          CALL MPI_Bcast(yyvar,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
          adum = 'xxvar-mpi'
          CALL CHECK(adum,ierr)
      
          CALL MPI_Bcast(ang,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
          adum = 'ang-mpi'
          CALL CHECK(adum,ierr)
      
      !      CALL MPI_Bcast(aamn,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      !      adum = 'aamn-mpi'
      !      CALL CHECK(adum,ierr)
      
      !      CALL MPI_Bcast(aavar,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      !      adum = 'aavar-mpi'
      !      CALL CHECK(adum,ierr)
      
    !       mcb = 9500                                                        !WSU 4/4/08
      
    !       ALLOCATE (xxx(mcb), STAT= IERR)
    !       ALLOCATE (yyy(mcb), STAT= IERR)
    !       ALLOCATE (ap(mcb), STAT= IERR)
    !       ALLOCATE (bp(mcb), STAT= IERR)
    !       ALLOCATE (x1(mcb), STAT= IERR)
    !       ALLOCATE (x2(mcb), STAT= IERR)
    !   !      ALLOCATE (x3(mcb), STAT= IERR)
    !       ALLOCATE (y1(mcb), STAT= IERR)
    !       ALLOCATE (y2(mcb), STAT= IERR)
    !   !      ALLOCATE (y3(mcb), STAT= IERR)
    !       ALLOCATE (wu(mcb), STAT= IERR)
    !       ALLOCATE (lu(mcb), STAT= IERR)
    !       ALLOCATE (lm(mcb), STAT= IERR)
    !       ALLOCATE (hhu(mcb), STAT= IERR)
    !       ALLOCATE (wm(mcb), STAT= IERR)
    !       ALLOCATE (iiub(mcb), STAT= IERR)
    !       !ALLOCATE (mm(mcb), STAT= IERR) !Naum never used
    !       ALLOCATE (zzi(mcb), STAT= IERR)
      
      !     ALLOCATE (xxx(mcb,ncbar_partition), STAT= IERR)
      !     ALLOCATE (yyy(mcb,ncbar_partition), STAT= IERR)
      !     ALLOCATE (ap(mcb,ncbar_partition), STAT= IERR)
      !     ALLOCATE (bp(mcb,ncbar_partition), STAT= IERR)
      !     ALLOCATE (x1(mcb,ncbar_partition), STAT= IERR)
      !     ALLOCATE (x2(mcb,ncbar_partition), STAT= IERR)
      ! !      ALLOCATE (x3(mcb), STAT= IERR)
      !     ALLOCATE (y1(mcb,ncbar_partition), STAT= IERR)
      !     ALLOCATE (y2(mcb,ncbar_partition), STAT= IERR)
      ! !      ALLOCATE (y3(mcb), STAT= IERR)
      !     ALLOCATE (wu(mcb,ncbar_partition), STAT= IERR)
      !     ALLOCATE (lu(mcb,ncbar_partition), STAT= IERR)
      !     ALLOCATE (lm(mcb,ncbar_partition), STAT= IERR)
      !     ALLOCATE (hhu(mcb,ncbar_partition), STAT= IERR)
      !     ALLOCATE (wm(mcb,ncbar_partition), STAT= IERR)
      !     ALLOCATE (iiub(mcb,ncbar_partition), STAT= IERR)
      !   !ALLOCATE (mm(mcb), STAT= IERR) !Naum never used
      !     ALLOCATE (zzi(mcb,ncbar_partition), STAT= IERR)
      
      
      ! This is set to one for testing!
      
          jtest = 0
          
      ! Assuming 1 compound bar per processor

      tmpwxub = 0.0
      tmpmxlub = 0.0
      tmpnumub = 0
      DO i_partition = 1, ncbar_process(iproc+1)
        !   kk = iproc * ncbar_partition + i_partition
        IF (iproc<nc_determine) THEN
            kk=iproc*(ncbar_partition+1)+i_partition
        ELSE
            ! kk1=myY*ncbar_partition+i_nb-1
            kk= nc_determine*(ncbar_partition+1)+&
                  (iproc-nc_determine)*ncbar_partition+i_partition
        END IF
          IF (kk <= nocb) THEN
              mcb = 9500                                                        !WSU 4/4/08
              iiii = kk-1
              ALLOCATE (xxx(mcb), STAT= IERR)
              ALLOCATE (yyy(mcb), STAT= IERR)
              ALLOCATE (ap(mcb), STAT= IERR)
              ALLOCATE (bp(mcb), STAT= IERR)
              ALLOCATE (x1(mcb), STAT= IERR)
              ALLOCATE (x2(mcb), STAT= IERR)
          !      ALLOCATE (x3(mcb), STAT= IERR)
              ALLOCATE (y1(mcb), STAT= IERR)
              ALLOCATE (y2(mcb), STAT= IERR)
          !      ALLOCATE (y3(mcb), STAT= IERR)
              ALLOCATE (wu(mcb), STAT= IERR)
              ALLOCATE (lu(mcb), STAT= IERR)
              ALLOCATE (lm(mcb), STAT= IERR)
              ALLOCATE (hhu(mcb), STAT= IERR)
              ALLOCATE (wm(mcb), STAT= IERR)
              ALLOCATE (iiub(mcb), STAT= IERR)
              !ALLOCATE (mm(mcb), STAT= IERR) !Naum never used
              ALLOCATE (zzi(mcb), STAT= IERR)
      
              zzi = 0.0
              xxi = 0     
              wa = (2*wcb(kk))/5.0                                              !WSU, 4/4/08                                                     
              wb = wcb(kk)                                                      !WSU, 4/4/08
              wc = (2*wcb(kk))/3.0                                              !WSU, 4/4/08
              l1 = lcb(kk)/3.0                                                  !WSU, 4/4/08
              l2 = lcb(kk)/3.0                                                  !WSU, 4/4/08
              l3 = lcb(kk)/3.0                                                  !WSU, 4/4/08
              xi = wcb(kk)/2.0                                                  !WSU, 4/4/08 
              yi = 0.0-50                                                       !WSU, 4/4/08
              zstart = z(kk)-humn/1.2                                             !WSU, 4/4/08
              lc = lcb(kk)                                                      !WSU, 4/4/08
              zx = z(kk)-5.0                                                  !WSU, 4/4/08
      !      xi = xc                                                           !WSU, 4/4/08
      !      yi = yc-lumn/3.0                                                 !WSU, 4/4/08
      !      zstart = z(kk)-humn/2.0                                           !WSU, 4/4/08
      !      lc = yc+lcb(kk)                                                  !WSU, 4/4/08
      !      lx = xc-((wcb(kk)/2.0)+wumn)                                     !WSU, 4/4/08
      !      rx = xc+((wcb(kk)/2.0)+wumn)                                     !WSU, 4/4/08
      !      zx = zstart-(hcb(kk)+humn/2.0)                                   !WSU, 4/4/08 
              angcb1 = ATAN(l1/(wb/2.0-wa/2.0))                                  !WSU 06-19-08
              delango = ABS((ang-angcb1)/10.0)                                 !WSU 06-19-08
              delangs = ABS((ang-angcb1)/7.0)                                 !WSU 06-19-08
          
              l = 0
              m = 1
              zzi(m) = zstart 
      ! Loop with each layer
      
              k = 0
          9 IF (zzi(m) > zx)THEN
              k = k+1
              yii = yi
              lx = xi-wa/2.0                                                  !WSU, 4/4/08                                                
              rx = xi+wa/2.0                                                  !WSU, 4/4/08
              ninfl = 2                                                       !WSU, 4/4/08
              xxx(m) = xi
              yyy(m) = yi
              idir = 0                                                         !WSU 05-13-08
              il = 0                                                             !WSU 06-19-08
              ir = 0                                                            !WSU 06-19-08
       10   IF (xxx(m) >= lx)THEN
             il = il+1                                                       !WSU 06-19-08 
              IF(xxx(m) == xi)THEN
                iiub(m) = 1
              ELSE
                iiub(m) = 2
              END IF
       21     CONTINUE
              CALL norgen(w,wumn,wuvar,is1,is2,is3,jtest)
              IF (w <= 0.0) GOTO 21
       22     CONTINUE
              CALL norgen(l,lumn,luvar,is1,is2,is3,jtest)
              IF (l <= 0.0) GOTO 22
       19     CONTINUE
              CALL norgen(h,humn,huvar,is1,is2,is3,jtest)
              IF (h <= 0.0) GOTO 19
              wu(m) = w !70.0
              lu(m) = l !300.0
              hhu(m) = h !1.5
              IF (m == 1)THEN
                GOTO 98
              ENDIF
              wm(m) = wu(m)
              IF (wm(m) < wm(m-1))THEN
                wm(m) = wm(m-1)
              ENDIF
                lm(m) = lu(m)
              IF (lm(m) < lm(m-1))THEN
                lm(m) = lm(m-1)
              ENDIF
       98     CONTINUE
              x1(m) = xxx(m)
              y1(m) = yyy(m)
              ango = ang-il*delango                                          !WSU 06-19-08
              angs = ang-il*delangs                                           !WSU 06-19-08
              atii = -(TAN(angs))                                             !WSU 06-19-08
              IF(atii >= 0.0)THEN                                             !WSU 06-19-08
               atii = -(atii)                                                 !WSU 06-19-08
              ENDIF                                                           !WSU 06-19-08
              xx = lu(m)/TAN(ango)                                             !WSU 06-19-08 
              IF (idir==1)THEN                                              !WSU 05-13-08
               x2(m) = xxx(m)+xx                                             !WSU 05-13-08
              ELSE                                                          !WSU 05-13-08
               x2(m) = xxx(m)-xx
              ENDIF                                                         !WSU 05-13-08
              y2(m) = yyy(m)+lu(m)                                          !wsu 4/4/08
      !          x3(m) = xxx(m)                                               !wsu 4/4/08
      !          y3(m) = yyy(m)+lu(m)                                         !wsu 4/4/08
              ap(m) = -1/atii                                                !WSU 06-19-08
              bp(m) = y1(m)-ap(m)*x1(m)
              
      !vf          WRITE(1555,*)m,x1(m),y1(m),x2(m),y2(m)  
      !                  x3(m),y3(m),iinf,xxi,ap(m),bp(m),wu(m),lu(m),iiub(m),hhu(m),zzi(m)
      !vf          WRITE(1357,*)ap(m),bp(m)
              zzi(m) = zstart
              m = m+1
              xxx(m) = xxx(m-1)-(wu(m-1)/2.0)
              CALL norgen(yy,yymn,yyvar,is1,is2,is3,jtest)
      !          WRITE(1556,*)m
      !          yy = 8.0
              yyy(m) = yyy(m-1)-yy
              GOTO 10
            ELSE
              xxx(m) = xi+(wu(m-1)/2.0)
              yyy(m) = yi-yy
       50     IF (xxx(m) <= rx)THEN
                ir =ir+1                                                    !WSU 06-19-08
                iiub(m) = 3
       51       CONTINUE
                CALL norgen(w,wumn,wuvar,is1,is2,is3,jtest)
                IF (w <= 0.0) GOTO 51
       52       CONTINUE
                CALL norgen(l,lumn,luvar,is1,is2,is3,jtest)
                IF (l <= 0.0) GOTO 52
      119       CONTINUE
                CALL norgen(h,humn,huvar,is1,is2,is3,jtest)
                IF (h <= 0.0) GOTO 119
                 wu(m) = w !70.0
                lu(m) = l !300.0
                hhu(m) = h !1.5
                IF (m == 1)THEN
                  GOTO 199
                ENDIF
                IF (wm(m) < wm(m-1))THEN
                  wm(m) = wm(m-1)
                ENDIF
                lm(m) = lu(m)
                IF (lm(m) < lm(m-1))THEN
                  lm(m) = lm(m-1)
                ENDIF
      199       CONTINUE
                x1(m) = xxx(m)
                y1(m) = yyy(m)
                ango = ang-ir*delango                                       !WSU 06-19-08
                angs = ang-ir*delangs                                      !WSU 06-19-08
                atii = TAN(angs)                                           !WSU 06-19-08
                IF (atii <= 0.0)THEN                                       !WSU 06-19-08
                 atii = -(atii)                                            !WSU 06-19-08
                ENDIF
                xx = lu(m)/TAN(ango)                                        !WSU 05-13-08
                IF (idir==1)THEN                                           !WSU 05-13-08
                 x2(m) = xxx(m)-xx                                         !WSU 05-13-08
                ELSE                                                       !WSU 05-13-08
                 x2(m) = xxx(m)+xx
                ENDIF             
                y2(m) = yyy(m)+lu(m)                                       !WSU, 4/4/08 
      !           x3(m) = x(m)                                                  !WSU, 4/4/08
      !           y3(m) = y(m)+lu(m)                                            !WSU, 4/4/08
                ap(m) = -1/atii
                bp(m) = y1(m)-ap(m)*x1(m)
      !vf            WRITE(1555,*)m,x1(m),y1(m),x2(m),y2(m)
      !vf            WRITE(1357,*)ap(m),bp(m)  
      !                x3(m),y3(m),iinf,xxi,ap(m),bp(m),wu(m),lu(m),iiub(m),hhu(m),zzi(m)
                zzi(m) = zstart
                m = m+1
                xxx(m) = xxx(m-1)+(wu(m-1)/2.0)
                CALL norgen(yy,yymn,yyvar,is1,is2,is3,jtest)
      !vf            WRITE(155,*)'yy=',yy
      !           yy = 8.0
                yyy(m) = yyy(m-1)-yy
                GOTO 50
              ELSE
                IF (yyy(m) < (lc-lumn/2.0)) THEN                          !WSU, 4/4/08                               
                  IF (yi < lc/4.0)THEN                                     !WSU, 5/13/08
                    yyy(m) = yi+lu(m-1)/2.0
                  ELSE
                    yyy(m) = yi+lu(m-1)/3.5
                  END IF
                  IF (yyy(m) <=l1) THEN                                    !WSU, 4/4/08
                xa = ((wb/2.0-wa/2.0)/l1)*yyy(m)                               !WSU, 4/4/08
                IF ((zzi(m)-z(kk)) < hcb(kk))THEN                              !WSU, 4/4/08
                 lx = xi-(wa/2.0+xa)                                          !WSU, 4/4/08
                 rx = xi+(wa/2.0+xa)                                          !WSU, 4/4/08
                 GO TO 44                                                     !WSU, 4/4/08
                ELSE                                                         !WSU, 4/4/08
                 lx = xi-(wa/2.0+xa)-50                                       !WSU, 4/4/08
                 rx = xi+(wa/2.0+xa)+50                                       !WSU, 4/4/08
                 GO TO 44                                                     !WSU, 4/4/08
                ENDIF                                                        !WSU, 4/4/08
              ELSE                                                         !WSU, 4/4/08
                idir = 1                                                     !WSU 06-19-08
                IF (yyy(m) <=(l1+l2)) THEN                                    !WSU, 4/4/08
                 ya = (l1+l2)-yyy(m)                                           !WSU, 4/4/08
                 xa = ((wb/2.0-wc/2.0)/l2)*ya                                !WSU, 4/4/08
                 IF ((zzi(m)-z(kk)) < hcb(kk))THEN                            !WSU, 4/4/08
                  lx = xi-(wa/2.0+xa)                                          !WSU, 4/4/08
                  rx = xi+(wa/2.0+xa)                                          !WSU, 4/4/08
                  GO TO 44                                                     !WSU, 4/4/08
                 ELSE                                                         !WSU, 4/4/08
                  lx = xi-(wc/2.0+xa)-50                                       !WSU, 4/4/08
                  rx = xi+(wc/2.0+xa)+50                                       !WSU, 4/4/08
                  GO TO 44                                                    !WSU, 4/4/08 
                 ENDIF                                                      !WSU, 4/4/08
                ELSE                                                        !WSU, 4/4/08
                 ya = yyy(m)-(l1+l2)                                           !WSU, 4/4/08
                         xa = ((wc/2.0)/l3)*ya                              !WSU, 4/4/08
                IF ((zzi(m)-z(kk)) < hcb(kk))THEN                             !WSU, 4/4/08
                 lx = xi-(wa/2.0-xa)                                          !WSU, 4/4/08
                 rx = xi+(wa/2.0-xa)                                          !WSU, 4/4/08
                 GO TO 44                                                     !WSU, 4/4/08
                ELSE                                                         !WSU, 4/4/08
                 lx = xi-(wc/2.0-xa)-50                                       !WSU, 4/4/08
                 rx = xi+(wc/2.0-xa)+50                                       !WSU, 4/4/08
                 GO TO 44                                                    !WSU, 4/4/08
               ENDIF                                                       !WSU, 4/4/08
              ENDIF                                                       !WSU, 4/4/08
             ENDIF                                                        !WSU, 4/4/08
         44  continue                                                      !WSU, 4/4/08
                  yi = yyy(m)
                  xxx(m) = xi
                  ir = 0                                                   !WSU 06-19-08
                  il = 0                                                    !WSU 06-19-08
                  GOTO 10
                ELSE
                  GOTO 20
                END IF
              END IF
            END IF
       20   CONTINUE
        
            zzi(m) = zstart-(hhu(m-1)/2.0)                             !WSU, 4/4/08
            zstart = zzi(m)
            yi = yii
            GOTO 9
          END IF
      
          j = 1
          ! noub = m-1
          tmp_noub = m-1
          noub_array(kk)=tmp_noub
          noub = m-1
          IF (tmp_noub > tmpnumub) THEN
            tmpnumub = tmp_noub
          ENDIF

          IF (wm(m-1) > tmpwxub) THEN
            tmpwxub = wm(m-1)
          ENDIF

          IF (lm(m-1) > tmpmxlub) THEN
            tmpmxlub = lm(m-1)
          ENDIF
          WRITE(*,*) 'The number of unit bars after calculation ', kk, tmp_noub
          na =2                                                      !WSU 06-19-08
          ! ALLOCATE (xinfl(noub,na), STAT=IERR)
          ! ALLOCATE (yinfl(noub,na), STAT=IERR)
          ! ALLOCATE (a(noub,na), STAT=IERR)
          ! ALLOCATE (b(noub,na), STAT=IERR)
          ! ALLOCATE (width(noub), STAT=IERR)
          ! ALLOCATE (iub(noub), STAT=IERR)
          ! ALLOCATE (length(noub), STAT=IERR)
          ! ALLOCATE (hu(noub), STAT=IERR)
          ! ALLOCATE (zi(noub), STAT=IERR)
          WRITE(intform(3:3),'(I1)') icount_digits(iiii)
          filename(1:11) = 'length.out.'
          WRITE(filename(12:),intform) kk-1
          ! IF (ibin == 1) THEN
          OPEN(15,FILE=filename,FORM='unformatted',STATUS='unknown')
         
            ! OPEN(15,FILE=filename,FORM='formatted',STATUS='unknown')
    
          ALLOCATE (xinfl_1(noub_array(kk),na),STAT=IERR)
          ALLOCATE (yinfl_1(noub_array(kk),na),STAT=IERR)
          ALLOCATE (a_1(noub_array(kk),na),STAT=IERR)
          ALLOCATE (b_1(noub_array(kk),na),STAT=IERR)
          ALLOCATE (width_1(noub_array(kk)),STAT=IERR)
          ALLOCATE (length_1(noub_array(kk)),STAT=IERR)
          ALLOCATE (hu_1(noub_array(kk)),STAT=IERR)
          ALLOCATE (zi_1(noub_array(kk)),STAT=IERR)
          ALLOCATE (iub_1(noub_array(kk)),STAT=IERR)
      
          DO n = 1,noub_array(kk) 
      !       xinfl(n,1) =  x1(n)
      !       xinfl(n,2) =  x2(n)
      ! !       xinfl(j,n,3) =  x3(n)                                                !WSU, 4/4/08
      !       yinfl(n,1) =  y1(n)
      !       yinfl(n,2) =  y2(n)
      ! !       yinfl(j,n,3) =  y3(n)                                                 !WSU, 4/4/08
      !       a(n,1) = ap(n)
      !       b(n,1) = bp(n)
      !       width(n) = wu(n)
      !       length(n) = lu(n)
      !       iub(n) = iiub(n)
      !       hu(n) = hhu(n)
      !       zi(n) = zzi(n)

            xinfl_1(n,1)=x1(n)
            xinfl_1(n,2)=x2(n)
            yinfl_1(n,1)=y1(n)
            yinfl_1(n,2)=y2(n)
            a_1(n,1)=ap(n)
            b_1(n,1)=bp(n)
            width_1(n)=wu(n)
            length_1(n)=lu(n)
            hu_1(n)=hhu(n)
            zi_1(n)=zzi(n)
            iub_1(n)=iiub(n)
          ENDDO
          ! IF (ibin == 1) THEN
            WRITE(15) xinfl_1(1:noub_array(kk),1:na)
            WRITE(15) yinfl_1(1:noub_array(kk),1:na)
            WRITE(15) a_1(1:noub_array(kk),1:na)
            WRITE(15) b_1(1:noub_array(kk),1:na)
            WRITE(15) length_1(1:noub_array(kk))
            WRITE(15) width_1(1:noub_array(kk))
            WRITE(15) hu_1(1:noub_array(kk))
            WRITE(15) zi_1(1:noub_array(kk))
            WRITE(15) iub_1(1:noub_array(kk))
          ! ELSE
            ! WRITE(15) xinfl_1(1:noub_array(kk),1:na)
            ! WRITE(15) yinfl_1(1:noub_array(kk),1:na)
            ! WRITE(15) a_1(1:noub_array(kk),1:na)
            ! WRITE(15) b_1(1:noub_array(kk),1:na)
            ! WRITE(15) length_1(1:noub_array(kk))
            ! WRITE(15) width_1(1:noub_array(kk))
            ! WRITE(15) hu_1(1:noub_array(kk))
            ! WRITE(15) zi_1(1:noub_array(kk))
            ! WRITE(15) iub_1(1:noub_array(kk))
          ! ENDIF
          CLOSE(15)
      ! Can deallocate a bunch of arrays here.
      
      
      !      IF (kk /= 1) THEN
      !        IF (mm(kk) < mm(kk-1))THEN
      !           mm(kk) = mm(kk-1)
      !        ENDIF
      !      ENDIF
              DEALLOCATE(xxx, STAT=IERR)
              DEALLOCATE(yyy, STAT=IERR)
      
              DEALLOCATE(lm, STAT=IERR)
              DEALLOCATE(wm, STAT=IERR)
      
              DEALLOCATE(x1, STAT=IERR)
              DEALLOCATE(x2, STAT=IERR)
              DEALLOCATE(y1, STAT=IERR)
              DEALLOCATE(y2, STAT=IERR)
              DEALLOCATE(ap, STAT=IERR)
              DEALLOCATE(bp, STAT=IERR)
              DEALLOCATE(wu, STAT=IERR)
              DEALLOCATE(lu, STAT=IERR)
              DEALLOCATE(iiub, STAT=IERR)
              DEALLOCATE(hhu, STAT=IERR)
              DEALLOCATE(zzi, STAT=IERR)
  
              DEALLOCATE(xinfl_1,STAT=IERR)
              DEALLOCATE(yinfl_1,STAT=IERR)
              DEALLOCATE(a_1,STAT=IERR)
              DEALLOCATE(b_1,STAT=IERR)
              DEALLOCATE(width_1,STAT=IERR)
              DEALLOCATE(length_1,STAT=IERR)
              DEALLOCATE(hu_1,STAT=IERR)
              DEALLOCATE(zi_1,STAT=IERR)
              DEALLOCATE(iub_1,STAT=IERR)
          ENDIF
          ENDDO
      ! DEALLOCATE
      
      ! This is probably not needed...check to see if used in merge
      ! It is needed; otherwise, the following calculations will be 
      ! wrong.
        !   DEALLOCATE(xxx, STAT=IERR)
        !   DEALLOCATE(yyy, STAT=IERR)
  
        !   ! DEALLOCATE(lm, STAT=IERR)
        !   ! DEALLOCATE(wm, STAT=IERR)
  
        !   DEALLOCATE(x1, STAT=IERR)
        !   DEALLOCATE(x2, STAT=IERR)
        !   DEALLOCATE(y1, STAT=IERR)
        !   DEALLOCATE(y2, STAT=IERR)
        !   DEALLOCATE(ap, STAT=IERR)
        !   DEALLOCATE(bp, STAT=IERR)
        !   DEALLOCATE(wu, STAT=IERR)
        !   DEALLOCATE(lu, STAT=IERR)
        !   DEALLOCATE(iiub, STAT=IERR)
        !   DEALLOCATE(hhu, STAT=IERR)
        !   DEALLOCATE(zzi, STAT=IERR)

        !   CALL MPI_AllReduce(noub,maxUB,1,MPI_INTEGER,MPI_MAX, & !Naum
        !                      MPI_COMM_WORLD,ierr)                         
        !   adum = 'maxUB-mpi'
        !   CALL CHECK(adum,ierr)

          CALL MPI_AllReduce(tmpnumub,maxUB,1,MPI_INTEGER,MPI_MAX, & !Naum
                             MPI_COMM_WORLD,ierr)                         
          adum = 'maxUB-mpi'
          CALL CHECK(adum,ierr)
      
        !   CALL MPI_AllReduce(wm(m-1),wxub,1,MPI_REAL,MPI_MAX, &
        !                     MPI_COMM_WORLD,ierr)
        !   adum = 'wxub-mpi'
        !   CALL CHECK(adum,ierr)

          CALL MPI_AllReduce(tmpwxub,wxub,1,MPI_REAL,MPI_MAX, &
                            MPI_COMM_WORLD,ierr)
          adum = 'wxub-mpi'
          CALL CHECK(adum,ierr)
      
        !  CALL MPI_AllReduce(lm(m-1),mxlub,1,MPI_REAL,MPI_MAX, &
        !                     MPI_COMM_WORLD,ierr)
        !  adum = 'mxlub-mpi'
          CALL MPI_AllReduce(tmpmxlub,mxlub,1,MPI_REAL,MPI_MAX, &
          MPI_COMM_WORLD,ierr)
          adum = 'mxlub-mpi'

         WRITE(*,*) 'The maximum of noub is ', maxUB
         WRITE(*,*) 'The max wub is ', tmpwxub, wxub
         WRITE(*,*) 'The max lub is ', tmpmxlub, mxlub
         CALL CHECK(adum,ierr)     
        !  WRITE(55,*)maxUB,wm(m-1),lm(m-1)
         WRITE(55,*)maxUB,tmpwxub,tmpmxlub
          CLOSE(9)
          CLOSE(55) !Naum
          DEALLOCATE(lm, STAT=IERR)
          DEALLOCATE(wm, STAT=IERR)
          RETURN
        END SUBROUTINE arcgen

  !***********************************************************************************************************************************
  ! Code converted using TO_F90 by Alan Miller
  ! Date: 2007-07-28  Time: 10:03:26
  
  ! Program TRGENA:  This program generates a series
  ! of arcs with specified centerpoints and inflection
  ! points represents the centerlines of n-number of unit bars/compound bars.
  ! Modified on August 28,06.
  !   Input files
  !     Output from arcgen (arc info)
  
  !   Output files:
  !     1. Directly use as an input file in pixdat.f
  !     2. width, length of UB, will be used in creating plane equations (UBCH.f).
  !     iub = 1 (middle); =2(left); =3(right)
  !  June 1: add compound bar
  !********c*********c*********c*********c*********c*********c*********c**
        SUBROUTINE trgena
  
            USE prcmpi
            USE bar
      
            INCLUDE 'mpif.h'
      
            REAL*4, ALLOCATABLE                               :: xx(:,:)
      
            REAL*4, ALLOCATABLE                               :: xc(:,:)
            REAL*4, ALLOCATABLE                               :: yc(:,:)
            REAL*4, ALLOCATABLE                               :: r(:,:)
      
            REAL*4, ALLOCATABLE                               :: alpha(:,:)
            REAL*4, ALLOCATABLE                               :: larc(:,:)
            REAL*4, ALLOCATABLE                               :: lsum(:,:) !Naum OK here
            REAL*4, ALLOCATABLE                               :: zm(:)     !05/13/08
            
            REAL*4, ALLOCATABLE                               :: minx(:)
            REAL*4, ALLOCATABLE                               :: miny(:)
      
            REAL*4                                            :: acon
            REAL*4                                            :: bcon
            REAL*4                                            :: xmid
            REAL*4                                            :: ymid
            REAL*4                                            :: arg
            REAL*4                                            :: apbis
            REAL*4                                            :: bpbi
            REAL*4                                            :: dcon
      
            INTEGER*4, ALLOCATABLE                            :: narc(:)
            INTEGER*4, ALLOCATABLE                            :: iinf(:,:)
            INTEGER*4                                         :: ninfl
      
            INTEGER*4                                         :: j
            INTEGER*4                                         :: na
            INTEGER*4                                         :: n
          INTEGER*4      									:: i_partition 
          INTEGER*4											:: n_partition
          
          !Arrays for reading the data from the output files of arcgen subroutine
          REAL*4, ALLOCATABLE 								:: xinfl_1(:,:)
          REAL*4, ALLOCATABLE               :: yinfl_1(:,:)
          REAL*4, ALLOCATABLE								:: a_1(:,:)
          REAL*4, ALLOCATABLE								:: b_1(:,:)
          REAL*4, ALLOCATABLE								:: width_1(:)
          REAL*4, ALLOCATABLE								:: length_1(:)
          REAL*4, ALLOCATABLE								:: hu_1(:)
          REAL*4, ALLOCATABLE								:: zi_1(:)
          INTEGER*2, ALLOCATABLE							:: iub_1(:)
      
            CHARACTER (LEN=64)                                :: filename
          CHARACTER (LEN=64)                                :: filename1
            CHARACTER (LEN=15)                                :: adum
      
      
          1 FORMAT(1X,a)
         10 FORMAT(a12)
         15 FORMAT(1X,4(f15.5,2X),f15.3,2X,f15.9,2X,2(f15.5,2X),i2,2X, f15.5)
         20 FORMAT(1X,2(f10.5,2X),i2,2X,f10.5)
         25 FORMAT(1X,f10.5)
         31 FORMAT(4(i4,2X),3(f15.5,2X,f15.5,2X),2(i4,2X),4(f15.5,2X),  &
                   i2,2X,2(f13.2,2X))
         32 FORMAT(4(i4,2X),2X,3(f15.5,2X))
         33 FORMAT(2(i4,2X))
         34 FORMAT(2(i4,2X),f15.5)
        101 FORMAT(1X,/)
      
            IF (iproc == 0) THEN
              WRITE(*,101)
              WRITE(*,1) 'Arc Generation Program TRGENA: '
              WRITE(*,1) '********************************* '
              WRITE(*,101)
            ENDIF
            
          !Start the loop for each component bar
          
            DO i_partition = 0, ncbar_process(iproc+1)-1
                ! iiii = iproc * ncbar_partition + i_partition
                IF (iproc<nc_determine) THEN
                    iiii=iproc*(ncbar_partition+1)+i_partition
                ELSE
                    ! kk1=myY*ncbar_partition+i_nb-1
                    iiii= nc_determine*(ncbar_partition+1)+&
                          (iproc-nc_determine)*ncbar_partition+i_partition
                END IF
                IF (iiii < nocb) THEN
                  n_partition = iiii
                  na = 2
                  WRITE(*,*) 'The number of unit bars in trgena is ', iiii, noub_array(n_partition+1)                                                       !WSU 06-19-08
                  ALLOCATE (xc(noub_array(n_partition+1),na), STAT=IERR)
                  ALLOCATE (yc(noub_array(n_partition+1),na), STAT=IERR)
                  ALLOCATE (xx(noub_array(n_partition+1),na), STAT=IERR)
                  ALLOCATE (lsum(noub_array(n_partition+1),na), STAT=IERR)
                  ALLOCATE (narc(noub_array(n_partition+1)), STAT=IERR)
                  ALLOCATE (r(noub_array(n_partition+1),na), STAT=IERR)
                  ALLOCATE (alpha(noub_array(n_partition+1),na), STAT=IERR)
                  ALLOCATE (larc(noub_array(n_partition+1),na), STAT=IERR)
                  ALLOCATE (iinf(noub_array(n_partition+1),na), STAT=IERR)
            !vf      ALLOCATE (iinf(ncbar,noub,na), STAT=IERR)
                  ALLOCATE (zm(noub_array(n_partition+1)), STAT=IERR)                             !WSU, 5/13/08
                  ALLOCATE (minx(noub_array(n_partition+1)), STAT=IERR)
                  ALLOCATE (miny(noub_array(n_partition+1)), STAT=IERR)
      
                  WRITE(*,*) 'Do the calculation of trgena for component bar: ',n_partition
      
                  WRITE(intform(3:3),'(I1)') icount_digits(iiii)
                  filename(1:11) = 'trgena.out.'
                  WRITE(filename(12:),intform) n_partition !iproc
                  IF (ibin == 1) THEN
                      OPEN(8,FILE=filename,FORM='unformatted',STATUS='unknown')
                  ELSE
                      OPEN(8,FILE=filename,FORM='formatted',STATUS='unknown')
                  ENDIF
              !Allocate the memory space for arrays from arcgen subroutine
                  ALLOCATE (xinfl_1(noub_array(n_partition+1),na),STAT=IERR)
                  ALLOCATE (yinfl_1(noub_array(n_partition+1),na),STAT=IERR)
                  ALLOCATE (a_1(noub_array(n_partition+1),na),STAT=IERR)
                  ALLOCATE (b_1(noub_array(n_partition+1),na),STAT=IERR)
                  ALLOCATE (width_1(noub_array(n_partition+1)),STAT=IERR)
                  ALLOCATE (length_1(noub_array(n_partition+1)),STAT=IERR)
                  ALLOCATE (hu_1(noub_array(n_partition+1)),STAT=IERR)
                  ALLOCATE (zi_1(noub_array(n_partition+1)),STAT=IERR)
                  ALLOCATE (iub_1(noub_array(n_partition+1)),STAT=IERR)
              
              !Read the data generated by arcgen from the file
                  filename1(1:11) = 'length.out.'
                  WRITE(filename1(12:),intform) n_partition
                  OPEN(15,FILE=filename1,FORM='unformatted',STATUS='unknown')
                  READ(15)xinfl_1(1:noub_array(n_partition+1),1:na)
                  READ(15)yinfl_1(1:noub_array(n_partition+1),1:na)
                  READ(15)a_1(1:noub_array(n_partition+1),1:na)
                  READ(15)b_1(1:noub_array(n_partition+1),1:na)
                  READ(15)length_1(1:noub_array(n_partition+1))
                  READ(15)width_1(1:noub_array(n_partition+1))
                  READ(15)hu_1(1:noub_array(n_partition+1))
                  READ(15)zi_1(1:noub_array(n_partition+1))
                  READ(15)iub_1(1:noub_array(n_partition+1))
                  CLOSE(15)
    
                  iinf = 0
              
                  xc = 0.0
                  yc = 0.0
                  xx = 0.0
                  lsum = 0.0
                  narc = 0
                  r = 0.0
                  alpha = 0.0
                  larc = 0.0
                  iinf = 0
                  zm = 0.0
                  minx = 0.0
                  miny = 0.0
            
      
        ! Loop over each CB
      
                  j = 1
                  ninfl = 2                                                  !WSU, 4/4/08
      
          !  iiii = iproc
      
            
                !WRITE(intform(3:3),'(I1)') icount_digits(iiii)
                !filename(1:11) = 'trgena.out.'
                !WRITE(filename(12:),intform) n_partition !iproc
                !IF (ibin == 1) THEN
                  !	OPEN(8,FILE=filename,FORM='unformatted',STATUS='unknown')
                !ELSE
                  !	OPEN(8,FILE=filename,FORM='formatted',STATUS='unknown')
                !ENDIF
      
                  DO  n = 1,noub_array(n_partition+1)
                      lsum (n,1) = 0.0
                      narc(n) = ninfl - 1
                      zm(n) = zi_1(n)+hu_1(n)                                 !WSU, 6/19/08       
                      minx(n)=0.0
                      miny(n)=0.0
            !     Loop over each arc
                
                      DO  i = 1, narc(n)
      
            ! Find centerpoint coordinates
      
                        
                        IF (iinf(n,i) == 1) THEN
                    
            ! Calculate centerpoint for vertical perpendicular
                    
                            xmid = (xinfl_1(n,i+1)-xinfl_1(n,i))/2.0 + xinfl_1(n,i)
                            ymid = (yinfl_1(n,i+1)-yinfl_1(n,i))/2.0 + yinfl_1(n,i)
                            acon = (yinfl_1(n,i+1)-yinfl_1(n,i))/(xinfl_1(n,i+1) &
                                  -xinfl_1(n,i))
                            apbis = -1.0/acon
                            bpbis = ymid - apbis*xmid
                            xc(n,i) = xx(n,i)
                            yc(n,i) = apbis*xc(n,i) + bpbis
                
                        ELSE IF (ABS(yinfl_1(n,i+1)-yinfl_1(n,i)) < 0.000001) THEN
                    
            ! Calculate centerpoint for horizontal connector
            ! (vertical perpendicular bisector)
                    
                            xmid = (xinfl_1(n,i+1)-xinfl_1(n,i))/2.0 + xinfl_1(n,i)
                            ymid = yinfl_1(n,i)
                            xc(n,i) = xmid
                            yc(n,i) = a_1(n,i)*xc(n,i) + b_1(n,i)
                  
                        ELSE
                      
            ! Calculate centerpoint for both lines well defined
                    
                            xmid = (xinfl_1(n,i+1)-xinfl_1(n,i))/2.0 + xinfl_1(n,i)
                            ymid = (yinfl_1(n,i+1)-yinfl_1(n,i))/2.0 + yinfl_1(n,i)
                            acon = (yinfl_1(n,i+1)-yinfl_1(n,i))/(xinfl_1(n,i+1) &
                                -xinfl_1(n,i))
                            apbis = -1.0/acon
                            bpbis = ymid - apbis*xmid
                            xc(n,i) = (bpbis-b_1(n,i))/(a_1(n,i)-apbis)
                            yc(n,i) = a_1(n,i)*xc(n,i) + b_1(n,i)
                    
                        ENDIF
                  
            ! Find arc radius
      
                        
                        r(n,i) = xydist(xc(n,i),yc(n,i),xinfl_1(n,i), &
                                  yinfl_1(n,i))
                  
            ! Find central angle of arc
      
                        
                        dcon = xydist(xinfl_1(n,i),yinfl_1(n,i),xinfl_1(n,i+1),  &
                            yinfl_1(n,i+1))
                        arg = 1.0 - 0.5*((dcon/r(n,i))**2.0)
                        alpha(n,i) = ACOS(arg)
      
            ! Find length of arc segment, and add to sum of lengths
                  
                        larc(n,i) = alpha(n,i)*r(n,i)
                        lsum(n,i+1) = lsum(n,i) + larc(n,i)
                  
            ! Find slope and intercept of next perpendicular
                  
                        diff = ABS(xinfl_1(n,i+1) - xc(n,i))
                          
                        IF (diff < 0.000001) THEN
                            iinf(n,i+1) = 1
                            xx(n,i+1) = xinfl_1(n,i+1)
                            a_1(n,i+1) = 0.0
                            b_1(n,i+1) = 0.0
                        ELSE
                            iinf(n,i+1) = 0
                            xx(n,i+1) = 0
                            a_1(n,i+1) = (yinfl_1(n,i+1)-yc(n,i))/ &
                                    (xinfl_1(n,i+1)-xc(n,i))
                            b_1(n,i+1) = yinfl_1(n,i+1) - a_1(n,i+1)*xinfl_1(n,i+1)
                      END IF
            !vf     WRITE(888,*)xc(j,n,i), yc(j,n,i),r(j,n,i),alpha(j,n,i), larc(j,n,i), lsum(j,n,i)
                    END DO
                
            ! ! Write out parameters for the final inflection point, and the width
                
                  END DO
      
                  IF (ibin == 1) THEN
                      WRITE(8) xc(1:noub_array(n_partition+1),1:na)
                      WRITE(8) yc(1:noub_array(n_partition+1),1:na)
                    ! WRITE(8) xinfl(1:noub,1:na,1:ncbar_partition)
                    ! WRITE(8) yinfl(1:noub,1:na,1:ncbar_partition)
      
                      WRITE(8) xinfl_1(1:noub_array(n_partition+1),1:na)
                      WRITE(8) yinfl_1(1:noub_array(n_partition+1),1:na)
      
                      WRITE(8) r(1:noub_array(n_partition+1),1:na)
                      WRITE(8) alpha(1:noub_array(n_partition+1),1:na)
                      WRITE(8) larc(1:noub_array(n_partition+1),1:na)
                      WRITE(8) lsum(1:noub_array(n_partition+1),1:na)
                      WRITE(8) iinf(1:noub_array(n_partition+1),1:na)
                      WRITE(8) xx(1:noub_array(n_partition+1),1:na)
                      WRITE(8) minx(1:noub_array(n_partition+1))
                      WRITE(8) miny(1:noub_array(n_partition+1))                                !WSU, 6/19/08
            ! !       WRITE(8) 0.0                                              !WSU, 6/19/08
                    ! WRITE(8) zi(1:noub,1:ncbar_partition)                                   !WSU, 6/19/08
                    ! WRITE(8) width(1:noub,1:ncbar_partition)                                !WSU, 6/19/08
                    ! WRITE(8) length(1:noub,1:ncbar_partition)                               !WSU, 6/19/08
                      WRITE(8) zi_1(1:noub_array(n_partition+1))                                   !WSU, 6/19/08
                      WRITE(8) width_1(1:noub_array(n_partition+1))                                !WSU, 6/19/08
                      WRITE(8) length_1(1:noub_array(n_partition+1))   
      
                      WRITE(8) zm(1:noub_array(n_partition+1))                                   !WSU, 6/19/08
                      WRITE(8) hu_1(1:noub_array(n_partition+1))
            ! ! Write out the height
      
                    ! WRITE(8) hu(1:noub,1:ncbar_partition)
                  ELSE
                      WRITE(8,*) xc(1:noub_array(n_partition+1),1:na)
                      WRITE(8,*) yc(1:noub_array(n_partition+1),1:na)
                    ! WRITE(8,*) xinfl(1:noub,1:na,1:ncbar_partition)
                    ! WRITE(8,*) yinfl(1:noub,1:na,1:ncbar_partition)
      
                      WRITE(8) xinfl_1(1:noub_array(n_partition+1),1:na)
                      WRITE(8) yinfl_1(1:noub_array(n_partition+1),1:na)
      
                      WRITE(8,*) r(1:noub_array(n_partition+1),1:na)
                      WRITE(8,*) alpha(1:noub_array(n_partition+1),1:na)
                      WRITE(8,*) larc(1:noub_array(n_partition+1),1:na)
                      WRITE(8,*) lsum(1:noub_array(n_partition+1),1:na)
                      WRITE(8,*) iinf(1:noub_array(n_partition+1),1:na)
                      WRITE(8,*) xx(1:noub_array(n_partition+1),1:na)
                      WRITE(8,*) minx(1:noub_array(n_partition+1))
                      WRITE(8,*) miny(1:noub_array(n_partition+1))                                !WSU, 6/19/08
            ! !       WRITE(8,*) 0.0                                              !WSU, 6/19/08
                    ! WRITE(8,*) zi(1:noub,1:ncbar_partition)                                   !WSU, 6/19/08
                    ! WRITE(8,*) width(1:noub,1:ncbar_partition)                                !WSU, 6/19/08
                    ! WRITE(8,*) length(1:noub,1:ncbar_partition)                               !WSU, 6/19/08
                      WRITE(8,*) zi_1(1:noub_array(n_partition+1))                                   !WSU, 6/19/08
                      WRITE(8,*) width_1(1:noub_array(n_partition+1))                                !WSU, 6/19/08
                      WRITE(8,*) length_1(1:noub_array(n_partition+1))   
      
                      WRITE(8,*) zm(1:noub_array(n_partition+1))
                      WRITE(8,*) hu_1(1:noub_array(n_partition+1))
                    CLOSE(8)                                   !WSU, 6/19/08
                  ENDIF
      
                  WRITE(*,*) 'Free the memory of the arrays for the component bar ', n_partition
                  DEALLOCATE (xc, STAT=IERR)
                  DEALLOCATE (yc, STAT=IERR)
                  DEALLOCATE (xx, STAT=IERR)
                  !DEALLOCATE (xinfl, STAT=IERR)
                  !DEALLOCATE (yinfl, STAT=IERR)
                  DEALLOCATE (lsum, STAT=IERR)
                  DEALLOCATE (narc, STAT=IERR)
                  DEALLOCATE (r, STAT=IERR)
                  DEALLOCATE (alpha, STAT=IERR)
                  DEALLOCATE (larc, STAT=IERR)
                  DEALLOCATE (iinf, STAT=IERR)
                  DEALLOCATE (zm, STAT=IERR)
                  DEALLOCATE (minx, STAT=IERR)
                  DEALLOCATE (miny, STAT=IERR)
                  WRITE(*,*) 'Finish freeing the memory of the arrays for the component bar ', n_partition
            ! ! Write out the height
      
                    ! WRITE(8,*) hu(1:noub,1:ncbar_partition)
                  ! ENDIF
                  DEALLOCATE(xinfl_1, STAT=IERR)
                  DEALLOCATE(yinfl_1, STAT=IERR)
                  DEALLOCATE(a_1, STAT=IERR)
                  DEALLOCATE(b_1, STAT=IERR)
                  DEALLOCATE(width_1, STAT=IERR)
                  DEALLOCATE(length_1, STAT=IERR)
                  DEALLOCATE(hu_1,STAT=IERR)
                  DEALLOCATE(zi_1,STAT=IERR)
                  DEALLOCATE(iub_1,STAT=IERR)
                  WRITE(*,*) 'Finish freeing the memory of the array from the last subroutine ', n_partition
                  
              ENDIF
            END DO
      
      
      
            RETURN
            END SUBROUTINE trgena


   !***********************************************************************
   
  ! Code converted using TO_F90 by Alan Miller
  ! Date: 2007-07-28  Time: 11:24:32
  
  ! Subroutine UBCH establishes the dataset which describes the
  ! characteristic geometry of a unit bar(middle, left, and right).
  ! created by Guin, September 21,06..creates n number of ubs
  !     iub = 1 (middle); =2(left); =3(right)
  ! 4/24/07   change height, h1(top) and h2(bottom) h2=h1(given)-h(distribution)
  ! may9/07...add bottom boundaries
  !
  !    Output file: Plane equation (used as input file in pixdat.f
  !***********************************************************************
            SUBROUTINE ubch
  
                USE prcmpi
                USE bar
          
                INCLUDE 'mpif.h'
          
                REAL*4                                            :: thi1
                REAL*4                                            :: thi2
                REAL*4                                            :: thi3
                REAL*4                                            :: thi4
                REAL*4                                            :: thi5
                REAL*4                                            :: tript(3,3)
          
                REAL*4, ALLOCATABLE                               :: aaa(:,:)
                REAL*4, ALLOCATABLE                               :: bbb(:,:)
                REAL*4, ALLOCATABLE                               :: ccc(:,:)
                REAL*4, ALLOCATABLE                               :: ddd(:,:)
          
                REAL*4                                            :: l
                REAL*4                                            :: l1
                REAL*4                                            :: l2
                REAL*4                                            :: capl
                REAL*4                                            :: w
          
                REAL*4, ALLOCATABLE                               :: ab(:,:)
                REAL*4, ALLOCATABLE                               :: bb(:,:)
                REAL*4, ALLOCATABLE                               :: cb(:,:)
                REAL*4, ALLOCATABLE                               :: db(:,:)
        
            !Arrays for reading the data from the output files of arcgen subroutine
                REAL*4, ALLOCATABLE 								              :: xinfl_1(:,:)
                REAL*4, ALLOCATABLE                               :: yinfl_1(:,:)
                REAL*4, ALLOCATABLE								                :: a_1(:,:)
                REAL*4, ALLOCATABLE								                :: b_1(:,:)
                REAL*4, ALLOCATABLE								                :: width_1(:)
                REAL*4, ALLOCATABLE								                :: length_1(:)
                REAL*4, ALLOCATABLE								                :: hu_1(:)
                REAL*4, ALLOCATABLE								                :: zi_1(:)
                REAL*4, ALLOCATABLE                               :: zhmax_1(:)
                REAL*4, ALLOCATABLE                               :: zhmin_1(:)
                INTEGER*2, ALLOCATABLE							              :: iub_1(:)
          
                REAL*4                                            :: zh
                REAL*4                                            :: zhmn
                REAL*4                                            :: zhvar
                REAL*4                                            :: th
                REAL*4                                            :: h1
                REAL*4                                            :: h2
                REAL*4                                            :: det
          
                INTEGER*4                                         :: ii
                INTEGER*4                                         :: n
                INTEGER*4                                         :: jj
                INTEGER*4                                         :: j
                INTEGER*4                                         :: nup
                INTEGER*4                                         :: nbp
                INTEGER*4                                         :: nhp
                INTEGER*4                                         :: m
                INTEGER*4                                         :: kk
                INTEGER*4											                    :: i_partition
                INTEGER*4											                    :: n_partition
                INTEGER*4                                         :: n_partition1
                INTEGER*4                                         :: na
              
          
                INTEGER*4,ALLOCATABLE                             :: nop(:)
          
                INTEGER*4                                         :: jtest
          
                CHARACTER (LEN=64)                                :: filename
                CHARACTER (LEN=64)                                :: filename1
                CHARACTER (LEN=64)                                :: filename2
                CHARACTER (LEN=15)								:: adum
          
          !     REAL*4, ALLOCATABLE :: length(:,:),width(:,:),hu(:,:),hi(:,:)
          
              1 FORMAT(1X,a)
             10 FORMAT(a12)
             15 FORMAT(1X,6(f8.6,2X))
             16 FORMAT(1X,7(f8.6,2X))
             20 FORMAT(1X,3(f11.8,3X))
             25 FORMAT(1X,4(f15.7,3X))
             30 FORMAT(1X,a,i1,a,f15.7)
             33 FORMAT(2(i4,2X),f15.5)
             34 FORMAT(2(i4,2X))
             35 FORMAT(1X,9(f12.8,2X))
            101 FORMAT(1X,/)
             32 FORMAT(4(i4,2X),2X,3(f15.5,2X))
          
                IF (iproc == 0) THEN
                  WRITE(*,101)
                  WRITE(*,101)
                  WRITE(*,1) 'Set Characteristic Geometry: ubch'
                  WRITE(*,1) '****************************'
                  WRITE(*,101)
                ENDIF
                
                OPEN(4,FILE='UBPLANE.dat',STATUS='old')
          
          ! Reading input files
          
                IF (iproc == 0) THEN
                  READ(4,*)zhmn,zhvar,th
                ENDIF
                CLOSE(4)
                CALL MPI_Bcast(zhmn,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
                adum = 'zhmn-mpi'
                CALL CHECK(adum,ierr)
        
                CALL MPI_Bcast(zhvar,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
                adum = 'zhvar-mpi'
                CALL CHECK(adum,ierr)
        
                CALL MPI_Bcast(th,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
                adum = 'th-mpi'
                CALL CHECK(adum,ierr)
          ! Starting loop
              !  iiii = iproc
                DO i_partition = 0, ncbar_process(iproc+1)-1
                  ! iiii = iproc * ncbar_partition + i_partition
                  IF (iproc<nc_determine) THEN
                      iiii=iproc*(ncbar_partition+1)+i_partition
                  ELSE
                      ! kk1=myY*ncbar_partition+i_nb-1
                      iiii= nc_determine*(ncbar_partition+1)+&
                            (iproc-nc_determine)*ncbar_partition+i_partition
                  END IF
                    IF (iiii < nocb) THEN
                        WRITE(*,*) 'The ubch from the process ', iiii
                        n_partition = iiii
                        na = 2
                        WRITE(*,*) 'The calculation of ubch for the component bar is ',n_partition
                        !For each component bar
                        ! ALLOCATE (zhmin(noub), STAT=IERR)
                        ! adum = 'zhmin'
                        ! CALL CHECK(adum,ierr)
                  
                        ! ALLOCATE (zhmax(noub), STAT=IERR)
                        ! adum = 'zhmax'
                        ! CALL CHECK(adum,ierr)

                        ALLOCATE (nop(noub_array(n_partition+1)), STAT=IERR)
                        adum = 'nop'
                        CALL CHECK(adum,ierr)
                  
                        ALLOCATE (zhmin_1(noub_array(n_partition+1)), STAT=IERR)
                        adum = 'zhmin_1'
                        CALL CHECK(adum,ierr)
                  
                        ALLOCATE (zhmax_1(noub_array(n_partition+1)), STAT=IERR)
                        adum = 'zhmax_1'
                        CALL CHECK(adum,ierr)
                  
                        ALLOCATE(aaa(noub_array(n_partition+1),22))
                        adum = 'aaa-ubch'
                        CALL CHECK(adum,ierr)
                  
                        ALLOCATE(bbb(noub_array(n_partition+1),22))
                        adum = 'bbb-ubch'
                        CALL CHECK(adum,ierr)
                  
                        ALLOCATE(ccc(noub_array(n_partition+1),22))
                        adum = 'ccc-ubch'
                        CALL CHECK(adum,ierr)
                  
                        ALLOCATE(ddd(noub_array(n_partition+1),22))
                        adum = 'ddd-ubch'
                        CALL CHECK(adum,ierr)
                  
                        ALLOCATE(ab(noub_array(n_partition+1),22))
                        adum = 'ab-ubch'
                        CALL CHECK(adum,ierr)
                  
                        ALLOCATE(bb(noub_array(n_partition+1),22))
                        adum = 'bb-ubch'
                        CALL CHECK(adum,ierr)
                  
                        ALLOCATE(cb(noub_array(n_partition+1),22))
                        adum = 'cb-ubch'
                        CALL CHECK(adum,ierr)
                  
                        ALLOCATE(db(noub_array(n_partition+1),22))
                        adum = 'db-ubch'
                        CALL CHECK(adum,ierr)
        
                        !Read the file including original zi, hu data
                        ALLOCATE (xinfl_1(noub_array(n_partition+1),na),STAT=IERR)
                        ALLOCATE (yinfl_1(noub_array(n_partition+1),na),STAT=IERR)
                        ALLOCATE (a_1(noub_array(n_partition+1),na),STAT=IERR)
                        ALLOCATE (b_1(noub_array(n_partition+1),na),STAT=IERR)
                        ALLOCATE (width_1(noub_array(n_partition+1)),STAT=IERR)
                        ALLOCATE (length_1(noub_array(n_partition+1)),STAT=IERR)
                        ALLOCATE (hu_1(noub_array(n_partition+1)),STAT=IERR)
                        ALLOCATE (zi_1(noub_array(n_partition+1)),STAT=IERR)
                        ALLOCATE (iub_1(noub_array(n_partition+1)),STAT=IERR)
        
                        WRITE(*,*) 'Finishing allocate the memory for the array in ubch subroutine and start&
                         to read in process ',n_partition
        
                    !Read the data generated by arcgen from the file
                        WRITE(intform(3:3),'(I1)') icount_digits(iiii)
                        filename1(1:11) = 'length.out.'
                        WRITE(filename1(12:),intform) n_partition
                        OPEN(15,FILE=filename1,FORM='unformatted',STATUS='unknown')
                        READ(15)xinfl_1(1:noub_array(n_partition+1),1:na)
                        READ(15)yinfl_1(1:noub_array(n_partition+1),1:na)
                        READ(15)a_1(1:noub_array(n_partition+1),1:na)
                        READ(15)b_1(1:noub_array(n_partition+1),1:na)
                        READ(15)length_1(1:noub_array(n_partition+1))
                        READ(15)width_1(1:noub_array(n_partition+1))
                        READ(15)hu_1(1:noub_array(n_partition+1))
                        READ(15)zi_1(1:noub_array(n_partition+1))
                        READ(15)iub_1(1:noub_array(n_partition+1))
                        CLOSE(15)

                        ! WRITE(intform(3:3),'(I1)') icount_digits(iiii)
                        filename2(1:11) = 'zhpara.out.'
                        WRITE(filename2(12:),intform) n_partition
                        OPEN(151,FILE=filename2,FORM='unformatted',STATUS='unknown')
        
        
                        WRITE(*,*) 'Finishing reading the data for the process ',n_partition
        
                        aaa = 0.0
                        bbb = 0.0
                        ccc = 0.0
                        ddd = 0.0
                  
                        ba = 0.0
                        bb = 0.0
                        bc = 0.0
                        bd = 0.0
        
                        ! WRITE(intform(3:3),'(I1)') icount_digits(iiii)
                        filename(1:9) = 'ubch.out.'
                        WRITE(filename(10:),intform) n_partition	!iproc
                        OPEN(8,FILE=filename,FORM='unformatted',STATUS='unknown')
          
          ! ! This is set to one for testing only!
          
                        jtest = 0
          
                        kk = 1
                        DO  n = 1,noub_array(n_partition+1)
                
                      125   CONTINUE
                            CALL norgen(zh,zhmn,zhvar,is1,is2,is3,jtest)
                            IF (zh <= 0.0) GO TO 125
                    ! !           zh = 0.1 
                                zhmin_1(n) = zi_1(n)                                         !WSU, 4/4/08
                                zhmax_1(n) = zhmin_1(n)+hu_1(n)                             !WSU, 4/4/08

                                ! zhmin(n) = zi_1(n)                                         !WSU, 4/4/08
                                ! zhmax(n) = zhmin(n)+hu_1(n)                             !WSU, 4/4/08  
                                   
                                nup = 17
                                nbp = 2
                                nhp = 1
          
                            IF(iub_1(n) == 1)THEN
                                nop(n) = 2
                            ELSE
                                nop(n) = 1
                            END IF
          
              ! ! Calculation of parameters (nup = 17)
                      
                            w = width_1(n)
                            capl = length_1(n)
              ! !vf        WRITE(1567,*)n,capl,w,zi(kk,n),hu(kk,n)
                          l1 = hu_1(n)/TAN(0.436332)
                          l2 = capl-l1
                      
                          r = w/2
                          thi1 = 0.2617994
                          thi2 = 2*thi1
                         thi3 = 3*thi1
                          thi4 = 4*thi1
                          thi5 = 5*thi1
                          l = r-l1
                
                      ! Plane 1(top):
                
                        tript(1,1) = 0.0
                        tript(1,2) = r
                      !      IF (iub(kk,n) == 3)THEN                                   !WSU 4/4/08
                      !      tript(1,3) = zhmax(kk,n)/4.0                                      !WSU 4/4/08
                      !      ELSE                                                      !WSU 4/4/08
                        tript(1,3) = zhmax_1(n)
                      !      END IF                                                    !WSU 4/4/08
                        tript(2,1) = r-l*SIN(thi5)
                        tript(2,2) = r-l*COS(thi5)
                        tript(2,3) = zhmin_1(n)
                        tript(3,1) = l1
                        tript(3,2) = r
                        tript(3,3) = zhmin_1(n)
                        jj = 1
                        CALL det3d(tript,det)
                        ddd(n,jj) = -det
                        CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
                
                      ! Plane 2(top):
          
                        DO  k = 1, 3
                            tript(3,k) = tript(1,k)
                        END DO
                        tript(1,1) = r-r*SIN(thi4)
                        tript(1,2) = r-r*COS(thi4)
                      !       IF (iub(kk,n) == 3)THEN                                          !WSU 4/4/08
                      !       tript(1,3) = zhmax(kk,n)/2.0                                             !WSU 4/4/08
                      !       ELSE                                                              !WSU 4/4/08
                        tript(1,3) = zhmax_1(n)
                      !       END IF                                                            !WSU 4/4/08
                        jj = 2
                        CALL det3d(tript,det)
                        ddd(n,jj) = -det
                        CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
                
                      ! Plane 3(top):
                
                        DO  k = 1, 3
                            tript(3,k) = tript(2,k)
                        END DO
                        tript(2,1) = r-l*SIN(thi3)
                        tript(2,2) = r-l*COS(thi3)
                        tript(2,3) = zhmin_1(n)
                        jj = 3
                        CALL det3d(tript,det)
                        ddd(n,jj) = -det
                        CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
                
                      ! Plane 4(top)
                
                        DO  k = 1, 3
                          tript(3,k) = tript(1,k)
                        END DO
                        tript(1,1) = r-r*SIN(thi2)
                        tript(1,2) = r-r*COS(thi2)
                      !      IF (iub(kk,n) == 3)THEN                                             !WSU 4/4/08
                      !       tript(1,3) = zhmax(kk,n)/1.33                                              !WSU 4/4/08
                      !       ELSE                                                                !WSU 4/4/08
                        tript(1,3) = zhmax_1(n)
                      !       END IF                                                              !WSU 4/4/08
                        jj = 4
                        CALL det3d(tript,det)
                        ddd(n,jj) = -det
                        CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
                
                      ! Plane 5(top)
          
                      DO  k = 1, 3
                          tript(3,k) = tript(2,k)
                      END DO
                      tript(2,1) = r-l*SIN(thi1)
                      tript(2,2) = r-l*COS(thi1)
                      tript(2,3) = zhmin_1(n)
                      jj = 5
                      CALL det3d(tript,det)
                      ddd(n,jj) = -det
                      CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
              
                    ! Plane 6(top)
              
                      DO  k = 1, 3
                          tript(3,k) = tript(1,k)
                      END DO
                      tript(1,1) = r
                      tript(1,2) = 0
                      tript(1,3) = zhmax_1(n)
                      jj = 6
                      CALL det3d(tript,det)
                      ddd(n,jj) = -det
                      CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
              
              
                    ! Plane 7(top)
              
                      DO  k = 1, 3
                          tript(3,k) = tript(2,k)
                      END DO
                      tript(2,1) = r+l*SIN(thi1)
                      tript(2,2) = r-l*COS(thi1)
                      tript(2,3) = zhmin_1(n)
                      jj = 7
                      CALL det3d(tript,det)
                      ddd(n,jj) = -det
                      CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
          
                    ! Plane 8(top)
              
                      DO  k = 1, 3
                        tript(3,k) = tript(1,k)
                      END DO
                      tript(1,1) = r+r*SIN(thi2)
                      tript(1,2) = r-r*COS(thi2)
                    !       IF (iub(kk,n) == 2)THEN                                                     !WSU 4/4/08
                    !       tript(1,3) = zhmax(kk,n)/1.33                                                      !WSU 4/4/08
                    !       ELSE                                                                        !WSU 4/4/08
                      tript(1,3) = zhmax_1(n)
                    !       END IF                                                                      !WSU 4/4/08
                      jj = 8
                      CALL det3d(tript,det)
                      ddd(n,jj) = -det
                      CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
              
              
                    ! Plane 9(top)
              
                      DO  k = 1, 3
                        tript(3,k) = tript(2,k)
                      END DO
                      tript(2,1) = r+l*SIN(thi3)
                      tript(2,2) = r-l*COS(thi3)
                      tript(2,3) = zhmin_1(n)
                      jj = 9
                      CALL det3d(tript,det)
                      ddd(n,jj) = -det
                      CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
              
                    !  Plane 10(top)
              
                      DO  k = 1, 3
                        tript(3,k) = tript(1,k)
                      END DO
                      tript(1,1) = r+r*SIN(thi4)
                      tript(1,2) = r-r*COS(thi4)
                    !       IF (iub(kk,n) == 2)THEN                                                   !WSU 4/4/08
                    !       tript(1,3) = zhmax(kk,n)/2.0                                                     !WSU 4/4/08
                    !       ELSE                                                                      !WSU 4/4/08
                      tript(1,3) = zhmax_1(n)
                    !       END IF                                                                   !WSU 4/4/08
                      jj = 10
                      CALL det3d(tript,det)
                      ddd(n,jj) = -det
                      CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
              
                    !  Plane 11(top)
              
                      DO  k = 1, 3
                        tript(3,k) = tript(2,k)
                      END DO
                      tript(2,1) = r+l*SIN(thi5)
                      tript(2,2) = r-l*COS(thi5)
                      tript(2,3) = zhmin_1(n)
                      jj = 11
                      CALL det3d(tript,det)
                      ddd(n,jj) = -det
                      CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
              
                    ! Plane 12(top)
              
                      DO  k = 1, 3
                        tript(3,k) = tript(1,k)
                      END DO
                      tript(1,1) = w
                      tript(1,2) = r 
                    !      IF (iub(kk,n) == 2)THEN                                          !WSU 4/4/08
                    !       tript(1,3) = zhmax(kk,n)/4.0                                             !WSU 4/4/08
                    !       ELSE                                                              !WSU 4/4/08
                      tript(1,3) = zhmax_1(n)
                    !       END IF                                                            !WSU 4/4/08
                      jj = 12
                      CALL det3d(tript,det)
                      ddd(n,jj) = -det
                      CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
              
                    ! Plane 13(top)
              
                      DO  k = 1, 3
                        tript(3,k) = tript(2,k)
                      END DO
                      tript(2,1) = r+l
                      tript(2,2) = r
                      tript(2,3) = zhmin_1(n)
                      jj = 13
                      CALL det3d(tript,det)
                      ddd(n,jj) = -det
                      CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
              
                    ! Plane 14(top)
              
                      tript(1,1) = r+l*SIN(thi1)
                      tript(1,2) = r-l*COS(thi1)
                      tript(1,3) = zhmin_1(n)
                      tript(3,1) = l1
                      tript(3,2) = r
                      tript(3,3) = zhmin_1(n)
                      jj = 14
                      CALL det3d(tript,det)
                      ddd(n,jj) = -det
                      CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
              
                    ! Plane 15(top)
              
                      DO  k = 1, 3
                        tript(2,k) = tript(3,k)
                      END DO
                      tript(1,1) = 0.0
                      tript(1,2) = r
                    !       IF (iub(kk,n) == 3)THEN                                                         !WSU 4/4/08 
                    !      tript(1,3) = zhmax(kk,n)/4.0                                                           !WSU 4/4/08 
                    !      ELSE                                                                            !WSU 4/4/08
                      tript(1,3) = zhmax_1(n)
                    !      END IF                                                                          !WSU 4/4/08
                      tript(3,1) = 0.0
                      tript(3,2) = capl
                      tript(3,3) = zhmax_1(n)
                      jj = 15
                      CALL det3d(tript,det)
                      ddd(n,jj) = -det
                      CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
              
                    ! Plane 16(top)
              
                      DO  k = 1, 3
                        tript(1,k) = tript(2,k)
                      END DO
                      tript(2,1) = r+l
                      tript(2,2) = r
                      tript(2,3) = zhmin_1(n)
                      jj = 16
                      CALL det3d(tript,det)
                      ddd(n,jj) = -det
                      CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
              
                    ! Plane 17(top)
              
                      DO  k = 1, 3
                        tript(1,k) = tript(2,k)
                      END DO
                      tript(2,1) = w
                      tript(2,2) = r
                    !       IF (iub(kk,n) == 2)THEN                                            !WSU 4/4/08
                    !       tript(2,3) = zhmax(kk,n)/4.0                                              !WSU 4/4/08
                    !       ELSE                                                               !WSU 4/4/08 
                      tript(2,3) = zhmax_1(n)
                    !       END IF                                                             !WSU 4/4/08
                      tript(3,1) = w
                      tript(3,2) = capl
                      tript(3,3) = zhmax_1(n)
                      jj = 17
                      CALL det3d(tript,det)
                      ddd(n,jj) = -det
                      CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
              
                    ! Plane 18(bottom) (nbp = 2)
              
                      tript(1,1) = r
                      tript(1,2) = 0.0
                      tript(1,3) = zhmax_1(n)
                      IF (iub_1(n) == 3)THEN
                        tript(2,1) = 0.0
                      ELSE
                        tript(2,1) = r
                      END IF
                      tript(2,2) = capl
                      tript(2,3) = zhmax_1(n)
                      tript(3,1) = 0.0
                      IF (iub_1(n) == 3)THEN
                        tript(3,2) = r
                        tript(3,3) = zhmax_1(n)/4.0
                      ELSE
                        tript(3,2) = capl
                        tript(3,3) = zhmax_1(n)
                      END IF
                      jj = 18
                      CALL det3d(tript,det)
                      ddd(n,jj) = -det
                      CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
              
                    ! Plane 19(bottom)
              
                      IF(iub_1(n) == 1)THEN
                        DO  k = 1, 3
                        tript(3,k) = tript(2,k)
                        END DO
                      ELSE
                        IF(iub_1(n) == 2)THEN
                        tript(3,1) = w
                        tript(3,2) = r
                        tript(3,3) = zhmax_1(n)/4.0
                        ELSE
                        tript(3,1) = r
                        tript(3,2) = capl
                        tript(3,3) = zhmax_1(n)
                        END IF
                      END IF
                      tript(2,1) = w
                      tript(2,2) = capl
                      tript(2,3) = zhmax_1(n)
                      jj = 19
                      CALL det3d(tript,det)
                      ddd(n,jj) = -det
                      CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
              
                    ! Horizontal plane (nhp=1)
              
                      tript(1,1) = 0.0
                      tript(1,2) = 0.0
                      tript(1,3) = zhmin_1(n)+zh                                !WSU 4/4/08
                      tript(2,1) = w
                      tript(2,2) = capl
                      tript(2,3) = zhmin_1(n)+zh                                 !WSU 4/4/08
                      tript(3,1) = 0.0
                      tript(3,2) = capl
                      tript(3,3) = zhmin_1(n)+zh                                 !WSU 4/4/08
                      jj = 20
                    !        IF (iub(kk,n) == 1)THEN
                    !          jj = 22
                    !        ELSE
                    !          jj = 21
                    !        END IF
                      CALL det3d(tript,det)
                      ddd(n,jj) = -det
                      CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
              
                    ! Oblique plane(21) parallel to left unit bar's bottom plane (nop = 1 or 2)
              
                      tript(1,1) = r
                      tript(1,2) = 0.0
                      tript(1,3) = zhmin_1(n)                                      !WSU 5/13/08
                      IF (iub_1(n) == 3)THEN
                        tript(2,1) = w
                      ELSE
                        tript(2,1) = r
                      END IF
                      tript(2,2) = capl
                      IF (iub_1(n) == 3)THEN
                        tript(2,3) = zhmax_1(n)                                    !WSU 5/13/08           
                        tript(3,1) = r
                      ELSE
                        tript(2,3) = zhmin_1(n)                                    !WSU 5/13/08
                        tript(3,1) = 0.0
                      END IF
                      tript(3,2) = capl
                      IF (iub_1(n) == 3)THEN
                        tript(3,3) = zhmin_1(n)                                    !WSU 5/13/08
                      ELSE
                        tript(3,3) = zhmax_1(n)                                    !WSU 5/13/08
                      END IF
                      jj = 21
                      CALL det3d(tript,det)
                      ddd(n,jj) = -det
                      CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
              
                    ! Oblique plane(21) parallel to left right bar's bottom plane
              
                      IF (iub_1(n) /= 1)THEN
                        GO TO 99
                      END IF
                      DO  k = 1, 3
                        tript(3,k) = tript(2,k)
                      END DO
                      tript(2,1) = w
                      tript(2,2) = capl
                      tript(2,3) = zhmax_1(n)                                      !WSU 5/13/08
                      jj = 22
                      CALL det3d(tript,det)
                      ddd(n,jj) = -det
                      CALL findabc(tript,aaa(n,jj),bbb(n,jj),ccc(n,jj))
              
                      99     CONTINUE
              
                    ! Bottom planes
                    ! Plane 1a(parallel):  (nup = 17)
              
                      tript(1,1) = 0.0
                      tript(1,2) = r
                    !       IF (iub(kk,n) == 3)THEN                                                      !WSU 4/4/08
                    !       tript(1,3) = (zhmax(kk,n)+th)/4.0                                                    !WSU 4/4/08
                    !       ELSE                                                                          !WSU 4/4/08
                      tript(1,3) = (zhmax_1(n)+th)
                    !       END IF                                                                        !WSU 4/4/08
                      tript(2,1) = r-l*SIN(thi5)
                      tript(2,2) = r-l*COS(thi5)
                      tript(2,3) = (zhmin_1(n)+th)
                      tript(3,1) = l1
                      tript(3,2) = r
                      tript(3,3) = (zhmin_1(n)+th)
                      ii = 1
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
                    ! Plane 2a(top):
              
                      DO  k = 1, 3
                        tript(3,k) = tript(1,k)
                      END DO
                      tript(1,1) = r-r*SIN(thi4)
                      tript(1,2) = r-r*COS(thi4)
                    !       IF (iub(kk,n) == 3)THEN                                                            !WSU 4/4/08
                    !      tript(1,3) = (zhmax(kk,n)+th)/2.0                                                         !WSU 4/4/08
                    !      ELSE                                                                               !WSU 4/4/08
                      tript(1,3) = (zhmax_1(n)+th)
                    !      END IF                                                                              !WSU 4/4/08
                      ii = 2
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
              
                    ! Plane 3a(top):
              
                      DO  k = 1, 3
                        tript(3,k) = tript(2,k)
                      END DO
                      tript(2,1) = r-l*SIN(thi3)
                      tript(2,2) = r-l*COS(thi3)
                      tript(2,3) = zhmin_1(n)+th
                      ii = 3
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
              
                    ! Plane 4a(top)
              
                      DO  k = 1, 3
                        tript(3,k) = tript(1,k)
                      END DO
                      tript(1,1) = r-r*SIN(thi2)
                      tript(1,2) = r-r*COS(thi2)
                    !      IF (iub(kk,n) == 3)THEN                                                          !WSU 4/4/08
                    !       tript(1,3) = (zhmax(kk,n)+th)/1.33                                                      !WSU 4/4/08
                    !       ELSE                                                                             !WSU 4/4/08
                      tript(1,3) = (zhmax_1(n)+th)
                    !       END IF                                                                           !WSU 4/4/08
                      ii = 4
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
                    ! Plane 5a(top)
              
                      DO  k = 1, 3
                        tript(3,k) = tript(2,k)
                      END DO
                      tript(2,1) = r-l*SIN(thi1)
                      tript(2,2) = r-l*COS(thi1)
                      tript(2,3) = zhmin_1(n)+th
                      ii = 5
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
              
                    ! Plane 6a(top)
              
                      DO  k = 1, 3
                        tript(3,k) = tript(1,k)
                      END DO
                      tript(1,1) = r
                      tript(1,2) = 0
                      tript(1,3) = zhmax_1(n)+th
                      ii = 6
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
                    ! Plane 7a(top)
              
                      DO  k = 1, 3
                        tript(3,k) = tript(2,k)
                      END DO
                      tript(2,1) = r+l*SIN(thi1)
                      tript(2,2) = r-l*COS(thi1)
                      tript(2,3) = zhmin_1(n)+th
                      ii = 7
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
              
                    ! Plane 8a(top)
              
                      DO  k = 1, 3
                        tript(3,k) = tript(1,k)
                      END DO
                      tript(1,1) = r+r*SIN(thi2)
                      tript(1,2) = r-r*COS(thi2)
                    !      IF (iub(kk,n) == 2)THEN                                                      !WSU 4/4/08
                    !       tript(1,3) = (zhmax(kk,n)+th)/1.33                                                    !WSU 4/4/08
                    !       ELSE                                                                          !WSU 4/4/08
                      tript(1,3) = zhmax_1(n)+th
                    !       END IF                                                                        !WSU 4/4/08
                      ii = 8
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
                    ! Plane 9a(top)
              
                      DO  k = 1, 3
                        tript(3,k) = tript(2,k)
                      END DO
                      tript(2,1) = r+l*SIN(thi3)
                      tript(2,2) = r-l*COS(thi3)
                      tript(2,3) = zhmin_1(n)+th
                      ii = 9
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
              
                    ! Plane 10a(top)
              
                      DO  k = 1, 3
                        tript(3,k) = tript(1,k)
                      END DO
                      tript(1,1) = r+r*SIN(thi4)
                      tript(1,2) = r-r*COS(thi4)
                    !       IF (iub(kk,n) == 2)THEN                                               !WSU 4/4/08
                    !       tript(1,3) = (zhmax(kk,n)+th)/2.0                                            !WSU 4/4/08   
                    !       ELSE                                                                  !WSU 4/4/08
                      tript(1,3) = zhmax_1(n)+th
                    !       END IF                                                                !WSU 4/4/08
                      ii = 10
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
                    ! Plane 11a(top)
              
                      DO  k = 1, 3
                        tript(3,k) = tript(2,k)
                      END DO
                      tript(2,1) = r+l*SIN(thi5)
                      tript(2,2) = r-l*COS(thi5)
                      tript(2,3) = zhmin_1(n)+th
                      ii = 11
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
                    ! Plane 12a(top)
              
                      DO  k = 1, 3
                        tript(3,k) = tript(1,k)
                      END DO
                      tript(1,1) = w
                      tript(1,2) = r
                    !      IF (iub(kk,n) == 2)THEN                                         !WSU 4/4/08
                    !       tript(1,3) = (zhmax(kk,n)+th)/4.0                                        !WSU 4/4/08
                    !       ELSE                                                             !WSU 4/4/08
                      tript(1,3) = zhmax_1(n)+th
                    !       END IF                                                              !WSU 4/4/08
                      ii = 12
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
                    ! Plane 13a(top)
              
                      DO  k = 1, 3
                        tript(3,k) = tript(2,k)
                      END DO
                      tript(2,1) = r+l
                      tript(2,2) = r
                      tript(2,3) = zhmin_1(n)+th
                      ii = 13
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
                    ! Plane 14a(top)
              
                      tript(1,1) = r+l*SIN(thi1)
                      tript(1,2) = r-l*COS(thi1)
                      tript(1,3) = zhmin_1(n)+th
                      tript(3,1) = l1
                      tript(3,2) = r
                      tript(3,3) = zhmin_1(n)+th
                      ii = 14
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
                    ! Plane 15a(top)
              
                      DO  k = 1, 3
                        tript(2,k) = tript(3,k)
                      END DO
                      tript(1,1) = 0.0
                      tript(1,2) = r
                    !      IF (iub(kk,n) == 3)THEN                                        !WSU 4/4/08
                    !      tript(1,3) = (zhmax(kk,n)+th)/4.0                                      !WSU 4/4/08
                    !      ELSE                                                            !WSU 4/4/08
                      tript(1,3) = zhmax_1(n)+th
                    !      END IF                                                          !WSU 4/4/08
                      tript(3,1) = 0.0
                      tript(3,2) = capl
                      tript(3,3) = zhmax_1(n)+th
                      ii = 15
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
                    ! Plane 16a(top)
              
                      DO  k = 1, 3
                        tript(1,k) = tript(2,k)
                      END DO
                      tript(2,1) = r+l
                      tript(2,2) = r
                      tript(2,3) = zhmin_1(n)+th
                      ii = 16
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
                    ! Plane 17a(top)
              
                      DO  k = 1, 3
                        tript(1,k) = tript(2,k)
                      END DO
                      tript(2,1) = w
                      tript(2,2) = r
                    !       IF (iub(kk,n) == 2)THEN                                             !WSU 4/4/08
                    !       tript(2,3) = (zhmax(kk,n)+th)/4.0                                           !WSU 4/4/08
                    !       ELSE                                                                 !WSU 4/4/08
                      tript(2,3) = zhmax_1(n)+th
                    !       END IF                                                               !WSU 4/4/08
                      tript(3,1) = w
                      tript(3,2) = capl
                      tript(3,3) = zhmax_1(n)+th
                      ii = 17
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
                    ! Plane 18a(bottom) (nbp = 2)
              
                      tript(1,1) = r
                      tript(1,2) = 0.0
                      tript(1,3) = zhmax_1(n)-th
                      IF (iub_1(n) == 3)THEN
                        tript(2,1) = 0.0
                      ELSE
                        tript(2,1) = r
                      END IF
                      tript(2,2) = capl
                      tript(2,3) = zhmax_1(n)-th
                      tript(3,1) = 0.0
                      IF (iub_1(n) == 3)THEN
                        tript(3,2) = r
                        tript(3,3) = zhmax_1(n)-th/4.0
                      ELSE
                        tript(3,2) = capl
                        tript(3,3) = zhmax_1(n)-th
                      END IF
                      ii = 18
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
                    ! Plane 19a(bottom)
              
                      IF(iub_1(n) == 1)THEN
                        DO  k = 1, 3
                        tript(3,k) = tript(2,k)
                        END DO
                      ELSE
                        IF(iub_1(n) == 2)THEN
                            tript(3,1) = w
                            tript(3,2) = r
                            tript(3,3) = (zhmax_1(n)-th)/4.0
                        ELSE
                            tript(3,1) = r
                            tript(3,2) = capl
                            tript(3,3) = zhmax_1(n)-th
                        END IF
                      END IF
                      tript(2,1) = w
                      tript(2,2) = capl
                      tript(2,3) = zhmax_1(n)-th
                      ii = 19
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
                    ! Horizontal plane (nhp = 1)
              
                      tript(1,1) = 0.0
                      tript(1,2) = 0.0
                      tript(1,3) = zhmin_1(n)+zh+th                                        !WSU 4/4/08
                      tript(2,1) = w
                      tript(2,2) = capl
                      tript(2,3) = zhmin_1(n)+zh+th                                          !WSU 4/4/08
                      tript(3,1) = 0.0
                      tript(3,2) = capl
                      tript(3,3) = zhmin_1(n)+zh+th                                         !WSU 4/4/08
                    !       IF (iub(kk,n) == 1)THEN
                    !         ii = 22
                    !       ELSE
                    !         ii = 21
                    !       END IF
                      ii = 20
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
              
                    ! Oblique plane(20a) parallel to left unit bar's bottom plane (nop = 1 or 2)
              
                      tript(1,1) = r
                      tript(1,2) = 0.0
                      tript(1,3) = zhmin_1(n)+th                                             !WSU 5/13/08
                      IF (iub_1(n) == 3)THEN
                        tript(2,1) = w
                      ELSE
                        tript(2,1) = r
                      END IF
                      tript(2,2) = capl
                      IF (iub_1(n) == 3)THEN
                        tript(2,3) = (zhmax_1(n)+th)                                           !WSU 5/13/08
                        tript(3,1) = r
                      ELSE
                        tript(2,3) = zhmin_1(n)+th                                           !WSU 5/13/08
                        tript(3,1) = 0.0
                      END IF
                      tript(3,2) = capl
                      IF (iub_1(n) == 3)THEN
                        tript(3,3) = zhmin_1(n)+th                                            !WSU 5/13/08
                      ELSE
                        tript(3,3) = (zhmax_1(n)+th)                                          !WSU 5/13/08
                      END IF
                      ii = 21
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
              
                    ! Oblique plane(21a) parallel to left right bar's bottom plane
              
                      IF (iub_1(n) /= 1)THEN
                        GO TO 199
                      END IF
                      DO  k = 1, 3
                        tript(3,k) = tript(2,k)
                      END DO
                      tript(2,1) = w
                      tript(2,2) = capl
                      tript(2,3) = zhmax_1(n)+th                                              !WSU 5/13/08
                      ii = 22
                      CALL det3d(tript,det)
                      db(n,ii) = -det
                      CALL findabc(tript,ab(n,ii),bb(n,ii),cb(n,ii))
                      199   CONTINUE
                  ! !vf        DO  j = 1, jj
                  ! !vf        WRITE(998,*) aaa(n,j), bbb(n,j), ccc(n,j),ddd(n,j)
                  ! !vf        END DO
                  ! !        write(8,*)nup,nbp,nop,nhp
                  ! !vf        DO  j = 1, ii
                  ! !vf        WRITE(998,*) ab(n,j), bb(n,j), cb(n,j), db(n,j)
                  ! !vf        END DO
          
                    END DO
          
              ! ! Write results
          
              ! ! Check on what needs to be written out here.
          
                      WRITE(8)nop(1:noub_array(n_partition+1))
          
              ! ! First write nup = 17 planes
                      WRITE(8)aaa(1:noub_array(n_partition+1),1:nup)
                      WRITE(8)bbb(1:noub_array(n_partition+1),1:nup)
                      WRITE(8)ccc(1:noub_array(n_partition+1),1:nup)
                      WRITE(8)ddd(1:noub_array(n_partition+1),1:nup)
                    
                      WRITE(8)ab(1:noub_array(n_partition+1),1:nup)
                      WRITE(8)bb(1:noub_array(n_partition+1),1:nup)
                      WRITE(8)cb(1:noub_array(n_partition+1),1:nup)
                      WRITE(8)db(1:noub_array(n_partition+1),1:nup)
                
              ! ! Next nbp = 2 planes 
                      istart = nup + 1
                      WRITE(8)aaa(1:noub_array(n_partition+1),istart:istart+nbp-1)
                      WRITE(8)bbb(1:noub_array(n_partition+1),istart:istart+nbp-1)
                      WRITE(8)ccc(1:noub_array(n_partition+1),istart:istart+nbp-1)
                      WRITE(8)ddd(1:noub_array(n_partition+1),istart:istart+nbp-1)
                
                      WRITE(8)ab(1:noub_array(n_partition+1),istart:istart+nbp-1)
                      WRITE(8)bb(1:noub_array(n_partition+1),istart:istart+nbp-1)
                      WRITE(8)cb(1:noub_array(n_partition+1),istart:istart+nbp-1)
                      WRITE(8)db(1:noub_array(n_partition+1),istart:istart+nbp-1)
          
              ! ! Next nhp = 1 plane
                      istart = istart + nbp
                      WRITE(8)aaa(1:noub_array(n_partition+1),istart:istart+nhp-1)
                      WRITE(8)bbb(1:noub_array(n_partition+1),istart:istart+nhp-1)
                      WRITE(8)ccc(1:noub_array(n_partition+1),istart:istart+nhp-1)
                      WRITE(8)ddd(1:noub_array(n_partition+1),istart:istart+nhp-1)
                
                      WRITE(8)ab(1:noub_array(n_partition+1),istart:istart+nhp-1)
                      WRITE(8)bb(1:noub_array(n_partition+1),istart:istart+nhp-1)
                      WRITE(8)cb(1:noub_array(n_partition+1),istart:istart+nhp-1)
                      WRITE(8)db(1:noub_array(n_partition+1),istart:istart+nhp-1)
                
              ! ! Next nop = 1 or 2 planes
                      istart = istart + nhp
                      WRITE(8)aaa(1:noub_array(n_partition+1),istart:istart+1)
                      WRITE(8)bbb(1:noub_array(n_partition+1),istart:istart+1)
                      WRITE(8)ccc(1:noub_array(n_partition+1),istart:istart+1)
                      WRITE(8)ddd(1:noub_array(n_partition+1),istart:istart+1)
                
                      WRITE(8)ab(1:noub_array(n_partition+1),istart:istart+1)
                      WRITE(8)bb(1:noub_array(n_partition+1),istart:istart+1)
                      WRITE(8)cb(1:noub_array(n_partition+1),istart:istart+1)
                      WRITE(8)db(1:noub_array(n_partition+1),istart:istart+1)
                      CLOSE(8)
        
                      WRITE(151)zhmax_1(1:noub_array(n_partition+1))
                      WRITE(151)zhmin_1(1:noub_array(n_partition+1))
                      CLOSE(151)
        
                      DEALLOCATE(nop, STAT=IERR)
                      adum = 'nop-de-ubch'
                      CALL CHECK(adum,ierr)
        
                      DEALLOCATE(zhmin_1, STAT=IERR)
                      adum = 'zhmin-de-ubch'
                      CALL CHECK(adum,ierr)
        
                      DEALLOCATE(zhmax_1, STAT=IERR)
                      adum = 'zhmax-de-ubch'
                      CALL CHECK(adum,ierr)
                
                      DEALLOCATE(aaa, STAT=IERR)
                      adum = 'aaa-de-ubch'
                      CALL CHECK(adum,ierr)
          
                      DEALLOCATE(bbb, STAT=IERR)
                      adum = 'bbb-de-ubch'
                      CALL CHECK(adum,ierr)
          
                      DEALLOCATE(ccc, STAT=IERR)
                      adum = 'ccc-de-ubch'
                      CALL CHECK(adum,ierr)
          
                      DEALLOCATE(ddd, STAT=IERR)
                      adum = 'ddd-de-ubch'
                      CALL CHECK(adum,ierr)
          
                      DEALLOCATE(ab, STAT=IERR)
                      adum = 'ab-de-ubch'
                      CALL CHECK(adum,ierr)
          
                      DEALLOCATE(bb, STAT=IERR)
                      adum = 'bb-de-ubch'
                      CALL CHECK(adum,ierr)
          
                      DEALLOCATE(cb, STAT=IERR)
                      adum = 'cb-de-ubch'
                      CALL CHECK(adum,ierr)
          
                      DEALLOCATE(db, STAT=IERR)
                      adum = 'db-de-ubch'
                      CALL CHECK(adum,ierr)
        
                      DEALLOCATE(xinfl_1, STAT=IERR)
                      DEALLOCATE(yinfl_1, STAT=IERR)
                      DEALLOCATE(a_1, STAT=IERR)
                      DEALLOCATE(b_1, STAT=IERR)
                      DEALLOCATE(width_1, STAT=IERR)
                      DEALLOCATE(length_1, STAT=IERR)
                      DEALLOCATE(hu_1,STAT=IERR)
                      DEALLOCATE(zi_1,STAT=IERR)
                      DEALLOCATE(iub_1,STAT=IERR)
                  ENDIF
                END DO
          
          
          
                RETURN
                END SUBROUTINE ubch

  !***********************************************************************
  ! Code converted using TO_F90 by Alan Miller
  ! Date: 2007-07-28  Time: 09:53:26
  ! Program TRGENA:  This program generates a description of a trough
  !                  set and creates a file in the format needed for
  !                  input to the program DISCR (which then generates
  !                  a discrete description of the trough).  The
  !                  description generated here consists of a series
  !                  of arcs with specified centerpoints and inflection
  !                  points.
  ! Old version where only one arc can be used
  ! 6/25/07 added cbc arcs..
  !***********************************************************************
                SUBROUTINE trgenao
  
                    USE prcmpi
                    USE bar
              
                    INCLUDE 'mpif.h'
              
                    REAL*4, ALLOCATABLE                               :: xinflc(:)
                    REAL*4, ALLOCATABLE                               :: yinflc(:)
                    REAL*4, ALLOCATABLE                               :: aaa(:)
                    REAL*4, ALLOCATABLE                               :: bbb(:)
                    REAL*4, ALLOCATABLE                               :: xx(:)
              
                    REAL*4                               :: widthcb
                    REAL*4, ALLOCATABLE                               :: xc(:)
                    REAL*4, ALLOCATABLE                               :: yc(:)
                    REAL*4, ALLOCATABLE                               :: r(:)  
              
                    REAL*4, ALLOCATABLE                               :: alpha(:)
                    REAL*4, ALLOCATABLE                               :: larc(:)
                    REAL*4, ALLOCATABLE                               :: lsum(:)
              
                    REAL*4                                            :: acon
                    REAL*4                                            :: bcon
                    REAL*4                                            :: xmid
                    REAL*4                                            :: ymid
                    REAL*4                                            :: arg
                    REAL*4                                            :: apbis
                    REAL*4                                            :: bpbi
                    REAL*4                                            :: dcon
              
                    INTEGER*2, ALLOCATABLE                            :: iinf(:)
                    INTEGER*2                            :: narc
              
                    INTEGER*4                                         :: i_partition
                  INTEGER*4											:: n_partition
                  
                    INTEGER*4                                         :: j
                    INTEGER*4                                         :: k
                    INTEGER*2                                         :: na
              
              
              
                    CHARACTER (LEN=64)                                :: filename
                    CHARACTER (LEN=64)                                :: filename1
                  CHARACTER (LEN=15)								:: adum
                    
              
                  1 FORMAT(1X,a)
                 10 FORMAT(a12)
                 15 FORMAT(1X,4(f15.5,2X),f15.3,2X,f15.9,2X,2(f15.5,2X),i2,2X,f15.5)
                 20 FORMAT(1X,2(f10.5,2X),i2,2X,f10.5)
                 25 FORMAT(1X,f10.5)
                 33 FORMAT(6(f13.2),i4,4(f13.2))
                101 FORMAT(1X,/)
              
                    IF (iproc == 0) THEN
                      WRITE(*,101)
                      WRITE(*,1) 'Trough Generation Program TRGENAO: '
                      WRITE(*,1) '********************************* '
                      WRITE(*,101)
                    ENDIF
              
                  !  iiii = iproc
                    DO i_partition = 0, ncbar_process(iproc+1) - 1
                    !   iiii = iproc * ncbar_partition + i_partition
                      IF (iproc<nc_determine) THEN
                        iiii=iproc*(ncbar_partition+1)+i_partition
                      ELSE
                        ! kk1=myY*ncbar_partition+i_nb-1
                        iiii= nc_determine*(ncbar_partition+1)+&
                              (iproc-nc_determine)*ncbar_partition+i_partition
                      ENDIF
                      IF (iiii < nocb) THEN
                          n_partition = iiii
                      !ENDIF 
                    !END DO
                          na =3
              
                          ALLOCATE (xinflc(na), STAT=IERR)
                          adum = 'xinflc'
                          CALL CHECK(adum,ierr)
                          xinflc = 0.0
              
                          ALLOCATE (yinflc(na), STAT=IERR)
                          adum = 'yinflc'
                          CALL CHECK(adum,ierr)
                          yinflc = 0.0
              
                          ALLOCATE (iinf(na), STAT=IERR)
                          adum = 'iinf'
                          CALL CHECK(adum,ierr)
                          iinf = 0
              
                          ALLOCATE (xx(na), STAT=IERR)
                          adum = 'xx'
                          CALL CHECK(adum,ierr)
                          xx = 0.0
              
                          ALLOCATE (aaa(na), STAT=IERR)
                          adum = 'aaa'
                          CALL CHECK(adum,ierr)
                          aaa = 0.0
              
                          ALLOCATE (bbb(na), STAT=IERR)
                          adum = 'bbb'
                          CALL CHECK(adum,ierr)
                          bbb = 0.0
              
              !      ALLOCATE (widthcb(ncbar), STAT=IERR)
              !      adum = 'widthcb'
              !      CALL CHECK(adum,ierr)
                          widthcb = 0.0
              
                          ALLOCATE (lsum(na), STAT=IERR)
                          adum = 'lsum'
                          CALL CHECK(adum,ierr)
                          lsum = 0.0
              
              !      ALLOCATE (narc(ncbar), STAT=IERR)
              !      adum = 'narc'
              !      CALL CHECK(adum,ierr)
                          narc = 0
              
                          ALLOCATE (xc(na), STAT=IERR)
                          adum = 'xc'
                          CALL CHECK(adum,ierr)
                          xc = 0.0
              
                          ALLOCATE (yc(na), STAT=IERR)    
                          adum = 'yc'
                          CALL CHECK(adum,ierr)
                          yc = 0.0
              
                          ALLOCATE (r(na), STAT=IERR)
                          adum = 'r'
                          CALL CHECK(adum,ierr)
                          r = 0.0
              
                          ALLOCATE (alpha(na), STAT=IERR)    
                          adum = 'alpha'
                          CALL CHECK(adum,ierr)
                          alpha = 0.0
              
                          ALLOCATE (larc(na), STAT=IERR)
                          adum = 'larc'
                          CALL CHECK(adum,ierr)
                          larc = 0
              
              ! Need to read as direct access; need to change write to this file
                          j = 1
                   ! nn = iproc+1
                    !DO i_partition = 1, ncbar_partition
                          !nn = iiii + 1 !iproc * ncbar_partition + (i_partition + 1)
                          !IF (nn <= nocb) THEN
                          WRITE(*,*) 'The calculation for trgenao calculation for the component bar ', n_partition
                          WRITE(intform(3:3),'(I1)') icount_digits(iiii)
                          filename(1:8) = 'cbc.out.'
                          WRITE(filename(9:),intform) n_partition      
                          OPEN(31,FILE=filename,FORM='unformatted',STATUS='old')
                          READ(31) (xinflc(k),k=1,3),(yinflc(k),k=1,3), & !Naum
                                  aaa(1),bbb(1),widthcb  !WSU 6/19/08
                          
                          close(31)
              
                          !WRITE(intform(3:3),'(I1)') icount_digits(iiii)
                          filename1(1:10) = 'cbarc.out.'
                          WRITE(filename1(11:),intform) n_partition !iproc
              
                          iinf = 0
                          xx = 0.0
                      !      aaa(j,1) = -0.5                                !WSU 6/19/08
                      !      bbb(j,1) = 22.5                                !WSU 6/19/08
                      !      widthcb(j) = 20.0                              !WSU 6/19/08
              
              
                          lsum (1) = 0.0
                          narc = 2
                        
                      ! Loop over each arc
                        
                          DO  i = 1, narc
                        
                      ! Find centerpoint coordinates
              
                            IF (iinf(i) == 1) THEN
                            
                      ! Calculate centerpoint for vertical perpendicular***
                            
                            xmid = (xinflc(i+1)-xinflc(i))/2.0 + xinflc(i)
                            ymid = (yinflc(i+1)-yinflc(i))/2.0 + yinflc(i)
                            acon = (yinflc(i+1)-yinflc(i))/(xinflc(i+1)- &
                                xinflc(i))
                            apbis = -1.0/acon
                            bpbis = ymid - apbis*xmid
                            xc(i) = xx(i)
                            yc(i) = apbis*xc(i) + bpbis
              
                            ELSE IF (ABS(yinflc(i+1)-yinflc(i)) < 0.000001) THEN
                          
                      ! Calculate centerpoint for horizontal connector***
                      ! (vertical perpendicular bisector)
                          
                            xmid = (xinflc(i+1)-xinflc(i))/2.0 + xinflc(i)
                            ymid = yinflc(i)
                            xc(i) = xmid
                            yc(i) = aaa(i)*xc(i) + bbb(i)
              
                            ELSE
                            
                      ! Calculate centerpoint for both lines well defined
                            
                            xmid = (xinflc(i+1)-xinflc(i))/2.0 + xinflc(i)
                            ymid = (yinflc(i+1)-yinflc(i))/2.0 + yinflc(i)
                            acon = (yinflc(i+1)-yinflc(i))/(xinflc(i+1)- &
                                xinflc(i))
                            apbis = -1.0/acon
                            bpbis = ymid - apbis*xmid
                            xc(i) = (bpbis-bbb(i))/(aaa(i)-apbis)
                            yc(i) = aaa(i)*xc(i) + bbb(i)
              
                            END IF
                        
                      ! Find arc radius
                        
                            r(i) = xydist(xc(i),yc(i),xinflc(i),yinflc(i))
                        
                      ! Find central angle of arc
                      !
                            dcon = xydist(xinflc(i),yinflc(i),xinflc(i+1), &
                              yinflc(i+1))
                            arg = 1.0 - 0.5*((dcon/r(i))**2.0)
                            alpha(i) = ACOS(arg)
                        
                      ! Find length of arc segment, and add to sum of lengths
                        
                            larc(i) = alpha(i)*r(i)
                            lsum(i+1) = lsum(i) + larc(i)
                        
                      ! Find slope and intercept of next perpendicular
                        
                            diff = ABS(xinflc(i+1) - xc(i))
                            IF (diff < 0.000001) THEN
                            iinf(i+1) = 1
                            xx(i+1) = xinflc(i+1)
                            aaa(i+1) = 0.0
                            bbb(i+1) = 0.0
                            ELSE
                            iinf(i+1) = 0
                            xx(i+1) = 0
                            aaa(i+1) = (yinflc(i+1)-yc(i))/(xinflc(i+1)-xc(i))
                            bbb(i+1) = yinflc(i+1) - aaa(i+1)*xinflc(i+1)
                            END IF
                            
                          END DO
                        
                      !  Write out parameters for the final inflection point, and the trough width
                        
                          !CLOSE(8)
                          
                          !iiii = iproc
                          !iiii = iproc * ncbar_partition + i - 1
                          !n_partition = iiii
              
                          WRITE(*,*) 'Start to output the data for the component bar ',n_partition
                          IF (ibin == 1) THEN
                            OPEN(8,FILE=filename1,FORM='unformatted',STATUS='unknown')
                            WRITE(8)xc
                            WRITE(8)yc
                            WRITE(8)xinflc
                            WRITE(8)yinflc
                            WRITE(8)r
                            WRITE(8)alpha
                            WRITE(8)lsum
                            WRITE(8)widthcb
                          ELSE
                            OPEN(8,FILE=filename1,FORM='formatted',STATUS='unknown')
                            WRITE(8,*)xc
                            WRITE(8,*)yc
                            WRITE(8,*)xinflc
                            WRITE(8,*)yinflc
                            WRITE(8,*)r
                            WRITE(8,*)alpha
                            WRITE(8,*)lsum
                            WRITE(8,*)widthcb
                          ENDIF
                          CLOSE(8)
                          WRITE(*,*) 'Finish outputting the data for the component bar ',n_partition
              
                          DEALLOCATE (xinflc, STAT=IERR)
                          adum = 'xinfl-deal'
                          CALL CHECK(adum,ierr)
                    
                          DEALLOCATE (yinflc, STAT=IERR)
                          adum = 'yinfl-deal'
                          CALL CHECK(adum,ierr)
                    
                          DEALLOCATE (iinf, STAT=IERR)
                          adum = 'iinf-deal'
                          CALL CHECK(adum,ierr)
                    
                          DEALLOCATE (xx, STAT=IERR)
                          adum = 'xx-deal'
                          CALL CHECK(adum,ierr)
                    
                          DEALLOCATE (aaa, STAT=IERR)
                          adum = 'aaa-deal'
                          CALL CHECK(adum,ierr)
                    
                          DEALLOCATE (bbb, STAT=IERR)
                          adum = 'bbb-deal'
                          CALL CHECK(adum,ierr)
                    
                     !     DEALLOCATE (widthcb, STAT=IERR)
                     !     adum = 'widthcb-deal'
                     !     CALL CHECK(adum,ierr)
                    
                          DEALLOCATE (lsum, STAT=IERR)
                          adum = 'lsum-deal'
                          CALL CHECK(adum,ierr)
                    
                     !     DEALLOCATE (narc, STAT=IERR)
                     !     adum = 'narc-deal'
                     !     CALL CHECK(adum,ierr)
                    
                          DEALLOCATE (xc, STAT=IERR)
                          adum = 'xc-deal'
                          CALL CHECK(adum,ierr)
                    
                          DEALLOCATE (yc, STAT=IERR)    
                          adum = 'yc-deal'
                          CALL CHECK(adum,ierr)
                    
                          DEALLOCATE (r, STAT=IERR)
                          adum = 'r-deal'
                          CALL CHECK(adum,ierr)
                    
                          DEALLOCATE (alpha, STAT=IERR)    
                          adum = 'alpha-deal'
                          CALL CHECK(adum,ierr)
                    
                          DEALLOCATE (larc, STAT=IERR)
                          adum = 'larc-deal'
                          CALL CHECK(adum,ierr)
                        ENDIF
                      END DO
              
              
                    RETURN
                    END SUBROUTINE trgenao

!***********************************************************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2007-07-27  Time: 17:15:36

! Program CBCCH establishes the dataset which describes the
! characteristic geometry of a Cross bar channel.  This
! version, does not allow the randomness into the pattern produced.
! Orginally created by, Tim Scheibe, Modified by Guin, Feb 2, 07
!     Total 5 planes
!***********************************************************************

        SUBROUTINE cbcch

            USE prcmpi
            USE bar
        
            INCLUDE 'mpif.h'
        
            REAL*4                                            :: tript(3,3)
            REAL*4                                            :: hmn                            !06-19-08
            REAL*4                                            :: hvar                           !06-19-08
            REAL*4                                            :: ho                           !06-19-08
            REAL*4                                            :: aaa(5)
            REAL*4                                            :: bbb(5)
            REAL*4                                            :: ccc(5)
            REAL*4                                            :: ddd(5)
            REAL*4                                            :: l
            REAL*4                                            :: w
            REAL*4                                            :: det
            REAL*4                                            :: h
        
            INTEGER*4                                         :: i_partition
            INTEGER*4											                    :: n_partition
            
            INTEGER*4                                         :: jj
            INTEGER*4                                         :: nn
        
            CHARACTER (LEN=64)                                :: filename
            CHARACTER (LEN=15)                                :: adum
        
            1 FORMAT(1X,a)
            10 FORMAT(a12)
            15 FORMAT(1X,6(f8.6,2X))
            16 FORMAT(1X,7(f8.6,2X))
            20 FORMAT(1X,3(f11.8,3X))
            25 FORMAT(1X,4(f15.7,3X))
            30 FORMAT(1X,a,i1,a,f15.7)
            35 FORMAT(1X,9(f12.8,2X))
        101 FORMAT(1X,/)
        
            IF (iproc == 0) THEN
                WRITE(*,101)
                WRITE(*,101)
                WRITE(*,1) 'Set Characteristic Geometry: cbcch'
                WRITE(*,1) '****************************'
                WRITE(*,101)
            ENDIF
        
            
            !OPEN(32,FILE='cbclw.out',FORM='unformatted',STATUS='old') !Naum
        
        ! Read input data
        
        ! This file needs to be read as direct accces; also written that way in crtpts
        
            jtest = 0
        
            !  nn = iproc+1
            DO i_partition = 1, ncbar_process(iproc+1)
            ! nn = iproc * ncbar_partition + i_partition
            IF (iproc<nc_determine) THEN
                nn=iproc*(ncbar_partition+1)+i_partition
            ELSE
                ! kk1=myY*ncbar_partition+i_nb-1
                nn= nc_determine*(ncbar_partition+1)+&
                      (iproc-nc_determine)*ncbar_partition+i_partition
            ENDIF
            IF (nn <= nocb) THEN 
                WRITE(*,*) 'The calculation of cbcch is from the component bar ',nn-1
                OPEN(32,FILE='cbclw.out',FORM='unformatted',STATUS='old')
                DO k1 = 1, nn  !Naum
                    READ(32) w,l,ho !Naum      
                enddo !Naum
                close(32)
                hmn = 0.5
                hvar = 0.001
            45 continue
                call norgen(h,hmn,hvar,is1,is2,is3,jtest)
            if (h<=0.0) goto 45
        !      h = 0.7
            aaa = 0.0
            bbb = 0.0
            ccc = 0.0
            ddd = 0.0
        ! Plane 1:
        
            tript(1,1) = 0.0
            tript(1,2) = 0.0
            tript(1,3) = ho                                                    !WSU 4/4/08
            tript(2,1) = tript(1,1)+w/6.0
            tript(2,2) = 0.0
            tript(2,3) = ho+h/1.33                                              !WSU 4/4/08
            tript(3,1) = 0.0
            tript(3,2) = l
            tript(3,3) = ho                                                      !WSU 4/4/08
            jj = 1
            CALL det3d(tript,det)
            ddd(jj) = -det
            CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
        
        ! Plane 2:
        
            DO  k = 1, 3
                tript(1,k) = tript(2,k)
            END DO
            tript(2,1) = tript(1,1)+w/4.0
            tript(2,2) = 0.0
            tript(2,3) = ho+h                                               !WSU 4/4/08
            tript(3,1) = tript(1,1)
            tript(3,2) = l
            tript(3,3) = ho+h/1.33                                          !WSU 4/4/08
            jj = 2
            CALL det3d(tript,det)
            ddd(jj) = -det
            CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
            
        ! Plane 3:
            
            DO  k = 1, 3
                tript(1,k) = tript(2,k)
            END DO
            tript(2,1) = tript(1,1)+w/6.0
            tript(2,2) = 0.0
            tript(2,3) = ho+h                                                !WSU 4/4/08
            tript(3,1) = tript(1,1)
            tript(3,2) = l
            tript(3,3) = ho+h                                                 !WSU 4/4/08
            jj = 3
            CALL det3d(tript,det)
            ddd(jj) = -det
            CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
            
        ! Plane 4
            
            DO  k = 1, 3
                tript(1,k) = tript(2,k)
            END DO
            tript(2,1) = tript(1,1)+w/4.0
            tript(2,2) = 0.0
            tript(2,3) = ho+h/1.33                                            !WSU 4/4/08
            tript(3,1) = tript(1,1)
            tript(3,2) = l
            tript(3,3) = ho+h                                                 !WSU 4/4/08
            jj = 4
            CALL det3d(tript,det)
            ddd(jj) = -det
            CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
            
        ! Plane 5
            
            DO  k = 1, 3
                tript(1,k) = tript(2,k)
            END DO
            tript(2,1) = w
            tript(2,2) = 0.0
            tript(2,3) = ho                                                    !WSU 4/4/08
            tript(3,1) = w                                                      !WSU 4/4/08
            tript(3,2) = l
            tript(3,3) = ho                                            !WSU 4/4/08
            jj = 5
            CALL det3d(tript,det)
            ddd(jj) = -det
            CALL findabc(tript,aaa(jj),bbb(jj),ccc(jj))
        !vf      DO  j = 1, jj
        !vf      WRITE(788,*) aaa(j), bbb(j), ccc(j), ddd(j)
        !vf      END DO 
        !  Write results
        
            !  iiii = iproc
            iiii = nn-1!iproc * ncbar_partition + i_partition - 1
            n_partition = iiii
            WRITE(intform(3:3),'(I1)') icount_digits(iiii)
            filename(1:10) = 'cbcch.out.'
            WRITE(filename(11:),intform) n_partition !iproc
            IF (ibin == 1) THEN
                OPEN(8,FILE=filename,FORM='unformatted',STATUS='unknown')
                WRITE(8)aaa
                WRITE(8)bbb
                WRITE(8)ccc
                WRITE(8)ddd
            ELSE
                OPEN(8,FILE=filename,FORM='formatted',STATUS='unknown')
                WRITE(8,*)aaa
                WRITE(8,*)bbb
                WRITE(8,*)ccc
                WRITE(8,*)ddd
            ENDIF
            close(8)
        ENDIF 
        END DO 
        
        RETURN
        END SUBROUTINE cbcch

!***********************************************************************
!  Date: 2007-07-30  Time: 08:40:35
 
!  Program linegen generates a set of points for each coset lines,
!  determines its angle relative to the bottom plane of the unit bar,
!  calculates the line parameters (slope, intercept),
!  determines the # of troughs to be generated along each coset,
!  establishes the position (x,y,z) of  a trough set and determines the length,
!  width and height of a trough set.
!  This code also generates and assigns textures to each trough set based 
!  on input proportion.

!  Input files:
!    1) linedata.dat: contains mean and varaince to generate values for parameters 
!       and input proportion
!    2) ubwl.out:contains length and width info. of unit bar. Output from trgena.f90
!    3) nub.out: contains the # of unit bars in a compound bar. Output from arcgen.out

!  Output files
!   1) lenstat.out-contains line parameters. (Input for merge.f90).
!   2) length.out: contains min and max of y,z of a line. Also contains length 
!      of a line (file used for checking purpose only).
!   3) angle.out: contains angle of each coset (file used for checking purpose only)
!   4) indicator.out: contains info. about texture assgined to each trough. 
!      (input for merge.f90)
!   5) location.out: contains location of a trough with respect to other troughs. 
!      (input for merge.f90)
!   6) coset.out: contains the # of cosets in each unit bar. 
!      (input for merge.f90)
!   7) line.out: contains the # of lines in each coset.(input for merge.f90)

!  A local coordinate system is used for this purpose. Some of the subroutines are from Scheibe (1993).

!   Last modified on June 4 2007
!*******************************************************************************************************************************
      SUBROUTINE linegen

      USE prcmpi
      USE bar

      INCLUDE 'mpif.h'

      REAL*4, ALLOCATABLE                               :: cthi(:,:,:)
      REAL*4, ALLOCATABLE                               :: thi(:,:,:)

      REAL*4, ALLOCATABLE                               :: mcapl(:)
      REAL*4, ALLOCATABLE                               :: ppdn(:,:,:)
      REAL*4, ALLOCATABLE                               :: slope(:,:,:)
      REAL*4, ALLOCATABLE                               :: yint(:,:,:)
      REAL*4, ALLOCATABLE                               :: llen(:,:)
      REAL*4, ALLOCATABLE                               :: cymin(:,:,:)
      REAL*4, ALLOCATABLE                               :: cymax(:,:,:)
      REAL*4, ALLOCATABLE                               :: czmin(:,:,:)
      REAL*4, ALLOCATABLE                               :: czmax(:,:,:)
      REAL*4, ALLOCATABLE                               :: pcymin(:,:,:)
      REAL*4, ALLOCATABLE                               :: pslope(:,:,:)
      REAL*4, ALLOCATABLE                               :: pyint(:,:,:)

      REAL*4, ALLOCATABLE                               :: tymin(:,:)
      REAL*4, ALLOCATABLE                               :: tymax(:,:)
      REAL*4, ALLOCATABLE                               :: tzmin(:,:)
      REAL*4, ALLOCATABLE                               :: tzmax(:,:)
      REAL*4, ALLOCATABLE                               :: tyint(:,:)
      REAL*4, ALLOCATABLE                               :: tslope(:,:)
      REAL*4, ALLOCATABLE                               :: tpymin(:,:)
      REAL*4, ALLOCATABLE                               :: tpymax(:,:)
      REAL*4, ALLOCATABLE                               :: tpyint(:,:)

      REAL*4, ALLOCATABLE                               :: avymax(:,:)
      REAL*4, ALLOCATABLE                               :: avymin(:,:) 
      REAL*4, ALLOCATABLE                               :: avzmax(:,:)
      REAL*4, ALLOCATABLE                               :: avzmin(:,:)
      REAL*4, ALLOCATABLE                               :: avpymin(:,:)

      REAL*4, ALLOCATABLE                               :: twidth(:,:,:)
      REAL*4, ALLOCATABLE                               :: txw(:)
      REAL*4, ALLOCATABLE                               :: txmin(:)
      REAL*4, ALLOCATABLE                               :: txmax(:)

      REAL*4, ALLOCATABLE                               :: slopeN(:,:)
      REAL*4, ALLOCATABLE                               :: yintN(:,:)
      REAL*4, ALLOCATABLE                               :: xminN(:,:)
      REAL*4, ALLOCATABLE                               :: xmaxN(:,:)
      REAL*4, ALLOCATABLE                               :: cyminN(:,:)
      REAL*4, ALLOCATABLE                               :: czmaxN(:,:)
      REAL*4, ALLOCATABLE                               :: twidthN(:,:)
      REAL*4, ALLOCATABLE                               :: ppdnN(:,:)
      REAL*4, ALLOCATABLE                               :: pslopeN(:,:)
      REAL*4, ALLOCATABLE                               :: pyintN(:,:)
      REAL*4, ALLOCATABLE                               :: delhN(:,:)
      REAL*4, ALLOCATABLE                               :: rllenN(:,:)

      REAL*4, ALLOCATABLE 								              :: xinfl_1(:,:)
      REAL*4, ALLOCATABLE                               :: yinfl_1(:,:)
      REAL*4, ALLOCATABLE								                :: a_1(:,:)
      REAL*4, ALLOCATABLE								                :: b_1(:,:)
      REAL*4, ALLOCATABLE								                :: width_1(:)
      REAL*4, ALLOCATABLE								                :: length_1(:)
      REAL*4, ALLOCATABLE								                :: hu_1(:)
      REAL*4, ALLOCATABLE								                :: zi_1(:)
      REAL*4, ALLOCATABLE                               :: zhmax_1(:)
      REAL*4, ALLOCATABLE                               :: zhmin_1(:) 
      INTEGER*2, ALLOCATABLE							              :: iub_1(:)

      REAL*4                                            :: xc
      REAL*4                                            :: dy
      REAL*4                                            :: delc
      REAL*4                                            :: cy
      REAL*4                                            :: cy1

      REAL*4                                            :: rpfs
      REAL*4                                            :: rpofg
      REAL*4                                            :: rpsg
      REAL*4                                            :: pfs
      REAL*4                                            :: pfsg
      REAL*4                                            :: pofg
      REAL*4                                            :: psg1
      REAL*4                                            :: psg
      REAL*4                                            :: rpbm
      REAL*4                                            :: rpm
      REAL*4                                            :: pbm
      REAL*4                                            :: pm


      REAL*4                                            :: dh
      REAL*4                                            :: dhm
      REAL*4                                            :: dhmax
      REAL*4                                            :: dhvar
      REAL*4                                            :: rfz
      REAL*4                                            :: rfzm
      REAL*4                                            :: rfzvar
      REAL*4                                            :: sety
      REAL*4                                            :: setym
      REAL*4                                            :: setyvar
      REAL*4                                            :: delx
      REAL*4                                            :: delxm
      REAL*4                                            :: delxvar

      REAL*4                                            :: rcymax
      REAL*4                                            :: nxlub
      REAL*4                                            :: twm
      REAL*4                                            :: twvar
      REAL*4                                            :: txc
      REAL*4                                            :: ttw
      REAL*4                                            :: tthi
      REAL*4                                            :: thim
      REAL*4                                            :: thivar
      REAL*4                                            :: rdiff
      REAL*4                                            :: rdiffm
      REAL*4                                            :: rdiffvar
	  !REAL*4                                            :: delhm
	  !REAL*4                                            :: twidthm
	  !REAL*4                                            :: llenm

      INTEGER*4, ALLOCATABLE                            :: nol(:,:)
      INTEGER*4, ALLOCATABLE                            :: ncs(:,:)
      INTEGER*4, ALLOCATABLE                            :: indt(:,:,:)
      INTEGER*4, ALLOCATABLE                            :: indtN(:,:)

      INTEGER*4                                         :: i
      INTEGER*4                                         :: ii
      INTEGER*4                                         :: j
      INTEGER*4                                         :: jj
      INTEGER*4                                         :: ij
      INTEGER*4                                         :: k
      INTEGER*4                                         :: ll

      INTEGER*4                                         :: n
      INTEGER*4                                         :: ncb
      INTEGER*4                                         :: mcb1

      INTEGER*4                                         :: ind
      INTEGER*4                                         :: nofg
      INTEGER*4                                         :: nbm
      INTEGER*4                                         :: nm
      INTEGER*4                                         :: nsg
      INTEGER*4                                         :: nsg1
      INTEGER*4                                         :: nfs
      INTEGER*4                                         :: notc
      INTEGER*4                                         :: zindt
      INTEGER*4                                         :: count1a
      INTEGER*4                                         :: tnfs
      INTEGER*4                                         :: tnofg
      INTEGER*4                                         :: tnsg
      INTEGER*4                                         :: tnbm
      INTEGER*4                                         :: tnm

      INTEGER*4                                         :: nlt
      INTEGER*4                                         :: ncob
      INTEGER*4                                         :: nltm
      INTEGER*4                                         :: ncobm
      INTEGER*4                                         :: ans

      INTEGER*4                                         ::nltn
      INTEGER*4                                         ::nltn1
      INTEGER*4                                         ::nltn2
      INTEGER*4                                         ::ncobn

      INTEGER*4                                         :: jtest


      INTEGER*4                                         :: ierr
      INTEGER*4                                         :: ncoset
      INTEGER*4                                         :: nlines
      INTEGER*4                                         :: nocl
      INTEGER*4                                         :: ndum
      INTEGER*4                                         :: nlow
      INTEGER*4                                         :: nhi
      INTEGER*4                                         ::icount
! 06-23-08
      INTEGER*4                                         :: itrough
      INTEGER*4                                         :: igl
      REAL*4                                            :: mcthi

      CHARACTER (LEN=12)                                :: filein
      CHARACTER (LEN=15)                                :: adum
      CHARACTER (LEN=64)                                :: filename
      CHARACTER (LEN=64)                                :: filename1
      CHARACTER (LEN=64)                                :: filename2

!      REAL*4                                            :: zindt
      REAL*4                                            :: zindtm
      REAL*4                                            :: zindtvar

      INTEGER*4                                         :: iiii
      INTEGER*4                                         :: na 

      !CHARACTER (LEN=15)                                :: adum

  1   FORMAT(1X,a)
 10   FORMAT(a12)
 15   FORMAT(1X,6(f8.6,2X))
 16   FORMAT(1X,7(f8.6,2X))
 20   FORMAT(1X,3(f11.8,3X))
 25   FORMAT(1X,4(f15.7,3X))
 30   FORMAT(1X,a,i1,a,f15.7)
 31   FORMAT(1X,i8)
 32   FORMAT(4(i4,2X),2X,3(f15.5,2X))
 34   FORMAT(2(i4,2X))
 35   FORMAT(1X,6(f12.8,2X))
 40   FORMAT(1X,i4,1X,i8,2X,i8,2X,i8,2X,i8)
 50   FORMAT(1X,a5,1X,f9.4,1X,i4,1X,i6)
 56   FORMAT(1X,i4,1X,i4,2X,i4,2X,i4,2X,i8)
 93   FORMAT(1X,f10.5,2X,i4)
 94   FORMAT(1X,i4,1X,i4,2X,i4,2X,i4)
 96   FORMAT(1X,f10.4,2X,f10.4)
 97   FORMAT(1X,f10.4)
 98   FORMAT(1X,f10.4)
 99   FORMAT(3(1X,i4),2X,f10.6)
100   FORMAT(1X,i4,1X,i4,2X,i8,2X,i4)
102   FORMAT(1X,i4,1X,i4,1X,i4)
101   FORMAT(1X,/)

      IF(iproc == 0) THEN
        WRITE(*,101)
        WRITE(*,101)
        WRITE(*,1) 'Set Characteristic Geometry: linegen'
        WRITE(*,1) '****************************'
        WRITE(*,101)
        
        OPEN(34,FILE='TSLOC.dat',STATUS='old') !Naum
        READ(34,*) h
        READ(34,*) dy,tdelc
        READ(34,*) ncoset
        READ(34,*) nlines
        READ(34,*) thim,thivar
        READ(34,*) rdiffm,rdiffvar
        READ(34,*) twm,twvar
        READ(34,*) delxm,delxvar
        ! READ(34,*) indtm, indtvar
        READ(34,*)  rfzm,rfzvar
        READ(34,*) pofg,psg,pfs,pbm,pm
        READ(34,*) nlow,nhi
        READ(34,*) dhm,dhvar,dhmax
        READ(34,*) setym,setyvar
        READ(34,*) mcthi
        nhio=nhi
        CLOSE(34)
      ENDIF      
 
      CALL MPI_BCAST(ncoset,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      adum = 'ncoset-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(nlines,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      adum = 'nlines-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(nlow,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      adum = 'nlow-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(nhi,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      adum = 'nhi-mpi'
      CALL CHECK(adum,ierr)


! Next broadcast reals

      CALL MPI_BCAST(h,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'h-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(dy,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'dy-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(tdelc,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'tdelc-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(thim,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'thim-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(thivar,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'thivar-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(rdiffm,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'rdiffm-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(rdiffvar,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'rdiffvar-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(twm,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'twm-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(twvar,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'twvar-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(delxm,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'delxm-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(delxvar,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'delxvar-mpi'
      CALL CHECK(adum,ierr)

      ! CALL MPI_BCAST(indtvar,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      ! adum = 'indtvar-mpi'
      ! CALL CHECK(adum,ierr)

      CALL MPI_BCAST(rfzm,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'rfzm-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(rfzvar,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'rfzvar-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(pofg,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'pofg-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(psg,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'psg-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(pfs,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'pfs-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(pbm,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'pbm-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(pm,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'pm-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(dhm,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'dhm-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(dhvar,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'dhvar-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(dhmax,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'dhmax-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(setym,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'setym-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(setyvar,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'setyvar-mpi'
      CALL CHECK(adum,ierr)

      CALL MPI_BCAST(mcthi,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      adum = 'mcthi-mpi'
      CALL CHECK(adum,ierr)

    !   igl=0
!      
      DO i_partition = 1,ncbar_process(iproc+1)
        ! iiii = iproc*ncbar_partition + i_partition-1
        IF (iproc<nc_determine) THEN
            iiii=iproc*(ncbar_partition+1)+i_partition-1
        ELSE
            ! kk1=myY*ncbar_partition+i_nb-1
            iiii= nc_determine*(ncbar_partition+1)+&
                  (iproc-nc_determine)*ncbar_partition+i_partition-1
        ENDIF
        IF (iiii < nocb) THEN
            WRITE(*,*) 'Calculate the linegen for component bar ', iiii
            nn = iiii + 1!iproc+1
            n_partition = iiii
    
            !Read the file including original zi, hu data
            WRITE(*,*) 'Allocate the arrays for component bar ', n_partition
            na = 2
            ALLOCATE (xinfl_1(noub_array(n_partition+1),na),STAT=IERR)
            ALLOCATE (yinfl_1(noub_array(n_partition+1),na),STAT=IERR)
            ALLOCATE (a_1(noub_array(n_partition+1),na),STAT=IERR)
            ALLOCATE (b_1(noub_array(n_partition+1),na),STAT=IERR)
            ALLOCATE (width_1(noub_array(n_partition+1)),STAT=IERR)
            ALLOCATE (length_1(noub_array(n_partition+1)),STAT=IERR)
            ALLOCATE (hu_1(noub_array(n_partition+1)),STAT=IERR)
            ALLOCATE (zi_1(noub_array(n_partition+1)),STAT=IERR)
            ALLOCATE (iub_1(noub_array(n_partition+1)),STAT=IERR)
    
    
            WRITE(*,*) 'Read data of arrays from previous sub in process ', n_partition
            WRITE(intform(3:3),'(I1)') icount_digits(iiii)
            filename1(1:11) = 'length.out.'
            WRITE(filename1(12:),intform) n_partition
            OPEN(15,FILE=filename1,FORM='unformatted',STATUS='unknown')
            READ(15)xinfl_1(1:noub_array(n_partition+1),1:na)
            READ(15)yinfl_1(1:noub_array(n_partition+1),1:na)
            READ(15)a_1(1:noub_array(n_partition+1),1:na)
            READ(15)b_1(1:noub_array(n_partition+1),1:na)
            READ(15)length_1(1:noub_array(n_partition+1))
            READ(15)width_1(1:noub_array(n_partition+1))
            READ(15)hu_1(1:noub_array(n_partition+1))
            READ(15)zi_1(1:noub_array(n_partition+1))
            READ(15)iub_1(1:noub_array(n_partition+1))
            CLOSE(15)
    
            ALLOCATE (zhmax_1(noub_array(n_partition+1)),STAT=IERR)
            ALLOCATE (zhmin_1(noub_array(n_partition+1)),STAT=IERR)
            ! WRITE(intform(3:3),'(I1)') icount_digits(iiii)
            filename2(1:11) = 'zhpara.out.'
            WRITE(filename2(12:),intform) n_partition
            OPEN(151,FILE=filename2,FORM='unformatted',STATUS='unknown')
            READ(151)zhmax_1(1:noub_array(n_partition+1))
            READ(151)zhmin_1(1:noub_array(n_partition+1))
            CLOSE(151)
            ! WRITE(*,*) 'The zhmin is ', zhmin_1(1),zhmin(1),iproc

            igl=0
        IF(igl == 1) THEN
            ! iiii = iproc
            ! WRITE(intform(3:3),'(I1)') icount_digits(iiii)
            filename(1:6) = 't.out.'
            WRITE(filename(7:),intform) n_partition
        IF(ibin == 1) THEN
            OPEN(8421,FILE=filename,FORM='unformatted',STATUS='unknown')
            filename(1:6) = 'w.out.'
            WRITE(filename(7:),intform) n_partition
            OPEN(8422,FILE=filename,FORM='unformatted',STATUS='unknown')
            filename(1:6) = 'l.out.'
            WRITE(filename(7:),intform) n_partition
            OPEN(8423,FILE=filename,FORM='unformatted',STATUS='unknown')
            filename(1:8) = 'lwt.out.'
            WRITE(filename(9:),intform) n_partition
            OPEN(8424,FILE=filename,FORM='unformatted',STATUS='unknown')
        ELSE
            OPEN(8421,FILE=filename,FORM='formatted',STATUS='unknown')
            filename(1:6) = 'w.out.'
            WRITE(filename(7:),intform) n_partition
            OPEN(8422,FILE=filename,FORM='formatted',STATUS='unknown')
            filename(1:6) = 'l.out.'
            WRITE(filename(7:),intform) n_partition
            OPEN(8423,FILE=filename,FORM='formatted',STATUS='unknown')
            filename(1:8) = 'lwt.out.'
            WRITE(filename(9:),intform) n_partition
            OPEN(8424,FILE=filename,FORM='formatted',STATUS='unknown')
            ENDIF
            ndh=0
            ntw=0
            nllen=0
        ENDIF

! Write out max values of length, width, thickness of a trough. Mar 2009 WSU
        ! WRITE(intform(3:3),'(I1)') icount_digits(iiii)
	      filename(1:12) = 'maxvals.out.'
        WRITE(filename(13:),intform) n_partition
        OPEN(42,FILE=filename,FORM='formatted',STATUS='unknown') 



!  Input files

!      OPEN(55,FILE='maxv.out',ACCESS='APPEND',STATUS='old') !nobody read it

!  These files are used as check on linegen; eliminate in large run

      IF (ichk == 1) THEN
        ! WRITE(intform(3:3),'(I1)') icount_digits(iiii)
        filename(1:11) = 'length.out.'
!        kproc=TRIM(iproc)
        WRITE(filename(12:),intform) n_partition
!        write(8732,*)filename
        IF(ibin == 1) THEN
          OPEN(15,FILE=filename,FORM='unformatted',STATUS='unknown')
        ELSE
          OPEN(15,FILE=filename,FORM='formatted',STATUS='unknown')
        ENDIF
        
        ! WRITE(intform(3:3),'(I1)') icount_digits(iiii)
        filename(1:10) = 'angle.out.'
        WRITE(filename(11:),intform) n_partition
        IF(ibin == 1) THEN
          OPEN(12,FILE=filename,FORM='unformatted',STATUS='unknown')
        ELSE
          OPEN(12,FILE=filename,FORM='formatted',STATUS='unknown')
        ENDIF
      ENDIF

!  Output files
      ! WRITE(intform(3:3),'(I1)') icount_digits(iiii)
      filename(1:12) = 'lenstat.out.'
      WRITE(filename(13:),intform) n_partition
      IF(ibin == 1) THEN
        OPEN(10,FILE=filename,FORM='unformatted',STATUS='unknown')
      ELSE
        OPEN(10,FILE=filename,FORM='formatted',STATUS='unknown')
      ENDIF

      ! WRITE(intform(3:3),'(I1)') icount_digits(iiii)
      filename(1:10) = 'coset.out.'
      WRITE(filename(11:),intform) n_partition
      IF(ibin == 1) THEN
        OPEN(80,FILE=filename,FORM='unformatted',STATUS='unknown')
      ELSE
        OPEN(80,FILE=filename,FORM='formatted',STATUS='unknown')
      ENDIF

      ! WRITE(intform(3:3),'(I1)') icount_digits(iiii)
      filename(1:9) = 'line.out.'
      WRITE(filename(10:),intform) n_partition
      IF(ibin == 1) THEN
        OPEN(82,FILE=filename,FORM='unformatted',STATUS='unknown')
      ELSE
        OPEN(82,FILE=filename,FORM='formatted',STATUS='unknown')
      ENDIF
      
      ! WRITE(intform(3:3),'(I1)') icount_digits(iiii)
      filename(1:11) = 'rdleng.out.'
      WRITE(filename(12:),intform) n_partition
      IF(ibin == 1) THEN
        OPEN(43,FILE=filename,FORM='unformatted',STATUS='unknown')
      ELSE
        OPEN(43,FILE=filename,FORM='formatted',STATUS='unknown')
      ENDIF

      ! WRITE(intform(3:3),'(I1)') icount_digits(iiii)
      filename(1:14) = 'indicator.out.'
      WRITE(filename(15:),intform) n_partition
      IF(ibin == 1) THEN
        OPEN(11, FILE=filename,FORM='unformatted',STATUS='unknown')
      ELSE
        OPEN(11, FILE=filename,FORM='formatted',STATUS='unknown')
      ENDIF
      
!vlf
      ! WRITE(intform(3:3),'(I1)') icount_digits(iiii)
      filename(1:11) = 'params.out.'
      WRITE(filename(12:),intform) n_partition
      OPEN(33, FILE=filename,FORM='unformatted',STATUS='unknown') !Naum
!vlf

!vlf
!      OPEN(33, FILE='params.out',FORM='unformatted',ACCESS='direct', &
!           RECL=16,STATUS='unknown')
!vlf
      
!  Read in height of trough, the mean and variance for generating some parameters,
!  Input proportion for fine sand (pfs), open framework gravel (pofg), sandy gravel (psg).

    !   IF (iproc == 0) THEN
    !     READ(34,*) h
    !     READ(34,*) dy,tdelc
    !     READ(34,*) ncoset
    !     READ(34,*) nlines
    !     READ(34,*) thim,thivar
    !     READ(34,*) rdiffm,rdiffvar
    !     READ(34,*) twm,twvar
    !     READ(34,*) delxm,delxvar
    !     READ(34,*) indtm, indtvar
    !     READ(34,*)  rfzm,rfzvar
    !     READ(34,*) pofg,psg,pfs,pbm,pm
    !     READ(34,*) nlow,nhi
    !     READ(34,*) dhm,dhvar,dhmax
    !     READ(34,*) setym,setyvar
    !     READ(34,*) mcthi
    !     nhio=nhi
        
    !   ENDIF
      !CLOSE(3)

! First broadcast integers

!  Assume one compound bar per processor; processor number+1 matches compound bar number
!  Will not work if more than one compound bar per processor

    !   nn = iproc+1
      
      i = 1

!  Read in the information about width and length of unit bars.

      ALLOCATE (mcapl(noub_array(n_partition+1)), STAT=IERR)
      adum = 'mcapl'
      CALL CHECK(adum,ierr)
       
      ALLOCATE (cthi(1,noub_array(n_partition+1),ncoset), STAT=IERR)
      adum = 'cthi'
      CALL CHECK(adum,ierr)
       
      ALLOCATE (thi(1,noub_array(n_partition+1),ncoset), STAT=IERR)
      adum = 'thi'
      CALL CHECK(adum,ierr)
       
      ALLOCATE (tymax(noub_array(n_partition+1),ncoset), STAT=IERR)
      adum = 'thymax'
      CALL CHECK(adum,ierr)
       
      ALLOCATE (tymin(noub_array(n_partition+1),ncoset), STAT=IERR)
      adum = 'tymin'
      CALL CHECK(adum,ierr)
       
      ALLOCATE (tzmin(noub_array(n_partition+1),ncoset), STAT=IERR)
      adum = 'tzmin'
      CALL CHECK(adum,ierr)
       
      ALLOCATE (tzmax(noub_array(n_partition+1),ncoset), STAT=IERR)
      adum = 'tzmax'
      CALL CHECK(adum,ierr)
       
      ALLOCATE (tpyint(noub_array(n_partition+1),ncoset), STAT=IERR)
      adum = 'tpyint'
      CALL CHECK(adum,ierr)
       
      ALLOCATE (tpymin(noub_array(n_partition+1),ncoset), STAT=IERR)
      adum = 'tpymin'
      CALL CHECK(adum,ierr)
       
      ALLOCATE (tpymax(noub_array(n_partition+1),ncoset), STAT=IERR)
      adum = 'tpymax'
      CALL CHECK(adum,ierr)
       
      ALLOCATE (tyint(noub_array(n_partition+1),ncoset), STAT=IERR)
      adum = 'tyint'
      CALL CHECK(adum,ierr)
       
      ALLOCATE (tslope(noub_array(n_partition+1),ncoset), STAT=IERR) 
      adum = 'tslope'
      CALL CHECK(adum,ierr)
       
      ALLOCATE (ncs(noub_array(n_partition+1),ncoset), STAT=IERR)
      adum = 'ncs'
      CALL CHECK(adum,ierr)

      ALLOCATE(nncoset(noub_array(n_partition+1)), STAT=IERR)
      adum = 'nncoset'
      CALL CHECK(adum,ierr)
       
      ALLOCATE (txw(nlines), STAT=IERR) 
      adum = 'txw'
      CALL CHECK(adum,ierr)
       
      ALLOCATE (txmin(nlines), STAT=IERR)
      adum = 'txmin'
      CALL CHECK(adum,ierr)
       
      ALLOCATE (txmax(nlines), STAT=IERR)
      adum = 'txmax'
      CALL CHECK(adum,ierr)

!  Initial angle assigned 22 deg or 0.384 radians (equal to angle of repose)
!  Generate decrements for subsequent coset lines

!  This is set to one for testing in other subroutines; needs to be zero here!
!        filename(1:13) = 'countind.out.'
!        WRITE(filename(14:),intform) iproc
!        OPEN(717,FILE=filename,FORM='formatted',STATUS='unknown') 
!	    filename(1:10) = 'count.out.'
!        WRITE(filename(11:),intform) iproc
!        OPEN(718,FILE=filename,FORM='formatted',STATUS='unknown') 

      jtest = 0
      itest=0
      ii = 1
      IF (itest == 1) THEN
        is1 = 23546
        is2 = 19465
        is3 = 2154
      ENDIF
      
      DO  k =1,noub_array(n_partition+1)
        mcapl(k)=1.00*length_1(k)
!     write(654,*)length(ii,k),ii,k
        DO  i=1,ncoset
  301     CALL norgen(tthi,thim,thivar,is1,is2,is3,jtest)
         ! IF(tthi <= 0.01) THEN
         !   tthi=thim
          IF(tthi< -0.00005.OR.tthi > 0.0001) THEN
            GO TO 301
          END IF
!         tthi=thim
          thi(ii,k,i)=tthi
!vf       if(ii==1.and.k==1)then
!     write(651,*)thi(ii,k,i),ii,k,i
!vf       endif
        END DO
      END DO
      IF (iproc == 0)WRITE(*,*)'Finished generating angle decrements'
      !write(7891,*)ncoset,iproc
! Determine each coset line angle

      nocs=0
      nq1 = 0
      ncb = 0
      
      loopCB:  DO ij=1,ncbar        
        mcb1=0
!        mxlub=300.0
!        wxub=70.0
        loopUB: DO i = 1,noub_array(n_partition+1)
          delc=tdelc
          mcb1 = mcb1 + 1
          ncb=0
          cy = 0.0
          loopCS:   DO j = 1, ncoset
            ncb=ncb+1
            k=ncb
            ncs(ij,i)=k
            IF(ncs(ij,i) > nocs) then
              nocs=ncs(ij,i)
            END IF
            cthi(ij,i,j)=delc - thi(ij,i,k)
!vf         write(572,*)ij,i,j
!vf         write(572,*)delc, thi(ij,i,k), cthi(ij,i,j)
            IF(cthi(ij,i,k) <= mcthi) THEN
              cthi(ij,i,k)=cthi(ij,i,k-1)
            END IF
 270        nxlub=1.0*mxlub  
!vf            if(ij==1.and.i==1) then
!vf            write(652,*)mxlub, nxlub,ij,i,j  
!vf            endif  
            IF(cy >= nxlub.OR.cy < dy) THEN
              ncb=ncb-1
              k=ncb
              ics=k
              ncs(ij,i)=k
              cycle loopUB
            END IF
            tymax(i,k) = cy + (h/TAN(cthi(ij,i,k)))
            tymin(i,k) = cy
            tzmax(i,k) = zhmax_1(i)
            tzmin(i,k) = zhmin_1(i)
            delc=cthi(ij,i,k)-0.00872
!            wxub=70.0

!  Determine slope and intercept of each line
!
            tslope(i,k) = (tymax(i,k) - tymin(i,k)) / &
                          (tzmin(i,k) - tzmax(i,k))
            tyint(i,k)= tymax(i,k) - tslope(i,k)*tzmin(i,k)
!
!    calculate the // line parameters
!

            delta = (dhm/COS(1.571428-cthi(ij,i,k)))
            tpymin(i,k)=tymin(i,k)+ delta
            tpyint(i,k)= tpymin(i,k) - tslope(i,k)*tzmax(i,k)
            tpymax(i,k)=tpyint(i,k)

!  Calculate next coset line using the ymin of the parallel line.
!  dcy is then used to offset the coset line.
            tdelta=delta
            diff1=tpymin(i,k)-tymin(i,k)
            dcy=rdiffm*diff1
            cy=tpymin(i,k)-dcy
            cy1=tymin(i,k)+0.5*diff1
            IF (cy < tymin(i,k).OR.cy > cy1) THEN
              cy=tymin(i,k)+diff1/2
            END IF
          END DO  loopCS
        END DO  loopUB
      END DO  loopCB 
     
      ttw=0.0
261   DO k =1,nlines
        txw(k)=twm
        ttw = txw(k)+ttw
        IF(k == 1) THEN
          txc = txw(1)/2
        END IF
        nocl = k
        txmin(k)=txc-(txw(k)/2.0)
        txmax(k)= txc+(txw(k)/2.0)
        IF (txmin(k) >= wxub) THEN
          nocl=k-1
          GOTO 262
        END IF
        txc=txmax(k)
      END DO   
       
262   DEALLOCATE(tymin, STAT=IERR)
      adum = 'tyminD'
      CALL CHECK(adum,ierr)

      DEALLOCATE(tymax, STAT=IERR)
      adum = 'tymaxD'
      CALL CHECK(adum,ierr)

      DEALLOCATE(tzmax, STAT=IERR)
      adum = 'tzmaxD'
      CALL CHECK(adum,ierr)

      DEALLOCATE(tyint, STAT=IERR)
      adum = 'tyintD'
      CALL CHECK(adum,ierr)

      DEALLOCATE(tslope, STAT=IERR)
      adum = 'tslopeD'
      CALL CHECK(adum,ierr)

      DEALLOCATE(tpymin, STAT=IERR)
      adum = 'tpyminD'
      CALL CHECK(adum,ierr)

      DEALLOCATE(tpymax, STAT=IERR)
      adum = 'tpymaxD'
      CALL CHECK(adum,ierr)

      DEALLOCATE(tpyint, STAT=IERR)
      adum = 'tpyintD'
      CALL CHECK(adum,ierr)

      DEALLOCATE(txw, STAT=IERR)
      adum = 'txwD'
      CALL CHECK(adum,ierr)

      DEALLOCATE(txmin, STAT=IERR)
      adum = 'txminD'
      CALL CHECK(adum,ierr)

      DEALLOCATE(txmax, STAT=IERR)
      adum = 'txmaxD'
      CALL CHECK(adum,ierr)

      DEALLOCATE(tzmin, STAT=IERR)
      adum = 'tzmin-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE(thi, STAT=IERR)
      adum = 'thi-deal'
      CALL CHECK(adum,ierr)


! May also want nocs to be an array if more than one compound bar per processor

      i = 1
      ALLOCATE (avymax(noub_array(n_partition+1),nocs), STAT=IERR)
      adum = 'avymax'
      CALL CHECK(adum,ierr)

      ALLOCATE (avymin(noub_array(n_partition+1),nocs), STAT=IERR)
      adum = 'avymin'
      CALL CHECK(adum,ierr)

      ALLOCATE (avzmin(noub_array(n_partition+1),nocs), STAT=IERR)
      adum = 'avzmin'
      CALL CHECK(adum,ierr)

      ALLOCATE (avzmax(noub_array(n_partition+1),nocs), STAT=IERR)
      adum = 'avzmax'
      CALL CHECK(adum,ierr)

      ALLOCATE (avpymin(noub_array(n_partition+1),nocs), STAT=IERR)
      adum = 'avpymin'
      CALL CHECK(adum,ierr)

      ALLOCATE (pcymin(noub_array(n_partition+1),nocs,nocl), STAT=IERR)
      adum = 'pcymin'
      CALL CHECK(adum,ierr)

      ALLOCATE (pslope(noub_array(n_partition+1),nocs,nocl), STAT=IERR)
      adum = 'pslope'
      CALL CHECK(adum,ierr)

      ALLOCATE (pyint(noub_array(n_partition+1),nocs,nocl), STAT=IERR)
      adum = 'pyint'
      CALL CHECK(adum,ierr)

      ALLOCATE (delh(noub_array(n_partition+1),nocs,nocl), STAT=IERR)
      adum = 'delh'
      CALL CHECK(adum,ierr)

      ALLOCATE (rllen(noub_array(n_partition+1),nocs,nocl), STAT=IERR)
      adum = 'rllen'
      CALL CHECK(adum,ierr)

      ALLOCATE (ppdn(noub_array(n_partition+1),nocs,nocl), STAT=IERR)
      adum = 'ppdn'
      CALL CHECK(adum,ierr)

      ALLOCATE (slope(noub_array(n_partition+1),nocs,nocl), STAT=IERR)
      adum = 'slope'
      CALL CHECK(adum,ierr)

      ALLOCATE (yint(noub_array(n_partition+1),nocs,nocl), STAT=IERR)
      adum = 'yint'
      CALL CHECK(adum,ierr)

      ALLOCATE (llen(noub_array(n_partition+1),nocs), STAT=IERR)
      adum = 'llen'
      CALL CHECK(adum,ierr)

      ALLOCATE (cymin(noub_array(n_partition+1),nocs,nocl), STAT=IERR)
      adum = 'cymin'
      CALL CHECK(adum,ierr)

      ALLOCATE (cymax(noub_array(n_partition+1),nocs,nocl), STAT=IERR)
      adum = 'cymax'
      CALL CHECK(adum,ierr)

      ALLOCATE (czmin(noub_array(n_partition+1),nocs,nocl), STAT=IERR)
      adum = 'czmin'
      CALL CHECK(adum,ierr)

      ALLOCATE (czmax(noub_array(n_partition+1),nocs,nocl), STAT=IERR)
      adum = 'czmax'
      CALL CHECK(adum,ierr)

      ALLOCATE (xmin(noub_array(n_partition+1),nocs,nocl), STAT=IERR)
      adum = 'xmin'
      CALL CHECK(adum,ierr)

      ALLOCATE (xmax(noub_array(n_partition+1),nocs,nocl), STAT=IERR)
      adum = 'xmax'
      CALL CHECK(adum,ierr)

      ALLOCATE (nol(noub_array(n_partition+1),nocs), STAT=IERR)
      adum = 'nol'
      CALL CHECK(adum,ierr)
      nol=0

      ALLOCATE (twidth(noub_array(n_partition+1),nocs,nocl), STAT=IERR)
      adum = 'twidth'
      CALL CHECK(adum,ierr)

      ALLOCATE(nline(noub_array(n_partition+1),nocs), STAT=IERR)
      adum = 'nline'
      CALL CHECK(adum,ierr)
      nline=0

      ALLOCATE (indt(noub_array(n_partition+1),nocs,nocl), STAT=IERR)
      adum = 'indt'
      CALL CHECK(adum,ierr)
      
      nl=0
      ij = 1
      jcm=0
      ! jcmarray(n_partition+1)=jcm
      DO ij=1,ncbar 
        DO i = 1,noub_array(n_partition+1)
        !  WRITE(intform(3:3),'(I1)') icount_digits(iiii)
        filename(1:11) = 'mlines.out.'
        WRITE(filename(12:),intform) n_partition
        OPEN(423,FILE=filename,FORM='formatted',STATUS='unknown')
           mcb1 = mcb1 + 1
           ncb=0
           cy = 0.0
           DO j = 1, ncs(ij,i)
             ncb=ncb+1
             k=ncb
 257         IF(cy >= mcapl(i).OR.cy < dy) THEN
               ncb=ncb-1
               k=ncb
               ics=k
               GO TO 260
             END IF
 280         CONTINUE
 
!  Determine min and max vlaues of each coset line
             avymax(i,k) = cy + (h/TAN(cthi(ij,i,k)))
             avymin(i,k) = cy
             avzmax(i,k) = zhmax_1(i)
             avzmin(i,k) = zhmin_1(i)
!            if(ij==1.and.i==1) then
!vf             write(561,*)ij,i,k
!vf             write(561,*)cy,h,cthi(ij,i,k)
!vf      write(561,*)avymax(i,k)
!vf       write(561,*)avymin(i,k)
!vf       write(561,*)avzmax(i,k)
!vf        write(561,*)avzmin(i,k)
!     endif
   
!   Generate average thickness(dh) of troughs in a coset boundary
      
      IF (itest == 1) THEN
        is1 = 23546
        is2 = 19465
        is3 = 2154
      ENDIF

 311         CALL norgen(dh,dhm,dhvar,is1,is2,is3,jtest)
           !     dh=dhm
             IF(dh < 0.1.or.dh>0.5) THEN
                GO TO 311
             END IF
           dh=dhm
!    calculate the // line parameters
      
             delta = (dh/COS(1.571428-cthi(ij,i,k)))
             avpymin(i,k)=avymin(i,k)+ delta
      
!  Calculate next coset line using the ymin of the parallel line.
!  dcy is then used to offset the coset line.
      
      IF (itest == 1) THEN
        is1 = 23546
        is2 = 19465
        is3 = 2154
      ENDIF
 302         CALL norgen(rdiff,rdiffm,rdiffvar,is1,is2,is3,jtest)
           !      rdiff=rdiffm
                IF(rdiff<=0.0.OR.rdiff > 0.2) THEN
               GO TO 302
             END IF
!              rdiff=rdiffm
             diff1=avpymin(i,k)-avymin(i,k)
             dcy=rdiff*diff1
             cy=avpymin(i,k)-dcy
!vf             write(561,*)cy, avpymin(i,k), dcy
             cy1=avymin(i,k)+0.5*diff1
             IF (cy < avymin(i,k)) THEN
               cy=avymin(i,k)+diff1/2
             END IF
!vf             write(561,*)cy,cy1, diff1
      
!  IF ymin exceeds the max length of unit bar along y, then stop generating coset lines for this unit bar.
             if(ij==1.and.i==1) then            
             endif
             IF(cy >= mcapl(i).OR.cy < dy) THEN
               ncb=ncb-1
               k=ncb
               GO TO 260
             END IF
           END DO
    
    
!     Generate width of troughs and establish min, max along x.
!     delx is used to offset the next trough along x
!
260        IF (ncb > nocs) THEN
             nocs=ncb
           ENDIF
!    
           ttwidth=0.0

           loop501:  DO  ii=1,ncb
             ttwidth=0.0
             DO jj=1,nocl

      IF (itest == 1) THEN
        is1 = 23546
        is2 = 19465
        is3 = 2154
      ENDIF

       419     CALL norgen(tw,twm,twvar,is1,is2,is3,jtest)
            !         tw=twm
               IF(tw<1.2.or.tw>5.6) THEN
                goto 419
               ENDIF 
               tw=twm
               IF (igl==1) THEN   
                 IF (ibin == 1) THEN 
                   write(8422)tw
                 ELSE
                   write(8422,*)tw
                 ENDIF
                 ntw=ntw+1
               ENDIF
               twidth(i,ii,jj) = tw
               ttwidth=twidth(i,ii,jj)+ttwidth
               IF(jj == 1) THEN
                 xc = twidth(i,ii,1)/2.0
               END IF
               nol(i,ii)=jj
                IF (nol(i,ii) > nl) then
                    nl=nol(i,ii)
                END IF
               xmin(i,ii,jj) = xc - twidth(i,ii,jj)/2.0
               xmax(i,ii,jj) = xc + twidth(i,ii,jj)/2.0

      IF (itest == 1) THEN
        is1 = 23546
        is2 = 19465
        is3 = 2154
      ENDIF

 701           CALL norgen(delx,delxm,delxvar,is1,is2,is3,jtest)
              !    delx=delxm
               IF(delx <= 0.OR.delx > 0.3) THEN
                 GO TO 701
               END IF
!               delx=delxm
!  If xmax exceeds the width of unit bar. then stop generating troughs along this coset.
!  proceed to next coset.
        
               IF (xmin(i,ii,jj) >= width_1(i)) THEN

!  maximum limit reached. STOP for this coset.

                 nol(i,ii)=jj-1
                 IF (nol(i,ii) > nl) THEN
                    nl=nol(i,ii)
                 END IF
                 
                 CYCLE loop501
               
               END IF
               xc = xmax(i,ii,jj)+delx
        
! Generate thickness for each trough along a coset line
        
      IF (itest == 1) THEN
        is1 = 23546
        is2 = 19465
        is3 = 2154
      ENDIF
312            CALL norgen(dh,dhm,dhvar,is1,is2,is3,jtest)
             !   dh=dhm
               IF(dh < 0.1.OR.dh > 0.5) THEN
                 GO TO 312
               END IF
               dh=dhm
               delh(i,ii,jj)=dh
               IF(igl==1) THEN
                 IF (ibin == 1) THEN
                   write(8421)dh
                 ELSE
                   write(8421,*)dh
                 ENDIF
                 ndh=ndh+1
               ENDIF
               w2t=twidth(i,ii,jj)/dh
               nr=nr+1
        
!  Calculate line parameters for each trough in a coset line
!  First generate a number (sety) to set off the trough from the cosest line
        
      IF (itest == 1) THEN
        is1 = 23546
        is2 = 19465
        is3 = 2154
      ENDIF

 313           CALL norgen(sety,setym,setyvar,is1,is2,is3,jtest)
              !      sety=setym
               IF(sety > 0.2) THEN
                 GO TO 313
               END IF
!               sety=setym
               cymin(i,ii,jj)=avymin(i,ii)-sety
               IF(ii /= 1) THEN
                 IF(cymin(i,ii,jj) <= avymin(i,ii-1)) THEN
                   cymin(i,ii,jj)=avymin(i,ii)
                 END IF
               END IF
               cymax(i,ii,jj) = cymin(i,ii,jj) + (h/TAN(cthi(ij,i,ii)))
               czmax(i,ii,jj) = zhmax_1(i)
               czmin(i,ii,jj) = zhmin_1(i)
               slope(i,ii,jj) = (cymax(i,ii,jj) - cymin(i,ii,jj))  &
                   /(czmin(i,ii,jj) - czmax(i,ii,jj))
               yint(i,ii,jj)= cymax(i,ii,jj) - slope(i,ii,jj) &
                              *czmin(i,ii,jj)
!    calculate the // line parameters
        
               pslope(i,ii,jj)=slope(i,ii,jj)
               delta = (delh(i,ii,jj)/COS(1.571428-cthi(ij,i,ii)))
               pcymin(i,ii,jj)=cymin(i,ii,jj)+ delta
               pyint(i,ii,jj)= pcymin(i,ii,jj) - slope(i,ii,jj) &
                               *czmax(i,ii,jj)
!
!  Calculate the perpendicular dist from a pt (miny, maxz) on the // line to the coset line
        
               anr= pcymin(i,ii,jj)-pslope(i,ii,jj)*czmax(i,ii,jj)- &
                    yint(i,ii,jj)
               bdr=SQRT(1+(pslope(i,ii,jj)**2))
               ppdn(i,ii,jj)=anr/bdr
            END DO
          END DO loop501

    
! Initialize for calculating proportions of textures for troughs
    
          ind=1
          mtemp=0
       mnd=1
       ind=1
       mcnum = mtemp + mnd
       mtemp=mcnum
       nntemp=0
       nnct=1
       ncobt=0
       tofg=0
       tsg=0
       tfs=0
       tnotc=0
       tnofg=0
       tnfs=0
       tnbm=0
       tnm=0
       tnsg=0
       tnotc=0
       notc=0
       ntfs=0
       nfs=0
       ntofg=0
       nofg=0
       ntsg=0
       nsg=0
       nltm=0
       ncobm=0
       
          DO  ii = 1,ncb
            ncob=ii
            nnct=nnct+nntemp
            n = nol(i,ii)
            nlt=1
          
            IF (itest == 1) THEN
              is1 = 23546
              is2 = 19465
              is3 = 2154
            ENDIF

            DO  jj = 1,nol(i,ii)             
             if(nol(i,ii)>nocl.or.nol(i,ii)<=0) THEN
              write(423,*)i,ii,jj,nol(i,ii)          
             endif
              notc=notc+1
              nlt=jj
              rllen(i,ii,jj) = SQRT((cymax(i,ii,jj)-cymin(i,ii,jj))**2 &
                  +(czmin(i,ii,jj)-czmax(i,ii,jj))**2)
              
!  Assign texture randomly

!              zindtm=3.0
!              zindtvar=5.0
 
!             CALL norgen (zindt, zindtm, zindtvar, is1, is2, is3,jtest) 
 251            CALL iungen(zindt,nlow,nhi,is1,is2,is3)             
              IF(zindt<1.OR.zindt>4)GOTO 251
              IF (zindt >= 1.AND.zindt<2)  THEN               
                tnfs=nfs+1
                !rpfs=(100*tnfs)/notc !Naum
                rpfs=(100*nfs)/notc !!Naum                
                IF(rpfs > pfs) THEN
                  GO TO 251
                END IF
                nfs=nfs+1
                indt(i,ii,jj) = 3
                GO TO 254
              ELSEIF (zindt >= 2.AND.zindt<3) THEN               
                tnofg=nofg+1
                !rpofg=(100*tnofg)/notc !Naum
                rpofg=(100*nofg)/notc !Naum                
                IF(rpofg > pofg) THEN
                  GO TO 251
                END IF
                nofg=nofg+1
                indt(i,ii,jj) = 4
                GO TO 252
              ELSEIF (zindt >= 3.AND.zindt<4) THEN               
                tnsg=nsg+1
                !rpsg=(100*tnsg)/notc
                rpsg=(100*nsg)/notc !Naum                
                IF(rpsg > psg) THEN
                  GO TO 251
                END IF
                indt(i,ii,jj) = 5
                nsg=nsg+1
                GO TO 254
!              ELSEIF (zindt >=4.0.AND.zindt<5.0) THEN
!                tnm=nm+1
!                rpm=(100*tnm)/notc
!                IF(rpm > pm) THEN
!                  GO TO 251
!                END IF
!                nm=nm+1
!                indt(i,ii,jj) = 21
!                GO TO 254
!              ELSEIF(zindt >=5.0.AND.zindt<6.0) THEN
!                tnsg=nsg+1
!                rpsg=(100*tnsg)/notc
!                IF(rpsg>psg) THEN
!                  GOTO 251
!                ENDIF
!                indt(i,ii,jj) = 5
!                nsg=nsg+1
!               GO TO 254
              ELSE
                IF(rpfs>=pfs.AND.rpofg>=pofg) THEN
                  indt(i,ii,jj) = 5
                  nsg=nsg+1
                 GOTO 254
                ENDIF
                GO TO 251
              END IF
              WRITE(*,*)'Generated indicators'

! Randomly generate reduction factor in thickness (rfz)
! for open framework gravel and fine sand troughs
        
      IF (itest == 1) THEN
        is1 = 23546
        is2 = 19465
        is3 = 2154
      ENDIF

 252          CALL norgen (rfz, rfzm, rfzvar, is1, is2, is3,jtest)
              !    rfz=rfzm
              IF (rfz <= 0.0.OR.rfz > 0.4) THEN
                GO TO 252
              END IF
!               rfz=rfzm
              tczmin=rfz*(abs(zhmin_1(i)-zhmax_1(i)))
              rcymax=slope(i,ii,jj)*tczmin+yint(i,ii,jj)
        
! Calculate reduced length for ofg and fine sand troughs
        
              rllen(i,ii,jj) =  SQRT((rcymax-cymin(i,ii,jj))**2  &
                  +(tczmin-czmax(i,ii,jj))**2)
        
          

 254          CONTINUE           
      
!  Calculate the total number of troughs
!  determine its location wrt to other troughs
        
              jcm = jcm+1
              
            END DO
            nline(i,ii) = nlt
!            nltn=maxval(nline)
!            nltn1=minval(nline)
!            IF(ntln>nocl.or.nltn1<0) then
!               write(423,*)nltn,nltn1,width(ij,i),ii
!            ENDIF
!           dth(ij,i)IF (ibin == 1) THEN
!             WRITE(82)nlt
!           ELSE
!             WRITE(82,'(I4)')nlt
!           ENDIF
            nntemp=nnct
            
            llen(i,ii)=SQRT((avymax(i,ii)-avymin(i,ii))**2  &
                +(avzmin(i,ii)-avzmax(i,ii))**2)
            IF(igl==1) THEN
              IF (ibin == 1) THEN
                write(8423)llen(i,ii)
              ELSE
                write(8423,*)llen(i,ii)
              ENDIF
              nllen=nllen+1
            ENDIF
! Comment this out for large runs

            IF (ichk == 1) THEN
              IF (ibin == 1) THEN
                WRITE(12)n_partition+1,i,ii,cthi(ij,i,ii)
                WRITE(15)n_partition+1,i,ii,n
                WRITE(15)avymin(i,ii),avymax(i,ii)
                WRITE(15)avzmin(i,ii),avzmax(i,ii)
                WRITE(15)llen(i,ii)
              ELSE
                WRITE(12,'(3(I4,2X),1PE13.6)')iproc+1,i,ii,cthi(ij,i,ii)
                WRITE(15,'(3(I4,2X))')n_partition+1,i,ii,n
                WRITE(15,*)avymin(i,ii),avymax(i,ii)
                WRITE(15,*)avzmin(i,ii),avzmax(i,ii)
                WRITE(15,*)llen(i,ii)
              ENDIF
            ENDIF
          END DO
          nncoset(i) = ncob       
        END DO
      nltm = MAXVAL(nline)
      ncobm = MAXVAL(nncoset)
      nltn2 = MAXVAL(nol)
      nltn1 = MINVAL(nol)
      ncobn = MINVAL(nncoset)
       filename(1:11) = 'nlines.out.'
        WRITE(filename(12:),intform) n_partition
        OPEN(422,FILE=filename,FORM='formatted',STATUS='unknown')
write(422,*)nltm,ncobm,nltn2,ncobn,nltn1
close(422)
close(423)
      END DO
      jcmarray(n_partition+1)=jcm
      ! jcm=0
! Write parameters for array sizes here; get maximum # first.

      nltm = MAXVAL(nol)
      ncobm = MAXVAL(nncoset)     
      delhm=maxval(delh)
      twidthm=maxval(twidth)
      llenm=maxval(llen)
      
      WRITE(42,*)delhm, twidthm, llenm
      
	   nn = n_partition+1 !iproc+1

!vlf
      WRITE(33)noub_array(n_partition+1),jcm,ncobm,nltm !Naum      
!vlf
!      WRITE(33,REC=nn)noub,jcm,ncobm,nltm
!vlf
      
      DEALLOCATE (avymax, STAT=IERR)
      adum = 'avymax-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (avymin, STAT=IERR)
      adum = 'avymin-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (avzmin, STAT=IERR)
      adum = 'avzmin-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (avzmax, STAT=IERR)
      adum = 'avzmax-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (avpymin, STAT=IERR)
      adum = 'avpymin-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (pcymin, STAT=IERR)
      adum = 'pcymin-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (llen, STAT=IERR)
      adum = 'llen-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (cymax, STAT=IERR)
      adum = 'cymax-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (czmin, STAT=IERR)
      adum = 'czmin-deal'
      CALL CHECK(adum,ierr)
      
!  Allocate new arrays for write

      ALLOCATE(slopeN(ncbar,jcm), STAT=IERR)
      adum = 'slopeN'
      CALL CHECK(adum,ierr)

      ALLOCATE(yintN(ncbar,jcm), STAT=IERR)
      adum = 'yintN'
      CALL CHECK(adum,ierr)

      ALLOCATE(xminN(ncbar,jcm), STAT=IERR)
      adum = 'xminN'
      CALL CHECK(adum,ierr)

      ALLOCATE(xmaxN(ncbar,jcm), STAT=IERR)
      adum = 'xmaxN'
      CALL CHECK(adum,ierr)

      ALLOCATE(cyminN(ncbar,jcm), STAT=IERR)
      adum = 'cyminN'
      CALL CHECK(adum,ierr)

      ALLOCATE(czmaxN(ncbar,jcm), STAT=IERR)
      adum = 'cymaxN'
      CALL CHECK(adum,ierr)

      ALLOCATE(twidthN(ncbar,jcm), STAT=IERR)
      adum = 'twidthN'
      CALL CHECK(adum,ierr)

      ALLOCATE(ppdnN(ncbar,jcm), STAT=IERR)
      adum = 'ppdnN'
      CALL CHECK(adum,ierr)

      ALLOCATE(pslopeN(ncbar,jcm), STAT=IERR)
      adum = 'pslopeN'
      CALL CHECK(adum,ierr)

      ALLOCATE(pyintN(ncbar,jcm), STAT=IERR)
      adum = 'pyintN'
      CALL CHECK(adum,ierr)

      ALLOCATE(indtN(ncbar,jcm), STAT=IERR)
      adum = 'indtN'
      CALL CHECK(adum,ierr)
      
      ALLOCATE(delhN(ncbar,jcm), STAT=IERR)
      adum = 'delhN'
      CALL CHECK(adum,ierr)
      
      ALLOCATE(rllenN(ncbar,jcm), STAT=IERR)
      adum = 'rllenN'
      CALL CHECK(adum,ierr)     
       WRITE(intform(3:3),'(I1)') icount_digits(iiii)
        filename(1:10) = 'count.out.'
        WRITE(filename(11:),intform) n_partition
        OPEN(718,FILE=filename,FORM='formatted',STATUS='unknown') 

      ij=1
      DO  i=1,noub_array(n_partition+1)
       tnbm=0
       tnm=0
       tnsg=0
       nbm=0
       nm=0
       nsg1=0
       rpsg=0
       rpbm=0
       rpm=0
       psg1=100-pm-pbm
       count1a=0
        DO ii=1,nncoset(i)
          DO jj=1,nol(i,ii)

           IF(indt(i,ii,jj)==5) THEN
            count1a=count1a+1
            271 continue
             CALL iungen(zindt,nlow,nhi,is1,is2,is3)
!              IF(zindt<1.OR.zindt>4)GOTO 271
              IF (zindt >=1.AND.zindt<2) THEN
                tnbm=nbm+1
                rpbm=(100*tnbm)/count1a
                IF(rpbm > pbm) THEN
                  GO TO 271
                END IF
                nbm=nbm+1
                indt(i,ii,jj) = 20
              ELSEIF (zindt >=2.AND.zindt<3) THEN
                tnm=nm+1
                rpm=(100*tnm)/count1a
!                IF(rpm > pm) THEN
!                  GO TO 271
!                END IF
                nm=nm+1
                indt(i,ii,jj) = 21
              ELSEIF(zindt >=3.AND.zindt<4) THEN
                tnsg=nsg1+1
                rpsg=(100*tnsg)/count1a
                IF(rpsg>psg1) THEN
                  GOTO 271
                ENDIF
                indt(i,ii,jj) = 5
                nsg1=nsg1+1
              ELSE 
                goto 271
              ENDIF
             endif
!            IF(iproc==0)write(671,*)i,ii,jj,indt(i,ii,jj)
!            IF(iproc==1)write(672,*)i,ii,jj,indt(i,ii,jj)
             EnDDO
           ENDDO
         write(718,*)i
         write(718,*)nsg1
         write(718,*)nbm
         write(718,*)nm
         write(718,*)count1a
         write(718,*)psg1
         write(718,*)pbm
         write(718,*)pm
!          IF(iproc==0) WRITE(677,*)nbm,nm,count1a,nsg1,i
!          IF(iproc==0) WRITE(678,*)nbm,nm,count1a,nsg1,i
          ENDDO
          
         close(718)
      nq1 = 0
      ij = 1
      icount=0
      DO  i=1,noub_array(n_partition+1)
        DO  j=1,nncoset(i)
!             nlineN(ij,i,j)=
          DO  k=1,nol(i,j)
             nq1 = nq1+1
             slopeN(ij,nq1) = slope(i,j,k)
             yintN(ij,nq1) = yint(i,j,k)
             xminN(ij,nq1) = xmin(i,j,k)
             xmaxN(ij,nq1) = xmax(i,j,k)
             cyminN(ij,nq1) = cymin(i,j,k)
             czmaxN(ij,nq1) = czmax(i,j,k)
             twidthN(ij,nq1) = twidth(i,j,k)
             ppdnN(ij,nq1) = ppdn(i,j,k)
             pslopeN(ij,nq1) = pslope(i,j,k)
             pyintN(ij,nq1) = pyint(i,j,k)
             indtN(ij,nq1) = indt(i,j,k)
             rllenN(ij,nq1) = rllen(i,j,k)
             delhN(ij,nq1) = delh(i,j,k)
!         if(ij==1.and.i==1) then
           
!vf          if(indt(i,j,k)==4) then
!vf          WRITE(540,*) slope(i,j,k),i,j,k
!vf          WRITE(540,*)yint(i,j,k),i,j,k
!vf          WRITE(540,*) xmin(i,j,k),i,j,k
!vf          WRITE(540,*)xmax(i,j,k),i,j,k
!vf          WRITE(540,*) cymin(i,j,k),i,j,k
          !WRITE(540,*) cymax(i,j,k),i,j,k
          !WRITE(540,*) czmin(i,j,k),i,j,k
!vf          WRITE(540,*) czmax(i,j,k),i,j,k
!vf          WRITE(540,*) twidth(i,j,k),i,j,k
!vf          WRITE(540,*) ppdn(i,j,k),i,j,k
!vf          WRITE(540,*) pslope(i,j,k),i,j,k
!vf          WRITE(540,*)pyint(i,j,k), i,j,k
!vf          WRITE(541,*) indt(i,j,k),i,j,k
!          if(indt(i,j,k)==4) then
!vf           write(539,*)indt(i,j,k)
!vf           write(539,*)i,j,k
!          icount=icount+1
!vf          endif
          !WRITE(542,*) nncoset(ij,i),ij,i
!vf          WRITE(543,*) nline(ij,i,j), ij,i,j
!         endif
           ENDDO
         ENDDO
!         write(717,*)i
!         write(717,*)nfs
!         write(717,*)nofg
!         write(717,*)nbm
!         write(717,*)nm
!         write(717,*)nsg
       ENDDO
       
      nltm = MAXVAL(nol)
      ncobm = MAXVAL(nncoset)       
!         close(717)
!vf       write(712,*)icount
!    Writeout results for length statistics

      IF(ibin == 1) THEN
        WRITE(80)nncoset(1:noub_array(n_partition+1))
!       WRITE(80)nncoset(1,j), j=1,noub)
!vf     write(800,*)nncoset(1,1:noub)
        WRITE(82)nol(1:noub_array(n_partition+1),1:ncobm)
        !WRITE(828,*)nol(1:noub,1:ncobm)
        WRITE(10)slopeN
        WRITE(10)yintN
        WRITE(10)xminN
        WRITE(10)xmaxN
        WRITE(10)cyminN
        WRITE(10)czmaxN
        WRITE(10)twidthN
        WRITE(10)ppdnN
        WRITE(10)pslopeN
        WRITE(10)pyintN
        WRITE(11)indtN
        WRITE(43)rllenN
        WRITE(43)delhN
!        WRITE(4304,*)rllenN
!        WRITE(4304,*)delhN
!vf        do jj=1,jcm
!vf         if(indtN(ij,jj)==4) then        
!vf          WRITE(538,*)indtN(ij,jj)
!vf         endif
!vf        enddo
       ! WRITE(*,*)'FINISHED WRITING LENGTH STAT'
!
!
!          WRITE(540,*) (slopeN(ij,jj),jj = 1,jcm)
!          WRITE(540,*)(yintN(ij,jj),jj = 1,jcm)
!          WRITE(540,*) (xminN(ij,jj),jj = 1,jcm)
!          WRITE(540,*)(xmaxN(ij,jj),jj = 1,jcm)
!          WRITE(540,*) (cyminN(ij,jj),jj = 1,jcm)
!          WRITE(540,*) (czmaxN(ij,jj),jj = 1,jcm)
!          WRITE(540,*) (twidthN(ij,jj),jj = 1,jcm)
!          WRITE(540,*) (ppdnN(ij,jj),jj = 1,jcm)
!          WRITE(540,*) (pslopeN(ij,jj),jj = 1,jcm)
!          WRITE(540,*)(pyintN(ij,jj), jj=1,jcm)
!          WRITE(541,*) (indtN(ij,jj), jj=1,jcm)
!          WRITE(542,*) (nncoset(1,i), i=1,noub)
!          WRITE(543,*) ((nline(1,i,j),i=1, noub), j=1,ncob)
        

      ELSE
        WRITE(80,*)(nncoset(i),i=1,noub_array(n_partition+1))
        WRITE(82,*)((nline(i,j),i=1,noub_array(n_partition+1)),j=1,ncob)
        WRITE(10,*)(slopeN(ij,jj),jj = 1,jcm)
        WRITE(10,*)(yintN(ij,jj),jj = 1,jcm)
        WRITE(10,*)(xminN(ij,jj),jj = 1,jcm)
        WRITE(10,*)(xmaxN(ij,jj),jj = 1,jcm)
        WRITE(10,*)(cyminN(ij,jj),jj = 1,jcm)
        WRITE(10,*)(czmaxN(ij,jj),jj = 1,jcm)
        WRITE(10,*)(twidthN(ij,jj),jj = 1,jcm)
        WRITE(10,*)(ppdnN(ij,jj),jj = 1,jcm)
        WRITE(10,*)(pslopeN(ij,jj),jj = 1,jcm)
        WRITE(10,*)(pyintN(ij,jj),jj = 1,jcm)
        WRITE(11,*)(indtN(ij,jj),jj = 1,jcm)
        WRITE(43,*)(rllenN(ij,jj),jj = 1,jcm)
        WRITE(43,*)(delhN(ij,jj),jj = 1,jcm)
      ENDIF

      DEALLOCATE (slope, STAT=IERR)
      adum = 'slope-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (yint, STAT=IERR)
      adum = 'yint-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (cymin, STAT=IERR)
      adum = 'cymin-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (czmax, STAT=IERR)
      adum = 'czmax-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (twidth, STAT=IERR)
      adum = 'twidth-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (ppdn, STAT=IERR)
      adum = 'ppdn-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (pslope, STAT=IERR)
      adum = 'pslope-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (pyint, STAT=IERR)
      adum = 'pyint-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (nol, STAT=IERR)
      adum = 'nol-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (slopeN, STAT=IERR)
      adum = 'slopeN-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (yintN, STAT=IERR)
      adum = 'yintN-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (cyminN, STAT=IERR)
      adum = 'cyminN-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (czmaxN, STAT=IERR)
      adum = 'czmaxN-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (twidthN, STAT=IERR)
      adum = 'twidthN-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (ppdnN, STAT=IERR)
      adum = 'ppdnN-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (pslopeN, STAT=IERR)
      adum = 'pslopeN-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (pyintN, STAT=IERR)
      adum = 'pyintN-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE (indtN, STAT=IERR)
      adum = 'indtN-deal'
      CALL CHECK(adum,ierr)
      
      DEALLOCATE (rllenN, STAT=IERR)
      adum = 'rllenN-deal'
      CALL CHECK(adum,ierr)
      
      DEALLOCATE (delhN, STAT=IERR)
      adum = 'delhN-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE(mcapl,STAT=IERR)
      adum = 'mcapl-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE(cthi,STAT=IERR)
      adum = 'cthi-deal'
      CALL CHECK(adum,ierr)
      DEALLOCATE(nncoset,STAT=IERR)
      adum = 'nncoset-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE(delh,STAT=IERR)
      adum = 'delh-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE(rllen,STAT=IERR)
      adum = 'delh-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE(xmin,STAT=IERR)
      adum = 'delh-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE(xmax,STAT=IERR)
      adum = 'delh-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE(xminN,STAT=IERR)
      adum = 'delh-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE(xmaxN,STAT=IERR)
      adum = 'delh-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE(nline,STAT=IERR)
      adum = 'delh-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE(indt,STAT=IERR)
      adum = 'delh-deal'
      CALL CHECK(adum,ierr)

      DEALLOCATE(ncs,STAT=IERR)
      adum = 'delh-deal'
      CALL CHECK(adum,ierr)
      
     IF(igl==1) then
      IF (ibin == 1) THEN
        WRITE(8424) ndh
        WRITE(8424)ntw
        WRITE(8424)nllen
      ELSE
        WRITE(8424,*)'ndh', ndh
        WRITE(8424,*)'ntw', ntw
        WRITE(8424,*)'nllen', nllen
      ENDIF
      CLOSE(8421)
      CLOSE(8422)
      CLOSE(8423)
      CLOSE(8424)
    endif     


      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
      CLOSE(15)
      CLOSE(17)
      CLOSE(20)
      CLOSE(33)
      CLOSE(80)
      CLOSE(82)

      CLOSE(42)
      CLOSE(43)
!vf      CLOSE(9781)  
      DEALLOCATE(xinfl_1, STAT=IERR)
      DEALLOCATE(yinfl_1, STAT=IERR)
      DEALLOCATE(a_1, STAT=IERR)
      DEALLOCATE(b_1, STAT=IERR)
      DEALLOCATE(width_1, STAT=IERR)
      DEALLOCATE(length_1, STAT=IERR)
      DEALLOCATE(hu_1,STAT=IERR)
      DEALLOCATE(zi_1,STAT=IERR)
      DEALLOCATE(iub_1,STAT=IERR)
      DEALLOCATE(zhmax_1,STAT=IERR)
      DEALLOCATE(zhmin_1,STAT=IERR)
    ENDIF
    ENDDO

999   CONTINUE

      RETURN

      END SUBROUTINE linegen

!***********************************************************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2007-07-28  Time: 15:27:45

!     Program TROUGH This code generates the plane coefficients for each trough 
!     based on length, width and height generated by linegen.for.
!     A set of 31 planes are generated for each trough set.
!     A local coordinate system (different from that of linegen.for) is 
!     used for this purpose.

!     Input files:
!     1) data.dat : contains ratio's for determining y and z coordinates. Also 
!        contains the height of trough set (if non-random).

!     Output files
!     1) out1.out: contains plane coefficients for each plane defining a trough 
!        set. (output for merge.for)

!     Last modification on June 04 2007

!***********************************************************************
      SUBROUTINE tplane
  
        USE prcmpi
        USE bar
  
        INCLUDE 'mpif.h'
        
        REAL*4                                            :: r1
        REAL*4                                            :: r2
        REAL*4                                            :: r3
        REAL*4                                            :: r4
        REAL*4                                            :: r5
        REAL*4                                            :: det
        REAL*4                                            :: xr
        REAL*4                                            :: zin
        REAL*4                                            :: thb
  
        REAL*4                                            :: rx1
        REAL*4                                            :: rx2
        REAL*4                                            :: rx3
        REAL*4                                            :: rx4
        REAL*4                                            :: rx5
        REAL*4                                            :: rx6
        REAL*4                                            :: rx7
        REAL*4                                            :: rx8
        REAL*4                                            :: rx9
        REAL*4                                            :: rx10
        REAL*4                                            :: rx11
        REAL*4                                            :: rx12
        REAL*4                                            :: rx13
        REAL*4                                            :: rx14
        REAL*4                                            :: rx15
        REAL*4                                            :: rx16
        REAL*4                                            :: rx17
        REAL*4                                            :: rx18
        REAL*4                                            :: rx19
  
        REAL*4                                            :: z1
        REAL*4                                            :: z2
        REAL*4                                            :: z3
        REAL*4                                            :: z4
        REAL*4                                            :: z5
        REAL*4                                            :: z6
        REAL*4                                            :: z7
        REAL*4                                            :: z8
        REAL*4                                            :: z9
        REAL*4                                            :: z10
        REAL*4                                            :: z11
        REAL*4                                            :: z12
        REAL*4                                            :: z13
        REAL*4                                            :: z14
        REAL*4                                            :: z15
        REAL*4                                            :: z16
        REAL*4                                            :: z17
        REAL*4                                            :: z18
        REAL*4                                            :: z19
  
        REAL*4                                            :: pz1
        REAL*4                                            :: pz2
        REAL*4                                            :: pz3
        REAL*4                                            :: pz4
        REAL*4                                            :: pz5
        REAL*4                                            :: pz6
        REAL*4                                            :: pz7
        REAL*4                                            :: pz8
        REAL*4                                            :: pz9
        REAL*4                                            :: pz10
        REAL*4                                            :: pz11
        REAL*4                                            :: pz12
        REAL*4                                            :: pz13
        REAL*4                                            :: pz14
        REAL*4                                            :: pz15
        REAL*4                                            :: pz16
        REAL*4                                            :: pz17
        REAL*4                                            :: pz18
        REAL*4                                            :: pz19
  
        REAL*4                                            :: y1
        REAL*4                                            :: y2
        REAL*4                                            :: y3
        REAL*4                                            :: y4
        REAL*4                                            :: y5
        REAL*4                                            :: y6
        REAL*4                                            :: y7
        REAL*4                                            :: y8
        REAL*4                                            :: y9
        REAL*4                                            :: y10
        REAL*4                                            :: y11
        REAL*4                                            :: y12
        REAL*4                                            :: y13
        REAL*4                                            :: y14
        REAL*4                                            :: y15
        REAL*4                                            :: y16
        REAL*4                                            :: y17
        REAL*4                                            :: y18
        REAL*4                                            :: y19
  !  	  REAL*4                                            :: delhm
  !	  REAL*4                                            :: twidthm
  !	  REAL*4                                            :: llenm
      REAL*4                                            :: tmin
  
        REAL*4, ALLOCATABLE                               :: aaa(:,:,:)
        REAL*4, ALLOCATABLE                               :: bbb(:,:,:)
        REAL*4, ALLOCATABLE                               :: ccc(:,:,:)
        REAL*4, ALLOCATABLE                               :: ddd(:,:,:)
        REAL*4                                            :: tript(3,3)
  
        INTEGER*4                                         :: i_partition
      INTEGER*4											:: n_partition
      
        INTEGER*4                                         :: nn
  
        INTEGER*4                                         :: jj
        INTEGER*4                                         :: jk
        INTEGER*4                                         :: ilength
        INTEGER*4                                         :: ierr
        INTEGER*4                                         :: nt
        INTEGER*4                                         :: nq1,nq2, nq3
        INTEGER*4                                         :: itrough
  
        INTEGER*4                                         :: ncbar_1
  ! 06-23-08
         REAL*4                                           :: rb
         REAL*4                                           :: dr1
         REAL*4                                           :: dr2
  !	   REAL*4                                           :: thb
  ! 
        CHARACTER (LEN=12)                                :: filein
        CHARACTER (LEN=12)                                :: fileot
        CHARACTER (LEN=64)                                :: filename
        CHARACTER (LEN=15)                                :: adum
  
  
   1    FORMAT(1X,a)
   10   FORMAT(a12)
   15   FORMAT(1X,6(f8.6,2X))
   16   FORMAT(1X,7(f8.6,2X))
   20   FORMAT(1X,3(f11.8,3X))
   25   FORMAT(1X,4(f15.7,3X))
   30   FORMAT(1X,a,i1,a,f15.7)
   34   FORMAT(2(i4,2X))
   35   FORMAT(4(1X,i4))
   94   FORMAT(1X,i4,1X,i4,2X,i4,2X,i4)
   96   FORMAT(1X,f10.4,2X,f10.4)
   97   FORMAT(1X,f10.4)
   98   FORMAT(1X,f10.4)
   100  FORMAT(1X,i4,2X,i4,2X,i4)
   101  FORMAT(1X,/)
  
        IF (iproc == 0) THEN
          WRITE(*,101)
          WRITE(*,101)
          WRITE(*,1) 'Set Characteristic Geometry: tplane'
          WRITE(*,1) '****************************'
          WRITE(*,101)
          !WRITE(*,1) 'Enter input filename: ' !Naum
          !READ(*,10) filein      !Naum
          !OPEN(3,FILE=filein,STATUS='old')   !Naum
          OPEN(35,FILE='TSPLANE.dat',STATUS='old') !Naum
          !READ(*,10) fileot  !Naum
          READ(35,*) r1, r2, r3, r4,r5
          READ(35,*)dr1, dr2
          READ(35,*)rb
          CLOSE(35)
          fileot(1:8) = 'OUT1.out' !Naum
          ilength = LEN_TRIM(fileot)
        ENDIF
  
        !CALL MPI_Bcast(ilength,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        !adum = 'ilength-mpi'
        !CALL CHECK(adum,ierr)
  
        !CALL MPI_Bcast(fileot,ilength,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
        !adum = 'fileot-mpi'
        !CALL CHECK(adum,ierr)
        !iiii = iproc
  !     Read input data
  
        ! IF (iproc == 0) THEN
        !   READ(35,*) r1, r2, r3, r4,r5
        !   READ(35,*)dr1, dr2
        !   READ(35,*)rb
        ! ENDIF
     
        itrough=1 !Naum
        IF(itrough ==0 ) THEN
          GOTO 999
        ENDIF
        tmin=0.0
  ! READ in max values of length, width, thickness of a trough. Mar 2009 WSU
  
  
         
        CALL MPI_Bcast(r1,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'r1'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(r2,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'r2'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(r3,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'r3'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(r4,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'r4'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(r5,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'r5'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(rb,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'thb'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(dr1,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'dr1'
        CALL CHECK(adum,ierr)
       
        CALL MPI_Bcast(dr2,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'dr2'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(itrough,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        adum = 'itrough'
        CALL CHECK(adum,ierr)
  
      CALL MPI_Bcast(tmin,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'tmin'
        CALL CHECK(adum,ierr)
  
      ! CALL MPI_Bcast(twidthm,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      !   adum = 'twidthm'
      !   CALL CHECK(adum,ierr)
  
      ! CALL MPI_Bcast(delhm,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      !   adum = 'delhm'
      !   CALL CHECK(adum,ierr)
  
      ! CALL MPI_Bcast(llenm,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      !   adum = 'llenm'
      !   CALL CHECK(adum,ierr)
  
        ! DO i_partition = 0, ncbar_partition - 1
        !   iiii = iproc * ncbar_partition + i_partition  
        !   IF (iiii < nocb) THEN
        !     n_partition = iiii
        !     WRITE(intform(3:3),'(I1)') icount_digits(iiii)
        !     filename = 'OUT1.out.'
        !     WRITE(filename(10:),intform) n_partition  !iproc
        !     IF(ibin == 1) THEN
        !       OPEN(8,FILE=filename,FORM='unformatted', STATUS='unknown')
          
        !     ELSE
        !       OPEN(8,FILE=filename,FORM='formatted', STATUS='unknown')
        !     ENDIF
        !   ENDIF
        ! END DO
  ! ! Allocate arrays
  
  !       nt = 31
  !       ALLOCATE(aaa(ncbar,jcm,nt), STAT=IERR)
  !       adum = 'aaa'
  !       CALL CHECK(adum,ierr)
  
  !       ALLOCATE(bbb(ncbar,jcm,nt), STAT=IERR)
  !       adum = 'bbb'
  !       CALL CHECK(adum,ierr)
  
  !       ALLOCATE(ccc(ncbar,jcm,nt), STAT=IERR)
  !       adum = 'ccc'
  !       CALL CHECK(adum,ierr)
  
  !       ALLOCATE(ddd(ncbar,jcm,nt), STAT=IERR)
  !       adum = 'ddd'
  !       CALL CHECK(adum,ierr)
            
  ! The x- coordinates for all troughs remain the same.
        DO i_partition = 0, ncbar_process(iproc+1) - 1
        !   iiii = iproc * ncbar_partition + i_partition 
        IF (iproc<nc_determine) THEN
            iiii=iproc*(ncbar_partition+1)+i_partition
        ELSE
            ! kk1=myY*ncbar_partition+i_nb-1
            iiii= nc_determine*(ncbar_partition+1)+&
                    (iproc-nc_determine)*ncbar_partition+i_partition
        ENDIF 
          IF (iiii < nocb) THEN
            n_partition = iiii
            WRITE(intform(3:3),'(I1)') icount_digits(iiii)
            filename = 'OUT1.out.'
            WRITE(filename(10:),intform) n_partition  !iproc
            IF(ibin == 1) THEN
              OPEN(8,FILE=filename,FORM='unformatted', STATUS='unknown')
          
            ELSE
              OPEN(8,FILE=filename,FORM='formatted', STATUS='unknown')
            ENDIF
            WRITE(*,*) 'The calculation for tsplane is ',n_partition
  
            ! Allocate arrays
  
            nt = 31
            ! ncbar_1 = ncbar_partition
            ALLOCATE(aaa(ncbar,jcmarray(n_partition+1),nt), STAT=IERR)
            adum = 'aaa'
            CALL CHECK(adum,ierr)
      
            ALLOCATE(bbb(ncbar,jcmarray(n_partition+1),nt), STAT=IERR)
            adum = 'bbb'
            CALL CHECK(adum,ierr)
      
            ALLOCATE(ccc(ncbar,jcmarray(n_partition+1),nt), STAT=IERR)
            adum = 'ccc'
            CALL CHECK(adum,ierr)
      
            ALLOCATE(ddd(ncbar,jcmarray(n_partition+1),nt), STAT=IERR)
            adum = 'ddd'
            CALL CHECK(adum,ierr)
            nq1 = 0
            nq2=0
            DO  ij=1,ncbar
            DO  i=1,1     
              DO  j=1,1
              DO  k=1,1
                nq1 = nq1+1
        !  Starting to assign coordinates
            
                xr=ABS(tmin-twidthm)/5
                rx7=0.0
                rx6=rx7
                rx5=rx7
                rx4=rx7+ABS(tmin-twidthm)/2
                rx8=rx7+xr
                rx15=rx8
                rx16=rx8
                rx9=rx8+xr
                rx14=rx9
                rx17=rx9
                rx10=rx9+xr
                rx13=rx10
                rx18=rx10
                rx11=rx10+xr
                rx12=rx11
                rx19=rx11
                rx1=rx7+ABS(tmin-twidthm)
                rx2=rx1
                rx3=rx1
        !    !       if(ij==1.and.i==1) then
        !    !         write(591,*)rx1, rx2, rx3, rx4
        !    !         write(591,*)rx5, rx6, rx7, rx8
        !    !         write(591,*)rx9, rx10, rx11, rx12
        !             write(591,*)rx13, rx14, rx15, rx16
        !             write(591,*)rx17, rx18, rx19
        !          endif
        ! The z coordinates for # 1,2,3,4,5,6 and 7 equals zero and do not change any of the troughs.
            
                z1=0.0
                z2=z1
                z3=z1
                z4=z1
                z5=z1
                z6=z1
                z7=z1
                pz1=z1
                pz2=pz1
                pz3=pz1
                pz4=pz1
                pz5=pz1
                pz6=pz1
                pz7=pz1
            
        !To assign the remaining z values
            
              zin=0.0
              z10=zin+ABS(delhm)
              thb=rb*delhm
              pz10=z10-thb
              z9= z10
              pz9= pz10
              z11 = r5*z10
              pz11 = z11-thb
              z8 = z11
              z13=dr1*z10
              pz8 = pz11
              pz13=dr1*pz10
              z14=z13
              z12=dr1*z11
              pz14=pz13
              pz12=dr1*pz11
              z15=z12
              z18=dr2*z10
              pz15=pz12
              pz18=dr2*pz10
              z17=z18
              z19=dr2*z11
              pz17=pz18
              pz19=dr2*pz11
              z16=z19
              pz16=pz19
              
              
          !  To assign the y values
              
                  y1 = 0.0
                  y2 = y1 + r2*llenm
                  y3 = y2 + r3*llenm
                  y4 = y3 + r4*llenm
                  y5=y3
                  y6=y2
                  y7=y1
                  y8=y1
                  y9=y1
                  y10=y1
                  y11=y1
                  y12=y2
                  y13=y2
                  y14=y2
                  y15=y2
                  y16=y3
                  y17=y3
                  y18=y3
                  y19=y3
          !      if(i==1.and.j==1.and.k==1) then
          !         write(599,*)delh(i,j,k)
          !         write(599,*)rb
          !         write(599,*)z1,z2,z3,z4
          !         write(599,*)z5,z6,z7
          !         write(599,*)z8,z9,z10,z11
          !         write(599,*)z12,z13,z14,z15
          !         write(599,*)z16,z17,z18,z19
          !         write(599,*)pz1,pz2,pz3,pz4
          !         write(599,*)pz5,pz6,pz7
          !         write(599,*)pz8,pz9,pz10,pz11
          !         write(599,*)pz12,pz13,pz14,pz15
          !         write(599,*)pz16,pz17,pz18,pz19
          !         write(599,*)
          !         write(599,*)y1, y2, y3, y4
          !         write(599,*)y5, y6, y7, y8
          !         write(599,*)y9, y10, y11, y12
          !         write(599,*)y13, y14, y15, y16
          !         write(599,*)y17, y18, y19
          !       endif
  
          ! For plane 1
  
                  tript(1,1) = rx1
                  tript(1,2) = y1
                  tript(1,3) = z1
                  tript(2,1) = rx2
                  tript(2,2) = y2
                  tript(2,3) = z2
                  tript(3,1) = rx12
                  tript(3,2) = y12
                  tript(3,3) = z12
                  jj = 1
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
          !              if(ij==1.and.i==1) then
          !               write(592,*)aaa(ij,nq1,jj)
          !               write(592,*)bbb(ij,nq1,jj)
          !               write(592,*)ccc(ij,nq1,jj)
          !               write(592,*)ddd(ij,nq1,jj)
          !              endif            
          ! For plane 2
  
                  tript(1,1) = rx11
                  tript(1,2) = y11
                  tript(1,3) = z11
                  tript(2,1) = rx12
                  tript(2,2) = y12
                  tript(2,3) = z12
                  tript(3,1) = rx13
                  tript(3,2) = y13
                  tript(3,3) = z13
                  jj = 2
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
          !vf               if(ij==1.and.i==1) then
          !vf               write(592,*)aaa(ij,nq1,jj)
          !vf               write(592,*)bbb(ij,nq1,jj)
          !vf               write(592,*)ccc(ij,nq1,jj)
          !vf               write(592,*)ddd(ij,nq1,jj)
          !vf              endif
                
          ! For plane 3
  
                  tript(1,1) = rx10
                  tript(1,2) = y10
                  tript(1,3) = z10
                  tript(2,1) = rx13
                  tript(2,2) = y13
                  tript(2,3) = z13
                  tript(3,1) = rx14
                  tript(3,2) = y14
                  tript(3,3) = z14
                  jj = 3
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
                
          ! For plane 4
              
                  tript(1,1) = rx9
                  tript(1,2) = y9
                  tript(1,3) = z9
                  tript(2,1) = rx14
                  tript(2,2) = y14
                  tript(2,3) = z14
                  tript(3,1) = rx15
                  tript(3,2) = y15
                  tript(3,3) = z15
                  jj = 4
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
                
          ! For plane 5
              
                  tript(1,1) = rx8
                  tript(1,2) = y8
                  tript(1,3) = z8
                  tript(2,1) = rx15
                  tript(2,2) = y15
                  tript(2,3) = z15
                  tript(3,1) = rx6
                  tript(3,2) = y6
                  tript(3,3) = z6
                  jj = 5
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
                
          ! For plane 6
              
                  tript(1,1) = rx2
                  tript(1,2) = y2
                  tript(1,3) = z2
                  tript(2,1) = rx3
                  tript(2,2) = y3
                  tript(2,3) = z3
                  tript(3,1) = rx19
                  tript(3,2) = y19
                  tript(3,3) = z19
                  jj = 6
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 7
              
                  tript(1,1) = rx12
                  tript(1,2) = y12
                  tript(1,3) = z12
                  tript(2,1) = rx19
                  tript(2,2) = y19
                  tript(2,3) = z19
                  tript(3,1) = rx18
                  tript(3,2) = y18
                  tript(3,3) = z18
                  jj = 7
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 8
              
                  tript(1,1) = rx13
                  tript(1,2) = y13
                  tript(1,3) = z13
                  tript(2,1) = rx18
                  tript(2,2) = y18
                  tript(2,3) = z18
                  tript(3,1) = rx17
                  tript(3,2) = y17
                  tript(3,3) = z17
                  jj = 8
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 9
              
                  tript(1,1) = rx14
                  tript(1,2) = y14
                  tript(1,3) = z14
                  tript(2,1) = rx17
                  tript(2,2) = y17
                  tript(2,3) = z17
                  tript(3,1) = rx16
                  tript(3,2) = y16
                  tript(3,3) = z16
                  jj = 9
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 10
              
                  tript(1,1) = rx15
                  tript(1,2) = y15
                  tript(1,3) = z15
                  tript(2,1) = rx16
                  tript(2,2) = y16
                  tript(2,3) = z16
                  tript(3,1) = rx5
                  tript(3,2) = y5
                  tript(3,3) = z5
                  jj = 10
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 11
              
                  tript(1,1) = rx3
                  tript(1,2) = y3
                  tript(1,3) = z3
                  tript(2,1) = rx4
                  tript(2,2) = y4
                  tript(2,3) = z4
                  tript(3,1) = rx19
                  tript(3,2) = y19
                  tript(3,3) = z19
                  jj = 11
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
          !vf               if(ij==1.and.i==1) then
          !vf               write(581,*)aaa(ij,nq1,jj)
          !vf               write(581,*)bbb(ij,nq1,jj)
          !vf               write(581,*)ccc(ij,nq1,jj)
          !vf               write(581,*)ddd(ij,nq1,jj)
          !vf              endif
  
              
              
          ! For plane 12
              
                  tript(1,1) = rx19
                  tript(1,2) = y19
                  tript(1,3) = z19
                  tript(2,1) = rx4
                  tript(2,2) = y4
                  tript(2,3) = z4
                  tript(3,1) = rx18
                  tript(3,2) = y18
                  tript(3,3) = z18
                  jj = 12
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 13
              
                  tript(1,1) = rx18
                  tript(1,2) = y18
                  tript(1,3) = z18
                  tript(2,1) = rx4
                  tript(2,2) = y4
                  tript(2,3) = z4
                  tript(3,1) = rx17
                  tript(3,2) = y17
                  tript(3,3) = z17
                  jj = 13
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 14
              
                  tript(1,1) = rx17
                  tript(1,2) = y17
                  tript(1,3) = z17
                  tript(2,1) = rx4
                  tript(2,2) = y4
                  tript(2,3) = z4
                  tript(3,1) = rx16
                  tript(3,2) = y16
                  tript(3,3) = z16
                  jj = 14
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 15
              
                  tript(1,1) = rx16
                  tript(1,2) = y16
                  tript(1,3) = z16
                  tript(2,1) = rx4
                  tript(2,2) = y4
                  tript(2,3) = z4
                  tript(3,1) = rx5
                  tript(3,2) = y5
                  tript(3,3) = z5
                  jj = 15
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! Plane 16. Upper boundary plane.
              
                  tript(1,1) = rx1
                  tript(1,2) = y1
                  tript(1,3) = z1
                  tript(2,1) = rx4
                  tript(2,2) = y4
                  tript(2,3) = z4
                  tript(3,1) = rx7
                  tript(3,2) = y7
                  tript(3,3) = z7
                  jj = 16
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          !   Parallel planes for trough
          !
          ! For plane 17 (pp1)
              
                  tript(1,1) = rx1
                  tript(1,2) = y1
                  tript(1,3) = pz1
                  tript(2,1) = rx2
                  tript(2,2) = y2
                  tript(2,3) = pz2
                  tript(3,1) = rx12
                  tript(3,2) = y12
                  tript(3,3) = pz12
          !             
                  jj = 17
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 18
              
                  tript(1,1) = rx11
                  tript(1,2) = y11
                  tript(1,3) = pz11
                  tript(2,1) = rx12
                  tript(2,2) = y12
                  tript(2,3) = pz12
                  tript(3,1) = rx13
                  tript(3,2) = y13
                  tript(3,3) = pz13
                  jj = 18
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 19
              
                  tript(1,1) = rx10
                  tript(1,2) = y10
                  tript(1,3) = pz10
                  tript(2,1) = rx13
                  tript(2,2) = y13
                  tript(2,3) = pz13
                  tript(3,1) = rx14
                  tript(3,2) = y14
                  tript(3,3) = pz14
                  jj = 19
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
          !vf               if(ij==1.and.i==1) then
          !vf               write(582,*)aaa(ij,nq1,jj)
          !vf               write(582,*)bbb(ij,nq1,jj)
          !vf               write(582,*)ccc(ij,nq1,jj)
          !vf               write(582,*)ddd(ij,nq1,jj)
          !vf              endif
  
          ! For plane 20
              
                  tript(1,1) = rx9
                  tript(1,2) = y9
                  tript(1,3) = pz9
                  tript(2,1) = rx14
                  tript(2,2) = y14
                  tript(2,3) = pz14
                  tript(3,1) = rx15
                  tript(3,2) = y15
                  tript(3,3) = pz15
                  jj = 20
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 21
              
                  tript(1,1) = rx8
                  tript(1,2) = y8
                  tript(1,3) = pz8
                  tript(2,1) = rx15
                  tript(2,2) = y15
                  tript(2,3) = pz15
                  tript(3,1) = rx6
                  tript(3,2) = y6
                  tript(3,3) = pz6
                  jj = 21
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 22
              
                  tript(1,1) = rx2
                  tript(1,2) = y2
                  tript(1,3) = pz2
                  tript(2,1) = rx3
                  tript(2,2) = y3
                  tript(2,3) = pz3
                  tript(3,1) = rx19
                  tript(3,2) = y19
                  tript(3,3) = pz19
                  jj = 22
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 23
              
                  tript(1,1) = rx12
                  tript(1,2) = y12
                  tript(1,3) = pz12
                  tript(2,1) = rx19
                  tript(2,2) = y19
                  tript(2,3) = pz19
                  tript(3,1) = rx18
                  tript(3,2) = y18
                  tript(3,3) = pz18
                  jj = 23
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 24
              
                  tript(1,1) = rx13
                  tript(1,2) = y13
                  tript(1,3) = pz13
                  tript(2,1) = rx18
                  tript(2,2) = y18
                  tript(2,3) = pz18
                  tript(3,1) = rx17
                  tript(3,2) = y17
                  tript(3,3) = pz17
                  jj = 24
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 25
              
                  tript(1,1) = rx14
                  tript(1,2) = y14
                  tript(1,3) = pz14
                  tript(2,1) = rx17
                  tript(2,2) = y17
                  tript(2,3) = pz17
                  tript(3,1) = rx16
                  tript(3,2) = y16
                  tript(3,3) = pz16
                  jj = 25
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 26
              
                  tript(1,1) = rx15
                  tript(1,2) = y15
                  tript(1,3) = pz15
                  tript(2,1) = rx16
                  tript(2,2) = y16
                  tript(2,3) = pz16
                  tript(3,1) = rx5
                  tript(3,2) = y5
                  tript(3,3) = pz5
                  jj = 26
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 27
              
                  tript(1,1) = rx3
                  tript(1,2) = y3
                  tript(1,3) = pz3
                  tript(2,1) = rx4
                  tript(2,2) = y4
                  tript(2,3) = pz4
                  tript(3,1) = rx19
                  tript(3,2) = y19
                  tript(3,3) = pz19
                  jj = 27
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 28
              
                  tript(1,1) = rx19
                  tript(1,2) = y19
                  tript(1,3) = pz19
                  tript(2,1) = rx4
                  tript(2,2) = y4
                  tript(2,3) = pz4
                  tript(3,1) = rx18
                  tript(3,2) = y18
                  tript(3,3) = pz18
                  jj = 28
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 29
              
                  tript(1,1) = rx18
                  tript(1,2) = y18
                  tript(1,3) = pz18
                  tript(2,1) = rx4
                  tript(2,2) = y4
                  tript(2,3) = pz4
                  tript(3,1) = rx17
                  tript(3,2) = y17
                  tript(3,3) = pz17
                  jj = 29
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 30
              
                  tript(1,1) = rx17
                  tript(1,2) = y17
                  tript(1,3) = pz17
                  tript(2,1) = rx4
                  tript(2,2) = y4
                  tript(2,3) = pz4
                  tript(3,1) = rx16
                  tript(3,2) = y16
                  tript(3,3) = pz16
                  jj = 30
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
              
          ! For plane 31
              
                  tript(1,1) = rx16
                  tript(1,2) = y16
                  tript(1,3) = pz16
                  tript(2,1) = rx4
                  tript(2,2) = y4
                  tript(2,3) = pz4
                  tript(3,1) = rx5
                  tript(3,2) = y5
                  tript(3,3) = pz5
                  jj = 31
                  CALL det3d(tript,det)
                  ddd(ij,nq1,jj) = -det
                  CALL findabc(tript,aaa(ij,nq1,jj),bbb(ij,nq1,jj), &
                        ccc(ij,nq1,jj))
        !vf               if(ij==1.and.i==1) then
        !vf               write(583,*)aaa(ij,nq1,jj)
        !vf               write(583,*)bbb(ij,nq1,jj)
        !vf               write(583,*)ccc(ij,nq1,jj)
        !vf               write(583,*)ddd(ij,nq1,jj)
        !vf              endif
  
      !vf      if(i==1.and.j==1) then 
      !vf     if(k==1.or.k==4.or.k==10) then
      !vf        do ik=1,31
      !vf           write(576,*)aaa(ij,nq1,ik), bbb(ij,nq1,ik), ccc(ij,nq1,ik),&
      !vf                       ddd(ij,nq1,ik)
      !vf        enddo 
      !vf       endif
      !vf      endif 
      !       endif
              END DO
            END DO
          END DO
         !   nq2=MAXVAL(nq1)
  
          IF (ibin == 1) THEN
            WRITE(8)(aaa(ij,1,ii),ii=1,16)
            WRITE(8)(bbb(ij,1,ii),ii=1,16)
            WRITE(8)(ccc(ij,1,ii),ii=1,16)
            WRITE(8)(ddd(ij,1,ii),ii=1,16)
            WRITE(8)(aaa(ij,1,ii),ii=17,31)
            WRITE(8)(bbb(ij,1,ii),ii=17,31)
            WRITE(8)(ccc(ij,1,ii),ii=17,31)
            WRITE(8)(ddd(ij,1,ii),ii=17,31)
          ELSE
            WRITE(8,*)(aaa(ij,1,ii),ii=1,16)
            WRITE(8,*)(bbb(ij,1,ii),ii=1,16)
            WRITE(8,*)(ccc(ij,1,ii),ii=1,16)
            WRITE(8,*)(ddd(ij,1,ii),ii=1,16)
            WRITE(8,*)(aaa(ij,1,ii),ii=17,31)
            WRITE(8,*)(bbb(ij,1,ii),ii=17,31)
            WRITE(8,*)(ccc(ij,1,ii),ii=17,31)
            WRITE(8,*)(ddd(ij,1,ii),ii=17,31)
          ENDIF
          WRITE(*,*)'DONE!!!'
          DEALLOCATE(aaa, STAT=IERR)
          adum = 'aaa-deal'
          CALL CHECK(adum,ierr)
    
          DEALLOCATE(bbb, STAT=IERR)
          adum = 'bbb-deal'
          CALL CHECK(adum,ierr)
    
          DEALLOCATE(ccc, STAT=IERR)
          adum = 'ccc-deal'
          CALL CHECK(adum,ierr)
    
          DEALLOCATE(ddd, STAT=IERR)
          adum = 'ddd-deal'
          CALL CHECK(adum,ierr)
          END DO
  
  
  
          ENDIF
        END DO
        999   CONTINUE
        RETURN
  
        END SUBROUTINE tplane


!***********************************************************************

!     Subroutine NORGEN uses a set of subroutines to generate a single
!     random variable from a normal distribution with mean and
!     variance as specified in the call statement.  The
!     subroutine TSEED should be called before the first
!     call of NORGEN in a main program, in order to set the
!     initial seed values from the clock.

!***********************************************************************

      SUBROUTINE norgen(x,mean,var,is1,is2,is3,jtest)

      REAL*4, INTENT(OUT)                               :: x
      REAL*4, INTENT(IN)                                :: mean
      REAL*4, INTENT(IN OUT)                            :: var
      REAL*4                                            :: xarr(2)
                                                        
      INTEGER*4, INTENT(IN OUT)                         :: is1
      INTEGER*4, INTENT(IN OUT)                         :: is2
      INTEGER*4, INTENT(IN OUT)                         :: is3
      INTEGER*4                                         :: n
      INTEGER*4                                         :: jtest

! Number of uniform random variables to be generated = 2

      n = 2

! First, AS183 is called to generate two U[0,1] numbers
! (and to modify the seed values)

      CALL as183(is1,is2,is3,n,xarr)

! Next, BOXMUL is called to perform the Box-Muller transformation,
! which transforms the uniform variables into standard normal
! variables.

      CALL boxmul(n,xarr,xarr)

! Then, transform the standard normal into the distribution with
! specified mean and variance.

      x = xarr(1)*SQRT(var) + mean

      IF (jtest == 1) THEN
        is1 = 23546
        is2 = 19465
        is3 = 2154
      ENDIF

      RETURN
      END SUBROUTINE norgen
!***********************************************************************

      SUBROUTINE tseed(is1,is2,is3)

! TSEED constructs three integer seed values for the random number
! generator AS183 using the internal clock of the computer (PC/AT).

! The arguments must be type INTEGER.

! TSEED requires subroutine GETTIM, which is supplied by PROFORT.LIB.

! Written 8/85 by Tom Black


      INTEGER*4, INTENT(OUT)                     :: is1
      INTEGER*4, INTENT(OUT)                     :: is2
      INTEGER*4, INTENT(OUT)                     :: is3

      REAL*4                                     :: r(3)

      REAL*8                                     :: dseed
      REAL*8                                     :: d2p31m
      REAL*8                                     :: d2p31

      DATA  d2p31m/2147483647.d0/
      DATA  d2p31/2147483648.d0/

      CALL gettim(ihr,imin,isec,ihsec)
      dseed=ihsec*982743+isec*1666+ihr*1

      DO  i=1,3
        dseed=DMOD(16807.d0*dseed,d2p31m)
        r(i)=dseed/d2p31
      END DO

      is1=INT(30000*r(1))
      is2=INT(30000*r(2))
      is3=INT(30000*r(3))

      RETURN
      END SUBROUTINE tseed

!***********************************************************************

      SUBROUTINE boxmul(n,unif,norml)

! BOXMUL performs the standard Box-Muller normal transfomation on the
! rectangular [0,1] random numbers in UNIF and places them in NORML.
! N is the number of generated numbers in array UNIF.  UNIF and NORML
! have assumed dimensions corresponding to those in the calling program.

! THE TEMP STRUCTURE PERMITS UNIF AND NORML TO BE THE SAME VECTOR.

! N is type INTEGER.

! Arrays corresponding to UNIF and NORML must be single-precision REAL
! in the calling program.

! BOXMUL requires no FUNCTIONs or SUBROUTINEs.

! Written 3/85 by Tom Black.  REVISED 4/86.


      INTEGER*4, INTENT(IN)                      :: n
      REAL*4, INTENT(IN)                         :: unif(*)
      REAL*4, INTENT(OUT)                        :: norml(*)


100   DO  i=1,n-1,2
        temp1=unif(i)
        temp2=unif(i+1)
        norml(i)=SQRT(0.0-2.0*LOG(temp1))*COS(6.28318*temp2)
        norml(i+1)=SQRT(0.0-2.0*LOG(temp1))*SIN(6.28318*temp2)
      END DO

      RETURN

      END SUBROUTINE boxmul

!***********************************************************************

      SUBROUTINE as183(ix,iy,iz,n,unif)

!  AS183 uses three multiplicative congruential random number generators in
!  'parallel' to generate random numbers rectagularly distributed on [0,1].
!  The method is taken from 'An Efficient and Portable Pseudo-random Number
!  Generator' by B.A. Wichmann and I.D. Hill, Applied Statistics, algorithm
!  number AS 183, 1982.

!  IX, IY and IZ are integer seeds with 1.LE.seed.LE.30000.  N is the number
!  of random numbers to be generated.  The N generated numbers are placed in
!  array UNIF, which has assumed dimension equal to the actual array arg-
!  ument in the calling program.

!  IX, IY, IZ, and N are type integer.  UNIF is single precision.

!  AS183 requires no FUNCTIONs or SUBROUTINEs.

!  Written 3/85 by Tom Black


      INTEGER*4, INTENT(OUT)                     :: ix
      INTEGER*4, INTENT(OUT)                     :: iy
      INTEGER*4, INTENT(OUT)                     :: iz
      INTEGER*4, INTENT(IN)                      :: n

      REAL*4, INTENT(OUT)                        :: unif(*)
      real*4                                     :: num1, num2


100   DO  i=1,n
        ix=171*MOD(ix,177)-2*(ix/177)
        iy=172*MOD(iy,176)-35*(iy/176)
        iz=170*MOD(iz,178)-63*(iz/178)
        IF(ix < 0) ix=ix+30269
        IF(iy < 0) iy=iy+30307
        IF(iz < 0) iz=iz+30323
        unif(i)=AMOD(FLOAT(ix)/30269.0+FLOAT(iy)/30307.0+ FLOAT(iz)/30323.0,1.0)
         if(unif(i)<=0.0) then
          num1=unif(i)
          num2=0.99
          unif(i)=max(num1, num2)
         endif

      END DO

      RETURN

      END SUBROUTINE as183

!***********************************************************************

      SUBROUTINE iungen(n,nlow,nhi,ix,iy,iz)


      INTEGER*4, INTENT(OUT)                     :: n
      INTEGER*4, INTENT(IN)                      :: nlow
      INTEGER*4, INTENT(IN)                      :: nhi
      INTEGER*4, INTENT(IN OUT)                  :: ix
      INTEGER*4, INTENT(IN OUT)                  :: iy
      INTEGER*4, INTENT(IN OUT)                  :: iz
      INTEGER*4                                  :: i
      INTEGER*4                                  :: ndiv,n1

      REAL*4 :: unif(1), df

      CALL as183(ix,iy,iz,1,unif)
      ndiv = nhi - nlow + 1
      df = 1.0/FLOAT(ndiv)
      i = INT(unif(1)/df)
      n = nlow + i
!      n=unif(1)
      RETURN

      END SUBROUTINE iungen
!***********************************************************************

      SUBROUTINE iungen1(n,nlow,nhi,ix,iy,iz)


      INTEGER*4, INTENT(OUT)                     :: n
      INTEGER*4, INTENT(IN)                      :: nlow
      INTEGER*4, INTENT(IN)                      :: nhi
      INTEGER*4, INTENT(IN OUT)                  :: ix
      INTEGER*4, INTENT(IN OUT)                  :: iy
      INTEGER*4, INTENT(IN OUT)                  :: iz
      INTEGER*4                                  :: i
      INTEGER*4                                  :: ndiv

      REAL*4 :: unif(1), df

      CALL as183(ix,iy,iz,1,unif)
      ndiv = nhi - nlow + 1
      df = 1.0/FLOAT(ndiv)
      i = INT(unif(1)/df)
      n = nlow + i

      RETURN

      END SUBROUTINE iungen1


!***********************************************************************

!     Subroutine FINDABC, given coordinates of three points in a plane,
!     calculates the coefficients A, B, and C which describe the
!     plane and places them in the corresponding arrays at the
!     locations specified by the index jj in the passing argument.
!***********************************************************************

      SUBROUTINE findabc(tript,a,b,c)


      REAL*4, INTENT(IN)                         :: tript(3,3)
      REAL*4, INTENT(OUT)                        :: a
      REAL*4, INTENT(OUT)                        :: b
      REAL*4, INTENT(OUT)                        :: c
      REAL*4                                     :: det
      REAL*4                                     :: trimod(3,3)

!  Initialize vars

      a = 0.0
      b = 0.0
      c = 0.0
      det = 0.0

      DO  k = 1, 3
        trimod(k,1) = tript(k,1)
        trimod(k,2) = tript(k,2)
        trimod(k,3) = 1.0
      END DO
      CALL det3d(trimod,det)
      c = det
      DO  k = 1, 3
        trimod(k,1) = tript(k,3)
        trimod(k,2) = tript(k,1)
      END DO
      CALL det3d(trimod,det)
      b = det
      DO  k = 1, 3
        trimod(k,1) = tript(k,2)
        trimod(k,2) = tript(k,3)
      END DO
      CALL det3d(trimod,det)
      a = det

      RETURN
      END SUBROUTINE findabc
!***********************************************************************

!     Subroutine DET3D calculates the determinant of a three-dimensional
!        array. (see Tuma, 1987, p.8)

!***********************************************************************

      SUBROUTINE det3d(a,det)

      REAL*4, INTENT(IN)                         :: a(3,3)
      REAL*4, INTENT(OUT)                        :: det
      REAL*8                                     :: cof11, cof21, cof31

! Initialize vars

      cof11 = 0.d0
      cof21 = 0.d0
      cof31 = 0.d0
      det = 0.0

! Calculate cofactors
      cof11 = a(2,2)*a(3,3) - a(3,2)*a(2,3)
      cof21 = -(a(1,2)*a(3,3) - a(3,2)*a(1,3))
      cof31 = a(1,2)*a(2,3) - a(2,2)*a(1,3)

! Calculate determinant

      det = a(1,1)*cof11 + a(2,1)*cof21 + a(3,1)*cof31

      RETURN
      END SUBROUTINE det3d

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

!***********************************************************************
!  Function check checks the status of error flags
!***********************************************************************

      SUBROUTINE check(adum,ierr)

      INTEGER*4 :: ierr
      CHARACTER (LEN=15) :: adum

      IF (ierr /= 0) THEN
        WRITE(*,*) 'Error in allocating array or MPI call:  ',trim(adum)
        CALL MPI_Finalize(ierr)
        STOP
      ENDIF

      END SUBROUTINE check

!***********************************************************************

! Function XYDIST calculates the distance in cartesian coordinates
! between two points  (x1,y1) and (x2,y2).

!***********************************************************************

      FUNCTION xydist(x1,y1,x2,y2)



      REAL*4, INTENT(IN OUT)                            :: x1
      REAL*4, INTENT(IN OUT)                            :: y1
      REAL*4, INTENT(IN OUT)                            :: x2
      REAL*4, INTENT(IN OUT)                            :: y2
      REAL*4 :: xydist, arg

      arg = (ABS(x1-x2))**2.0 + (ABS(y1-y2))**2.0
      xydist = SQRT(arg)

      RETURN
      END FUNCTION xydist

!***********************************************************************

! Function GetTim calculates current time. Naum

!***********************************************************************
      
      SUBROUTINE gettim(ihr,imin,isec,i100th) 
    
        integer(2), intent(out):: ihr, imin, isec, i100th 
        character(8):: sdate 
        character(10):: stime 
        call date_and_time(sdate,stime) 
        read(sTime,"(I2,I2,I2,1x,I3)") ihr, imin, isec, i100th 
     
     END SUBROUTINE GetTim 
