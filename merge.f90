!  This is the main calling program for MERGE.
!  Ramya Ramanathan, Pacific Northwest Natioanl Lab/Wright State !  University. Last modified 2012/04/17.
!**********************************************************

      MODULE prcmpi

        INTEGER*4                                       ::  ierr
        INTEGER*4                                       ::  nproc
        INTEGER*4                                       ::  iproc
        INTEGER*4                                       ::  proc_num
        INTEGER*4                                       ::  dims(0:2)
        INTEGER*4                                       ::  myX
        INTEGER*4                                       ::  myY
        INTEGER*4                                       ::  myZ      
        INTEGER*4                                       ::  lxproc = 1 !Naum
  !     INTEGER*4                                       ::  lyproc !71 !Naum
        INTEGER*4                                       ::  lyproc
        INTEGER*4                                       ::  lzproc = 1
        INTEGER*4                                       ::  comm3d
        INTEGER*4                                       ::  coords(0:2)
  
        INTEGER*4                                       :: localIndices
  
        INTEGER*4                                       ::  upProc
        INTEGER*4                                       ::  downProc
        INTEGER*4                                       ::  frontProc
        INTEGER*4                                       ::  backProc
        INTEGER*4                                       ::  leftProc
        INTEGER*4                                       ::  rightProc
  
        INTEGER*4                                       ::  nodeLeft
        INTEGER*4                                       ::  nodeRight
        INTEGER*4                                       ::  nodeFront
        INTEGER*4                                       ::  nodeBack
        INTEGER*4                                       ::  nodeUp
        INTEGER*4                                       ::  nodeDown
  
        INTEGER*4                                       ::  lnx
        INTEGER*4                                       ::  lny
        INTEGER*4                                       ::  lnz
  
        INTEGER*4                                       ::  ibin
        INTEGER*4                                       ::  itest
        INTEGER*4                                       ::  ndims
        INTEGER*2                                       :: iii !Naum
  
        LOGICAL                                         ::  periods(0:2)
        LOGICAL                                         ::  reorder
  
        INTEGER*4, ALLOCATABLE                          ::  nodeFront1(:)
        INTEGER*4, ALLOCATABLE                          ::  nodeBack1(:)
  
        INTEGER*4                                       ::  ncbar_partition
        INTEGER*4, ALLOCATABLE                          :: ncbar_process(:)
  
        END MODULE prcmpi
  
  !**********************************************************
  
        MODULE tpix
  
        REAL*4, ALLOCATABLE                               :: slope(:)
        REAL*4, ALLOCATABLE                               :: yint(:)
        REAL*4, ALLOCATABLE                               :: xmin(:)
        REAL*4, ALLOCATABLE                               :: xmax(:)
        REAL*4, ALLOCATABLE                               :: cymin(:)
        REAL*4, ALLOCATABLE                               :: czmax(:)
        REAL*4, ALLOCATABLE                               :: twidth(:)
        REAL*4, ALLOCATABLE                               :: pslope(:)
        REAL*4, ALLOCATABLE                               :: pyint(:)
        REAL*4, ALLOCATABLE                               :: ppdn(:)
        REAL*4, ALLOCATABLE                               :: llen(:)
        REAL*4, ALLOCATABLE                               :: delh(:)
  
  
        !REAL*4                                            :: delhm !Naum never used
        !REAL*4                                            :: twidthm !Naum never used
        !REAL*4                                            :: llenm !Naum never used
  
        INTEGER*4, ALLOCATABLE                            :: lcnum(:,:,:) !Naum
        INTEGER*4, ALLOCATABLE                            :: ncoset(:)
        INTEGER*4, ALLOCATABLE                            :: nline(:,:)
        INTEGER*4, ALLOCATABLE                            :: ind(:)
        END MODULE tpix
  
  !***********************************************************************
  
        PROGRAM MAIN
  
        USE prcmpi
  
        INCLUDE 'mpif.h'
  
        INTEGER*4                                          :: nn
        INTEGER*4                                          :: nx
        INTEGER*4                                          :: ny
        INTEGER*4                                          :: nz
        INTEGER*4                                          :: ans
                                                           
        REAL*4                                             :: xulc
        REAL*4                                             :: yulc
        REAL*4                                             :: zulc
        REAL*4                                             :: xdim
        REAL*4                                             :: ydim
        REAL*4                                             :: zdim
        REAL*4                                             :: xspdim
        REAL*4                                             :: yspdim
        REAL*4                                             :: zspdim
      
        INTEGER*4                                          :: jans       
        INTEGER*4                                          :: lans       
        CHARACTER (LEN=15)                                 :: adum
  
        INTEGER*4                                          :: i_proc
  
        INTEGER*4                                          :: nocb
        INTEGER*4                                          :: numGridsY
        INTEGER*4                                          :: i_cbs
  
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
          lyproc = nproc ! Naum 06.07.13 
  
        IF (iproc == 0 ) THEN
          CALL READ_GRID(ans,xulc,yulc,zulc,xdim,ydim,zdim,xspdim, &
                         yspdim,zspdim,nx,ny,nz)
          OPEN(55,FILE='maxv.out',STATUS='old')
          READ(55,*)nocb
          !nc_determine = MOD(nocb, nproc)
          !IF (nc_determine == 0) THEN
          ! ncbar_partition = nocb/nproc
          !ELSE
          !    ncbar_partition = nocb / nproc + 1
          ! ENDIF
          ! WRITE(*, *) 'The number for component bars ', nocb
          ! WRITE(*, *) 'The number for each cores is ', ncbar_partition   
          CLOSE(55)       
        ENDIF
  
        CALL MPI_Bcast(ans,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        adum = 'ans-mpi'
        CALL CHECK(adum,ierr)
  
       CALL MPI_Bcast(nx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        adum = 'nx-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(ny,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        adum = 'ny-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(nz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        adum = 'nz-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_BCAST(xulc,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'xulc-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_BCAST(yulc,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'yulc-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_BCAST(zulc,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'zulc-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_BCAST(xdim,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'xdim-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_BCAST(ydim,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'ydim-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_BCAST(zdim,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'zdim-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_BCAST(xspdim,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'xspdim-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_BCAST(yspdim,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'yspdim-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_BCAST(zspdim,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
        adum = 'zspdim-mpi'
        CALL CHECK(adum,ierr)
  
        nn = nx*ny*nz
  
        dims(0) = lxproc
        dims(1) = lyproc
        dims(2) = lzproc
        ndims = 3
  
  ! Establish a cartesian topology communicator
  
        periods = .false.
        reorder = .false.
  
        CALL MPI_Cart_Create(MPI_COMM_WORLD,ndims,dims,periods, &
                             reorder,comm3d,ierr)
        adum = 'mpiCartCrt'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Cart_Coords(comm3d,iproc,ndims,coords,ierr)
        adum = 'mpiCartCd'
        CALL CHECK(adum,ierr)
  
  
        myX = coords(0)
        myY = coords(1)
        myZ = coords(2)
    !      write(*,*) 'line 206',myX,myY,myZ !Naum
  ! Returns the left and right neighbors of the Cartesian map 
  
        CALL MPI_Cart_shift(comm3d,0,1,leftProc,rightProc,ierr)
        adum = 'mpishift-lr'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Cart_shift(comm3d,1,1,frontProc,backProc,ierr)
        adum = 'mpishift-fb'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Cart_shift(comm3d,2,1,downProc,upProc,ierr)
        adum = 'mpishift-du'
        CALL CHECK(adum,ierr)
  
  ! Determine the subdomain this processor is solving
  
        ALLOCATE (nodeFront1(nproc),STAT=IERR)
        ALLOCATE (nodeBack1(nproc),STAT=IERR)
        WRITE(*,*) 'The y-axis coordinate is ', myY
        Do i_proc = 0, nproc-1 
          nodeFront1(i_proc+1)=(ny/lyproc)*i_proc+1+min(i_proc,mod(ny,lyproc))
          nodeBack1(i_proc+1)=(ny/lyproc)*(i_proc+1)+min(i_proc+1,mod(ny,lyproc))
          ! WRITE(*,*) 'The starting node is ',  nodeFront1(i_proc+1)
          ! WRITE(*,*) 'The ending node is ', nodeBack1(i_proc+1)
        ENDDO
        localIndices = myY+1
        ! WRITE(*,*) 'The local staring node is ', nodeFront1(localIndices)
  
        nodeLeft= (nx/lxproc)*myX+1+min(myX,mod(nx,lxproc))
        nodeRight = (nx/lxproc)*(myX+1)+min(myX+1,mod(nx,lxproc))
  
        nodeFront = (ny/lyproc)*myY+1+min(myY,mod(ny,lyproc))
        nodeBack = (ny/lyproc)*(myY+1)+min(myY+1,mod(ny,lyproc))
  
        nodeUp = (nz/lzproc)*myZ+1+min(myZ,mod(nz,lzproc))
        nodeDown = (nz/lzproc)*(myZ+1)+min(myZ+1,mod(nz,lzproc))
        
        lnx = nodeRight-nodeLeft+1
        lny = nodeBack-nodeFront+1
        lnz = nodeDown-nodeUp+1
  
        CALL MPI_Barrier(comm3d,ierr)
        adum = 'mpibarrier'
        CALL CHECK(adum,ierr)
        !IF(iproc==0) WRITE(*,*)'DO you need merge? Enter 1(yes) or 0(no)'
        !IF(iproc==0) READ(*,*)jans
        !IF(iproc==0)WRITE(5555,*)jans
        jans = 1
  
        !CALL MPI_BCAST(jans,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        !adum = 'jans-mpi'
        !CALL CHECK(adum,ierr)
        IF(jans==0 ) GOTO 599 
        CALL mergeAll(ans,xulc,yulc,zulc,xdim,ydim,zdim,xspdim, &
                      yspdim,zspdim,nx,ny,nz)
  
        599 CONTINUE
  
  !      CALL permgen(nx,ny,nz)
  
  !      DEALLOCATE (nodeFront1,STAT=IERR)
  !      DEALLOCATE (nodeBack1,STAT=IERR)
        CALL MPI_FINALIZE(ierr)
        adum = 'mpiFinalize'
        CALL CHECK(adum,ierr)
  
        END
  
  !**********************************************************
  ! 
  ! Code converted using TO_F90 by Alan Miller
  ! Date: 2007-07-29  Time: 16:27:20
  
  !    Program merge creates a file which digitally represents a
  !    compound bar with 'n' number of unit bars and 'n' # of trough sets:
  !    Each UB requires seperate arc information and plane equations
  !    Indicators are given to each UB. Each trough set requires separate
  !    line information and plane equations. Requires tpix.inc
  
  !    Input files
  !    1) lenstat.out-contains line parameters. (output from linegen.for)
  !    2) xcord.out: contains  xmin, xmax for each trough (output from linegen.for)
  !    3) angle.out: cotnains angle of each coset (output from linegen.for)
  !    4) indicator.out: contains info. about texture assgined to each trough.
  !       (output from linegen.for)
  !    5) location.out: cotnains location of a trough with respect to other 
  !       troughs. (output from linegen.for)
  !    6) coset.out: contains # of cosets in each unit bar. (output from linegen.for)
  !    7) line.out: contains # of lines in each coset.(output from linegen.for)
  !    8) out1.out: contains coefficients of plane equations for trough sets
  
  !      Output files
  !    1)pixdat.out: contains x,y, z coordinate values and the indicator of that point.
  !      This file is required to visualize the units that are simulated and for assiging permeability.
  !    2) results.out: contains volume proportion of each texture that are assigned.
  
  !    Indicators used:
  !    0-outside compound bar
  !    1-inside the bottom bed region of compound bar
  !    2-outside unit bars but inside compound bars
  !    3-inside trough set having fine sand as textural material
  !    4-inside trough set having open framework gravel as textural material
  !    5-inside trough set having sandy gravel as textural material
  !    6-unfilled space inside unit bar
  !    7-inside unit bar but outside all trough sets
  
  
  !    6/26/90 Tim Scheibe
  !    Modified last 07/14/07
  
  !***********************************************************************
        SUBROUTINE mergeAll(ans,xulc,yulc,zulc,xdim,ydim,zdim,xspdim, &
                            yspdim,zspdim,nx,ny,nz)
  
        USE tpix
        USE prcmpi
        include 'mpif.h'
  
  !     VARIABLE DECLARATIONS
  
        REAL*4, ALLOCATABLE                               :: a(:,:)
        REAL*4, ALLOCATABLE                               :: b(:,:)
        REAL*4, ALLOCATABLE                               :: c(:,:)
        REAL*4, ALLOCATABLE                               :: d(:,:)
  
        REAL*4, ALLOCATABLE                               :: ba(:,:)
        REAL*4, ALLOCATABLE                               :: bb(:,:)
        REAL*4, ALLOCATABLE                               :: bc(:,:)
        REAL*4, ALLOCATABLE                               :: bd(:,:)
  
        REAL*4, ALLOCATABLE                               :: ha(:,:)
        REAL*4, ALLOCATABLE                               :: hb(:,:)
        REAL*4, ALLOCATABLE                               :: hc(:,:)
        REAL*4, ALLOCATABLE                               :: hd(:,:)
  
        REAL*4, ALLOCATABLE                               :: oa(:,:)
        REAL*4, ALLOCATABLE                               :: ob(:,:)
        REAL*4, ALLOCATABLE                               :: oc(:,:)
        REAL*4, ALLOCATABLE                               :: od(:,:)
  
        REAL*4, ALLOCATABLE                               :: a1(:,:)
        REAL*4, ALLOCATABLE                               :: b1(:,:)
        REAL*4, ALLOCATABLE                               :: c1(:,:)
        REAL*4, ALLOCATABLE                               :: d1(:,:)
  
        REAL*4, ALLOCATABLE                               :: ba1(:,:)
        REAL*4, ALLOCATABLE                               :: bb1(:,:)
        REAL*4, ALLOCATABLE                               :: bc1(:,:)
        REAL*4, ALLOCATABLE                               :: bd1(:,:)
  
        REAL*4, ALLOCATABLE                               :: ha1(:,:)
        REAL*4, ALLOCATABLE                               :: hb1(:,:)
        REAL*4, ALLOCATABLE                               :: hc1(:,:)
        REAL*4, ALLOCATABLE                               :: hd1(:,:)
  
        REAL*4, ALLOCATABLE                               :: oa1(:,:)
        REAL*4, ALLOCATABLE                               :: ob1(:,:)
        REAL*4, ALLOCATABLE                               :: oc1(:,:)
        REAL*4, ALLOCATABLE                               :: od1(:,:)
  
        REAL*4, ALLOCATABLE                               :: xc(:,:)
        REAL*4, ALLOCATABLE                               :: yc(:,:)
        REAL*4, ALLOCATABLE                               :: xinfl(:,:)
        REAL*4, ALLOCATABLE                               :: yinfl(:,:)
        REAL*4, ALLOCATABLE                               :: r(:,:)
        REAL*4, ALLOCATABLE                               :: alpha(:,:)
        REAL*4, ALLOCATABLE                               :: larc(:,:)
        REAL*4, ALLOCATABLE                               :: lsum(:,:)
        REAL*4, ALLOCATABLE                               :: xx(:,:)
        REAL*4, ALLOCATABLE                               :: width(:)
  ! 06-18-08      
        REAL*4, ALLOCATABLE                               :: length(:)
        REAL*4, ALLOCATABLE                               :: height(:)
  !   
        REAL*4, ALLOCATABLE                               :: xpr(:)
        REAL*4, ALLOCATABLE                               :: ypr(:)
        REAL*4                               :: xprc
        REAL*4                               :: yprc
  
        REAL*4                                            :: ca(15)
        REAL*4                                            :: cb(15)
        REAL*4                                            :: cc(15)
        REAL*4                                            :: cd(15)
  
        REAL*4                                            :: ca1(15)
        REAL*4                                            :: cb1(15)
        REAL*4                                            :: cc1(15)
        REAL*4                                            :: cd1(15)
  
        REAL*4                                            :: mla(3)
        REAL*4                                            :: mlb(3)
        REAL*4                                            :: mlc(3)
        REAL*4                                            :: mld(3)
  
        REAL*4                                            :: mra(3)
        REAL*4                                            :: mrb(3)
        REAL*4                                            :: mrc(3)
        REAL*4                                            :: mrd(3)
  
        REAL*4, ALLOCATABLE                               :: ta(:)
        REAL*4, ALLOCATABLE                               :: tb(:)
        REAL*4, ALLOCATABLE                               :: tc(:)
        REAL*4, ALLOCATABLE                               :: td(:)
  
        REAL*4, ALLOCATABLE                               :: pta(:)
        REAL*4, ALLOCATABLE                               :: ptb(:)
        REAL*4, ALLOCATABLE                               :: ptc(:)
        REAL*4, ALLOCATABLE                               :: ptd(:)
  
        REAL*4                                            :: x1
        REAL*4                                            :: y1
        REAL*4                                            :: widthc
        REAL*4                                            :: lengthc
        REAL*4                                            :: alphacb
        !REAL*4                                            :: indc1 ! Naum Should be integer
        integer*4                                         :: indc1 ! Naum Should be integer
        integer*2                                         :: indc1_t !Naum 
  
        REAL*4, ALLOCATABLE                               :: xcc(:)
        REAL*4, ALLOCATABLE                               :: ycc(:)
        REAL*4, ALLOCATABLE                               :: xinflc(:)
        REAL*4, ALLOCATABLE                               :: yinflc(:)
        REAL*4, ALLOCATABLE                               :: rc(:)
        REAL*4, ALLOCATABLE                               :: alphacbc(:)
        REAL*4, ALLOCATABLE                               :: lsumc(:)
        REAL*4                               :: widthcb
  
        REAL*4, ALLOCATABLE                               :: cba(:)
        REAL*4, ALLOCATABLE                               :: cbb(:)
        REAL*4, ALLOCATABLE                               :: cbc(:)
        REAL*4, ALLOCATABLE                               :: cbd(:)
  
        REAL*4                               :: xprcb
        REAL*4                               :: yprcb
  
        REAL*4                                            :: x
        REAL*4                                            :: y
        REAL*4                                            :: z
        REAL*4                                            :: dsubx
        REAL*4                                            :: dsuby
        REAL*4                                            :: dsubz
  ! 06-18-08
        REAL*4                                            :: azpr
        REAL*4                                            :: aypr
        REAL*4                                            :: axpr
  
        REAL*4                                            :: zd
        REAL*4                                            :: yd
        REAL*4                                            :: xd
  !
  
        REAL*4                                            :: xulc
        REAL*4                                            :: yulc
        REAL*4                                            :: zulc
        REAL*4                                            :: xdim
        REAL*4                                            :: ydim
        REAL*4                                            :: zdim
        REAL*4                                            :: xspdim
        REAL*4                                            :: yspdim
        REAL*4                                            :: zspdim
        REAL*4                                            :: zmax
  
        INTEGER*4                                         :: nx
        INTEGER*4                                         :: ny
        INTEGER*4                                         :: nz
  
        INTEGER*4                                         :: nxpp
        INTEGER*4                                         :: nypp
        INTEGER*4                                         :: nzpp
  
        INTEGER*4                                         :: ans
        INTEGER*4                                         :: i
        INTEGER*4                                         :: j
        INTEGER*4                                         :: k
        INTEGER*4                                         :: m
        !INTEGER*4                                         :: jk
        INTEGER*4                                         :: ij
        INTEGER*4                                         :: kk
        INTEGER*4                                         :: nq1
        INTEGER*4                                         :: narc1
        INTEGER*4                                         :: narcc1
        INTEGER*4                                         :: indts
        INTEGER*4                                         :: mcb
        INTEGER*4                                         :: mub
        INTEGER*4                                         :: mloc
        
        INTEGER*4                                         :: nocb
        INTEGER*4                                         :: ncbar
        INTEGER*4                                         :: nbc
        INTEGER*4                                         :: np
        INTEGER*4                                         :: nps
        INTEGER*4                                         :: npf
  
        INTEGER*4                                         :: ncrofg
        INTEGER*4                                         :: ncrfs
        INTEGER*4                                         :: ncrsg1
        INTEGER*4                                         :: ncrsg2
        INTEGER*4                                         :: ncrsg3
        INTEGER*4                                         :: nofgc
        INTEGER*4                                         :: nsgc1
        INTEGER*4                                         :: nsgc2
        INTEGER*4                                         :: nsgc3
        INTEGER*4                                         :: nctot
        INTEGER*4                                         :: nspc
        INTEGER*4                                         :: nsp
        INTEGER*4                                         :: ncitx
        INTEGER*4                                         :: nfsc
        INTEGER*4                                         :: nctotG
        INTEGER*4                                         :: nspcG
        INTEGER*4                                         :: ncitxG
        INTEGER*4                                         :: nfscG
        INTEGER*4                                         :: nofgcG
        INTEGER*4                                         :: nsgc1G
        INTEGER*4                                         :: nsgc2G
        INTEGER*4                                         :: nsgc3G
        INTEGER*4                                         :: ncrfsG
        INTEGER*4                                         :: ncrofgG
        INTEGER*4                                         :: ncrsg1G
        INTEGER*4                                         :: ncrsg2G
        INTEGER*4                                         :: ncrsg3G
  
        INTEGER*4, ALLOCATABLE                            :: iinf(:,:)
  
        INTEGER*4, ALLOCATABLE                            :: ierru(:)
        INTEGER*4                            :: ierrcb
        INTEGER*4                                         :: iflag 
        INTEGER*4                                         :: ierrc
  
  
        INTEGER*4, ALLOCATABLE                            :: numCB(:,:,:)
        INTEGER*4, ALLOCATABLE                            :: numUB(:,:,:)
        INTEGER*4, ALLOCATABLE                            :: icval(:,:,:)
        INTEGER*4, ALLOCATABLE                            :: lcval(:,:,:)
  
        INTEGER*4                                         :: nup
        INTEGER*4                                         :: nbp
        INTEGER*4                                         :: nhp
        INTEGER*4, ALLOCATABLE                            :: nop(:)
  
        INTEGER*4                                         :: noub
        INTEGER*4                                         :: notss
        INTEGER*4                                         :: nocs
        INTEGER*4                                         :: nocl
  ! 07-07-08
        
        INTEGER*4                                         ::na
        INTEGER*4                                         ::na1
  
        INTEGER*4                                         :: jjcb
        INTEGER*4                                         :: jjcbc
        INTEGER*4                                         :: iicb
        INTEGER*4                                         :: mlj
  !     INTEGER*4, ALLOCATABLE                            :: ind(:,:)
  !     INTEGER*4, ALLOCATABLE                            :: ind(:,:,:,:)
        INTEGER*4, ALLOCATABLE                            :: indc(:)
        
  ! 06-18-08
        INTEGER*4                                         :: nosf
        INTEGER*4                                         :: nosfG !Naum
        INTEGER*4                                         :: nnacb
        INTEGER*4                                         :: nbbcb
        INTEGER*4                                         :: nbbcbG !Naum
        INTEGER*4                                         :: nncb
        INTEGER*4                                         :: nncbG !Naum
        INTEGER*4                                         :: nopub
        INTEGER*4                                         :: nopubG
        INTEGER*4                                         :: nnub
        INTEGER*4                                         :: nnubG !Naum
  !     INTEGER*4                                         :: nspc
  !     INTEGER*4                                         :: nfsc
        INTEGER*4                                         :: nfsct
  !      INTEGER*4                                         :: nofgc
        INTEGER*4                                         :: nofgct
  !      INTEGER*4                                         :: nsgc
        INTEGER*4                                         :: nsgct1
        INTEGER*4                                         :: nsgct2
        INTEGER*4                                         :: nsgct3
  !      INTEGER*4                                         :: ncrfs
        INTEGER*4                                         :: ncrfst
  !      INTEGER*4                                         :: ncrofg
        INTEGER*4                                         :: ncrofgt
  !      INTEGER*4                                         :: ncrsg
        INTEGER*4                                         :: ncrsgt1
        INTEGER*4                                         :: ncrsgt2
        INTEGER*4                                         :: ncrsgt3
  !      INTEGER*4                                         :: nctot
        INTEGER*4                                         :: nctemp
  !      
  ! Declaration for space fill
  ! 06-18-08
        REAL*4, ALLOCATABLE                               :: minx(:)
        REAL*4, ALLOCATABLE                               :: miny(:)
        REAL*4, ALLOCATABLE                               :: minz(:)
        REAL*4, ALLOCATABLE                               :: maxx(:)
        REAL*4, ALLOCATABLE                               :: maxy(:)
        REAL*4, ALLOCATABLE                               :: maxz(:)
        REAL*4                                            :: mxlub
        REAL*4                                            :: wxub
        REAL*4                                            :: hxub
        REAL*4                                            :: mnx
        REAL*4                                            :: mny
        REAL*4                                            :: mnz
        REAL*4                                            :: mxx
        REAL*4                                            :: mxy
        REAL*4                                            :: mxz
        REAL*4                                            :: nupe
        REAL*4                                            :: nbpe
        REAL*4                                            :: nope
        REAL*4                                            :: nhpe 
        REAL*4                                            :: nspe
        REAL*4                                            :: nopubo
        REAL*4                                            :: nocube
        REAL*4                                            :: jcbe
        REAL*4                                            :: icbe
        REAL*4                                            :: icerr
  !      REAL*4                                            :: ifill
  !      REAL*4                                            :: ifub
        REAL*4                                            :: ifind
  !
  ! Variables for indicators
  ! 06-18-08
        INTEGER*4                                         :: indcf
        INTEGER*4                                         :: indcb
        INTEGER*4                                         :: indub
  !      INTEGER*4                                         :: indts
        INTEGER*4                                         :: iocb
        INTEGER*4                                         :: ioub
        INTEGER*4                                         :: iots
        INTEGER*4, ALLOCATABLE                            :: ibcval(:,:,:)
        INTEGER*4, ALLOCATABLE                            :: ibuval(:,:,:)
        INTEGER*4                                         :: indtsg1
        INTEGER*4                                         :: indtsg2
        INTEGER*4                                         :: indtsg3
        INTEGER*4                                         :: indtofg
        INTEGER*4                                         :: indcs
        INTEGER*4                                         :: indcofg
        INTEGER*4                                         :: indcsg1
        INTEGER*4                                         :: indcsg2
        INTEGER*4                                         :: indcsg3
        INTEGER*4                                         :: indus
        INTEGER*4                                         :: indusg1
        INTEGER*4                                         :: indusg2
        INTEGER*4                                         :: indusg3
        INTEGER*4                                         :: induofg
        INTEGER*4                                         :: indrcf
        INTEGER*4                                         :: indlcf
  
        INTEGER*4                                         :: mtind
  ! Space fill info
  !
        INTEGER*4                                         :: ifill
        INTEGER*4                                         :: ifub
  
        INTEGER*4                                         :: itrough
        INTEGER*4                                         :: tmub
  
  !
  ! Variables to generate specific unit bars
  !
        INTEGER*4                                         :: nsub
        INTEGER*4                                         :: neub
  !                
        CHARACTER (LEN=64)                                :: filename
        CHARACTER (LEN=4)                                 :: intform
        CHARACTER (LEN=15)                                :: adum
        DATA intform / '(I )' /
  ! 06-30-10
        
        INTEGER*4                                         :: iy,ix,iz
        INTEGER*4                                         :: ii,jj
        INTEGER*4                                         :: narc,nnacbG
        
        INTEGER*4                                         :: nc_determine
  
        INTEGER*4                                         :: i_nb
        INTEGER*4                                         :: tmpNodeFront1
        INTEGER*4                                         :: tmpNodeBack1
        INTEGER*4                                         :: tmpNodeFront12
        INTEGER*4                                         :: tmpNodeBack12
        INTEGER*4                                         :: kk1
        INTEGER*4                                         :: tmpDiff
  
        INTEGER*4                                         :: tmpnctot
        INTEGER*4                                         :: tmpnosf
        INTEGER*4                                         :: tmpnnacb
        INTEGER*4                                         :: tmpnbbcb
        INTEGER*4                                         :: tmpnncb
        INTEGER*4                                         :: tmpnopub
        INTEGER*4                                         :: tmpnnub
        INTEGER*4                                         :: tmpnfsc
        INTEGER*4                                         :: tmpnofgc
        INTEGER*4                                         :: tmpnsgc
        INTEGER*4                                         :: tmpnsgc2
        INTEGER*4                                         :: tmpnsgc1
        INTEGER*4                                         :: tmpncrfs
        INTEGER*4                                         :: tmpncrofg
        INTEGER*4                                         :: tmpncrsg3
        INTEGER*4                                         :: tmpncrsg2
        INTEGER*4                                         :: tmpncrsg1
        INTEGER*4                                         :: tmpnspc
  
  !     FORMAT STATEMENTS
  
      1 FORMAT(1X,a)
     10 FORMAT(a12)
     15 FORMAT(1X,4(f15.5,2X),f15.3,2X,f15.9,2X, 2(f15.5,2X),i2,2X,f15.5)
     20 FORMAT(1X,10(f6.4,1X))
     21 FORMAT(80I2)
     25 FORMAT(1X,4(f15.7,3X))
     26 FORMAT(1X,4(f25.12,3X))
     31 FORMAT(f7.3,x,f7.3,x,f7.3,x,i2)
     32 FORMAT (x, i8)
     34 FORMAT(2(i4,2X))
     35 FORMAT(1X,i2)
     36 FORMAT(5(1X,f7.2),1X,f7.2)
     37 FORMAT(i4,2X,5(f13.2,1X),i4)
     41 FORMAT(f7.3,x,f7.3,x,f7.3,2(x,f7.3),x,i2)
     43 FORMAT(x,f5.2,x,f5.2,x,f5.2,x,i3)
     40 FORMAT(1X,i4,1X,i8,2X,i8,2X,i8,2X,i8)
     61 FORMAT(4(1X,i4))
     69 FORMAT(1X,i4,1X,i4,2X,i4,2X,i4,2X,i8)
     45 FORMAT(1X,i4,2X,i4,2X,i4)
     46 FORMAT(1X,i4,2X,i4,2X,i4,2X,i4)
     59 FORMAT(1X, a20, i9)
     94 FORMAT(1X,i4,2X,i4,2X,i4)
     95 FORMAT(1X,i4,1X,i4,2X,i4,2X,i4)
     96 FORMAT(1X,f10.4,2X,f10.4)
     97 FORMAT(1X,f10.4)
    100 FORMAT(1X,i4,1X,i4,1X,i4)
    101 FORMAT(1X,/)
    102 FORMAT(1X,i4,1X,i4)
    103 FORMAT(1X,i4,1X,i4,1X,i4)
    105 FORMAT(1X,i4,1X,i4,2X,i8,2X,i4)
    145 FORMAT(1X,i4,2X,i4,2X,i4)
  
  ! FILENAME READS AND OPEN STATEMENTS
  
        ibin = 1
        ! ibin=0
        itrough = 1
        ifill = 0
        ifub = 0
        tmub=1
        IF (iproc == 0) THEN
  
          WRITE(*,101)
          WRITE(*,101)
          WRITE(*,101) 
          WRITE(*,1) 'Program mergeAll'
          WRITE(*,1) '**************'
          WRITE(*,101)
  
  ! Always ascii files
  
          OPEN(8,FILE='maxv.out',STATUS='old')
          OPEN(16,FILE='results.out',STATUS='unknown')
  !        OPEN(20,FILE='indunits.out',STATUS='old')
          !OPEN(42,FILE='maxvals.out',STATUS='old')
          !iii = iproc
          !WRITE(intform(3:3),'(I1)') icount_digits(iii) !Naum maxvals.out never used
          !filename(1:12) = 'maxvals.out.'
          !WRITE(filename(13:),intform) iproc
          !OPEN(42,FILE=filename,FORM='formatted',STATUS='old') 
  
  ! Write data output file header
  
          READ(8,*) nocb
          CLOSE(8)
          ncbar_partition = nocb/nproc
  ! Ask if you need trough sets
  
        !Write(*,*)'Do you want to simulate trough sets? Enter 1 or 0'
        !READ(*,*)itrough    
        !IF(iproc==0)WRITE(5555,*)'itrough',itrough
        itrough = 1
  
  
  ! Space fill info.
  !
        !Write(*,*)'Is unfilled space filled with sand or donor sets?'
        !Write(*,*)'Enter 0 or 1'
        !Read(*,*)ifill
        !IF(iproc==0)WRITE(5555,*)'ifill',ifill
        ifill = 0
  
        !IF (ifill==0) THEN
        ! WRITE(*,*)'Enter the indicator for unfilled space'
        ! READ(*,*)ifind
        !IF(iproc==0)WRITE(5555,*)'ifind',ifind
        !ENDIF
        ifind = 0 !Naum
  
        !If(ifill==1) then
        !  Write(*,*)'Enter the unit bar number to be used as a donor'
        !  Read(*,*)ifub
        !IF(iproc==0)WRITE(5555,*)'ifub',ifub
        !Endif
  !
  !
  ! Read in maxvals
  
  
         !READ(42,*)delhm, twidthm, llenm       
         !CLOSE(42)
  
  ! Hard code indicator for each unit type
          indcf     = 17
          indrcf    = 18
          indlcf    = 19
          indcb     = 1
          indub     = 2
          indts     = 6
          indtofg   = 7
          indtsg1   = 22
          indtsg2   = 23
          indtsg3   = 8
          indcs     = 9
          indcofg   = 10
          indcsg1   = 24
          indcsg2   = 25
          indcsg3   = 11
          indus     = 12
          induofg   = 13
          indusg1   = 26
          indusg2   = 27
          indusg3   = 14
          iocb      = 0
          ioub      = 15
          iots      = 16
  
  ! Read in indicator for each unit type
  
  !        READ(20,*)indcf     !17
  !        READ(20,*)indrcf    !18
  !        READ(20,*)indlcf    !19
  !        READ(20,*)indcb     !1
  !        READ(20,*)indub     !2
  !        READ(20,*)indts     !6
  !        READ(20,*)indtofg   !7
  !        READ(20,*)indtsg1   !22
  !        READ(20,*)indtsg2   !23
  !        READ(20,*)indtsg3   !8
  !        READ(20,*)indcs     !9
  !        READ(20,*)indcofg   !10
  !        READ(20,*)indcsg1   !24
  !        READ(20,*)indcsg2   !25
  !        READ(20,*)indcsg3   !11
  !        READ(20,*)indus     !12
  !        READ(20,*)induofg   !13
  !        READ(20,*)indusg1   !26
  !        READ(20,*)indusg2   !27
  !        READ(20,*)indusg3   !14
  !        READ(20,*)iocb      !0
  !        READ(20,*)ioub      !15
  !        READ(20,*)iots      !16
  
  
  !       WRITE(202,*)indcf
  !       write(202,*)indrcf
  !       write(202,*)indlcf
  !       write(202,*)indcb
  !       write(202,*)indub
  !       write(202,*)indts
  !       write(202,*)indtofg
  !       write(202,*)indtsg1
  !       write(202,*)indtsg2
  !       write(202,*)indtsg3
  !       write(202,*)indcs
  !       write(202,*)indcofg
  !       write(202,*)indcsg1
  !       write(202,*)indcsg2
  !       write(202,*)indcsg3
  !       write(202,*)indus
  !       write(202,*)induofg
  !       write(202,*)indusg1
  !       write(202,*)indusg2
  !       write(202,*)indusg3
  !       write(202,*)iocb
  !       write(202,*)ioub
  !       write(202,*)iots
   
  !        CLOSE(20)
  !        CLOSE(202)
          WRITE(*,*)'Read Indicators for each unit type'
  
  
        ENDIF
  
        CALL MPI_Bcast(nocb,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'nocb-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(ncbar_partition,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        adum = 'ans-mpi'
        CALL CHECK(adum, ierr)
   
        CALL MPI_Bcast(itrough,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'itrough-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(ifill,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'ifill-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(ifind,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'ifind-mpi'
        CALL CHECK(adum,ierr)
  
      CALL MPI_Bcast(ifub,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'ifub-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(indcf,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indcf-mpi'
        CALL CHECK(adum,ierr)
        
        CALL MPI_Bcast(indrcf,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indrcf-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(indlcf,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indlcf-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(indcb,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indcb-mpi'
        CALL CHECK(adum,ierr)
        
        CALL MPI_Bcast(indub,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indub-mpi'
        CALL CHECK(adum,ierr)
       
        CALL MPI_Bcast(indts,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indts-mpi'
        CALL CHECK(adum,ierr)
        
        CALL MPI_Bcast(indtofg,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indtofg-mpi'
        CALL CHECK(adum,ierr)
        
        CALL MPI_Bcast(indtsg3,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indtsg3-mpi'
        CALL CHECK(adum,ierr)
        
        CALL MPI_Bcast(indtsg2,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indtsg2-mpi'
        CALL CHECK(adum,ierr)
       
        CALL MPI_Bcast(indtsg1,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indtsg1-mpi'
        CALL CHECK(adum,ierr)
        
        CALL MPI_Bcast(indcs,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indcs-mpi'
        CALL CHECK(adum,ierr)
       
        CALL MPI_Bcast(indcofg,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indcofg-mpi'
        CALL CHECK(adum,ierr)
   
        CALL MPI_Bcast(indcsg3,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indcsg3-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(indcsg2,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indcsg2-mpi'
        CALL CHECK(adum,ierr)
   
        CALL MPI_Bcast(indcsg1,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indcsg1-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(indus,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indus-mpi'
        CALL CHECK(adum,ierr)
   
        CALL MPI_Bcast(induofg,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'induofg-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(indusg3,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indusg3-mpi'
        CALL CHECK(adum,ierr)
   
        CALL MPI_Bcast(indusg2,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indusg2-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(indusg1,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'indusg1-mpi'
        CALL CHECK(adum,ierr)
   
        CALL MPI_Bcast(iocb,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'iocb-mpi'
        CALL CHECK(adum,ierr)
  
        CALL MPI_Bcast(ioub,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'ioub-mpi'
        CALL CHECK(adum,ierr)
   
        CALL MPI_Bcast(iots,1,MPI_INTEGER,0,comm3d,ierr)
        adum = 'iots-mpi'
        CALL CHECK(adum,ierr)
       
        CALL MPI_Bcast(tmub,1, MPI_INTEGER,0,comm3d,ierr)
        adum = 'tmub-mpi'
        CALL CHECK(adum,ierr)
   
        CALL MPI_Bcast(mtind,1, MPI_INTEGER,0,comm3d,ierr)
        adum = 'mtind-mpi'
        CALL CHECK(adum,ierr)
  
        !CALL MPI_Bcast(delhm,1, MPI_REAL,0,comm3d,ierr) !Naum never used
        !adum = 'delhm-mpi'
        !CALL CHECK(adum,ierr)
   
        !CALL MPI_Bcast(twidthm,1, MPI_REAL,0,comm3d,ierr)
        !adum = 'twidthm-mpi'
        !CALL CHECK(adum,ierr)
  
        !CALL MPI_Bcast(llenm,1, MPI_REAL,0,comm3d,ierr)
        !adum = 'llenm-mpi'
        !CALL CHECK(adum,ierr)
  
        
  ! Get array sizes
  
  !vlf
  !      OPEN(33, FILE='params.out',FORM='unformatted',ACCESS='direct', &
  !           RECORDTYPE='FIXED',RECL=16,STATUS='unknown')
  !vlf
   
  
  ! Allocate Array sized for compound bars; domain decomposition
  
  ! Allocate for cbbch.out read
        ALLOCATE (ncbar_process(lyproc),STAT=IERR)
        ! Establish a cartesian topology communicator
        nc_determine = MOD(nocb, lyproc)
        ! ncbar_partition = nocb/nproc
        ! ncbar_process(iproc) = ncbar_partition
        ! IF (nc_determine > 0) THEN
        ncbar_process(myY+1) = ncbar_partition
        IF (myY < nc_determine) THEN
            ncbar_process(myY+1) = ncbar_partition + 1
        ENDIF
        WRITE(*,*) 'Process here is ', myY
        WRITE(*,*) 'The number of component bars in this process is ', ncbar_process(myY+1)
        tmpNodeFront1 = nodeFront1(localIndices)
  
        nfsc = 0
        nfscG = 0
  
        nofgc = 0
        nofgcG = 0
  
        nsgc3 = 0
        nsgc3G = 0
  
        nsgc2 = 0
        nsgc2G = 0
  
        nsgc1 = 0
        nsgc1G = 0
  
        nctot = 0
        nctotG = 0
  
        ncitx = 0
        ncitxG = 0
  
        ncrofg = 0
        ncrofgG = 0
  
        ncrfs = 0
        ncrfsG = 0
  
        ncrsg3 = 0
        ncrsg3G = 0
  
        ncrsg2 = 0
        ncrsg2G = 0
  
        ncrsg1 = 0
        ncrsg1G = 0
  
        nsp = 0
  
        nspc = 0
        nspcG = 0
  
        nnub = 0
        nnubG = 0
  
        nosf = 0
        nosfG = 0
        ! nopubo = 0
        ! nocube = 0
        nncb = 0
        nncbG = 0
  
        nbbcb=0
        nbbcbG=0
  
        nopub=0
        nopubG=0
       
        nnacb=0
        nnacbG=0
  
        ! nsgct3 =0
        ! nsgct2 =0
        ! nsgct1 =0
  
        ! ncrsgt3 =0
        ! ncrsgt2 =0
        ! ncrsgt1 =0
  
        DO i_nb = 1,ncbar_process(myY+1)
          IF (myY<nc_determine) THEN
            kk1=myY*(ncbar_partition+1)+i_nb-1
          ELSE
            ! kk1=myY*ncbar_partition+i_nb-1
            kk1= nc_determine*(ncbar_partition+1)+&
                  (myY-nc_determine)*ncbar_partition+i_nb-1
          END IF
          tmpNodeFront12=(lny/ncbar_process(myY+1))*(i_nb-1)+1+&
              min(i_nb-1,mod(lny,ncbar_process(myY+1)))
          tmpNodeBack12=(lny/ncbar_process(myY+1))*(i_nb)+&
              min(i_nb,mod(lny,ncbar_process(myY+1)))
          tmpDiff=tmpNodeBack12-tmpNodeFront12
          tmpNodeBack1 = tmpNodeFront1+tmpDiff
          WRITE(*,*) 'The start node for a component bar is ', tmpNodeFront1
          WRITE(*,*) 'The end node for a component bar is ', tmpNodeBack1
  
          nbc = 5
  
          ALLOCATE(cba(nbc), STAT=IERR)
         adum = 'cba'
          CALL CHECK(adum,ierr)
  
          ALLOCATE(cbb(nbc), STAT=IERR)
          adum = 'cbb'
         CALL CHECK(adum,ierr)
  
        ALLOCATE(cbc(nbc), STAT=IERR)
        adum = 'cbc'
        CALL CHECK(adum,ierr)
  
        ALLOCATE(cbd(nbc), STAT=IERR)
        adum = 'cbd'
        CALL CHECK(adum,ierr)
  
  ! Allocate for cbarc.out read
  
        na = 3
        ALLOCATE(xcc(na), STAT=IERR)
        adum = 'xcc'
        CALL CHECK(adum,ierr)
        
        ALLOCATE(ycc(na), STAT=IERR)
        adum = 'ycc'
        CALL CHECK(adum,ierr)
        
        ALLOCATE(xinflc(na), STAT=IERR)
        adum = 'xinflc'
        CALL CHECK(adum,ierr)
        
        ALLOCATE(yinflc(na), STAT=IERR)
        adum = 'yinflc'
        CALL CHECK(adum,ierr)
        
        ALLOCATE(rc(na), STAT=IERR)
        adum = 'rc'
        CALL CHECK(adum,ierr)
        
        ALLOCATE(alphacbc(na), STAT=IERR)
        adum = 'alphacbc'
        CALL CHECK(adum,ierr)
        
        ALLOCATE(lsumc(na), STAT=IERR)
        adum = 'lsumc'
        CALL CHECK(adum,ierr)
        
  !      ALLOCATE(widthcb(ncbar), STAT=IERR)
  !      adum = 'widthcb'
  !      CALL CHECK(adum,ierr)
  
  ! Allocate for x-y-z loop
  
  !      ALLOCATE(xprc(ncbar), STAT=IERR)
  !      adum = 'xprc'
  !      CALL CHECK(adum,ierr)
  
  !      ALLOCATE(yprc(ncbar), STAT=IERR)
  !      adum = 'yprc'
  !      CALL CHECK(adum,ierr)
  
  !      ALLOCATE(xprcb(ncbar), STAT=IERR)
  !      adum = 'xprcb'
  !      CALL CHECK(adum,ierr)
  
  !      ALLOCATE(yprcb(ncbar), STAT=IERR)
  !      adum = 'yprcb'
  !      CALL CHECK(adum,ierr)
  
  !      ALLOCATE(ierrcb(ncbar), STAT=IERR)
  !      adum = 'ierrcb'
  !      CALL CHECK(adum,ierr)
  
      !   ALLOCATE(numCB(nodeLeft:nodeRight,nodeFront:nodeBack, &
      !            nodeUp:nodeDown))
          ! ALLOCATE(numCB(nodeLeft:nodeRight,&
          !         nodeFront1(localIndices):nodeBack1(localIndices), &
          !         nodeUp:nodeDown))
          ALLOCATE(numCB(nodeLeft:nodeRight,&
            tmpNodeFront1:tmpNodeBack1, &
            nodeUp:nodeDown))
        adum = 'numCB'
        numCB = 0
  
      !   ALLOCATE(numUB(nodeLeft:nodeRight,nodeFront:nodeBack, &
      !            nodeUp:nodeDown))
          ! ALLOCATE(numUB(nodeLeft:nodeRight,&
          !         nodeFront1(localIndices):nodeBack1(localIndices), &
          !         nodeUp:nodeDown))
          ALLOCATE(numUB(nodeLeft:nodeRight,&
            tmpNodeFront1:tmpNodeBack1, &
            nodeUp:nodeDown))
        adum = 'numUB'
        numUB = 0
  
      !   ALLOCATE(icval(nodeLeft:nodeRight,nodeFront:nodeBack, &
      !            nodeUp:nodeDown))
        ! ALLOCATE(icval(nodeLeft:nodeRight,&
        !           nodeFront1(localIndices):nodeBack1(localIndices), &
        !           nodeUp:nodeDown))
        ALLOCATE(icval(nodeLeft:nodeRight,&
          tmpNodeFront1:tmpNodeBack1, &
          nodeUp:nodeDown))
        adum = 'icval'
        icval = 0
  
      !   ALLOCATE(lcval(nodeLeft:nodeRight,nodeFront:nodeBack, &
      !            nodeUp:nodeDown))
        ! ALLOCATE(lcval(nodeLeft:nodeRight,&
        !           nodeFront1(localIndices):nodeBack1(localIndices), &
        !           nodeUp:nodeDown))
        ALLOCATE(lcval(nodeLeft:nodeRight,&
          tmpNodeFront1:tmpNodeBack1, &
          nodeUp:nodeDown))
        adum = 'lcval'
        lcval = 0
  
      !   ALLOCATE(ibcval(nodeLeft:nodeRight,nodeFront:nodeBack, &
      !            nodeUp:nodeDown))
          ! ALLOCATE(ibcval(nodeLeft:nodeRight,&
          !         nodeFront1(localIndices):nodeBack1(localIndices), &
          !         nodeUp:nodeDown))
        ALLOCATE(ibcval(nodeLeft:nodeRight,&
          tmpNodeFront1:tmpNodeBack1, &
          nodeUp:nodeDown))
       adum = 'ibcval'
        ibcval = 0
  
      !   ALLOCATE(ibuval(nodeLeft:nodeRight,nodeFront:nodeBack, &
      !            nodeUp:nodeDown))
        ! ALLOCATE(ibuval(nodeLeft:nodeRight,&
        !           nodeFront1(localIndices):nodeBack1(localIndices), &
        !           nodeUp:nodeDown))
        ALLOCATE(ibuval(nodeLeft:nodeRight,&
          tmpNodeFront1:tmpNodeBack1, &
          nodeUp:nodeDown))
        adum = 'ibuval'
        ibuval = 0
  
  !  Initialize for calculating volume proportion, # of cells simulated etc.
  
        
        nfsct = 0
  
        nofgct = 0
  
        nsgc3t = 0
  
        nsgc2t= 0
  
        nsgc1t = 0
  
        ! nctot = 0
  
  
        ncitx = 0
        ncitxG = 0
  
        ncrofgt = 0
  
        ncrfst = 0
  
        ncrsg3t = 0
  
        ncrsg2t = 0
  
        ncrsg1t = 0
  
        nsp = 0
  
        ! nspc = 0
        ! nspcG = 0
  
        ! nnub = 0
  
        ! nosf = 0
  
        nopubo = 0
        nocube = 0
  
        ! nncb = 0
  
        ! nbbcb=0
  
        ! nopub=0
       
        ! nnacb=0
  
        nsgct3 =0
        nsgct2 =0
        nsgct1 =0
  
        ncrsgt3 =0
        ncrsgt2 =0
        ncrsgt1 =0
  
        tmpnctot = 0
        tmpnosf = 0
        tmpnnacb = 0
        tmpnbbcb = 0
        tmpnncb = 0
        tmpnopub = 0
        tmpnnub = 0
        tmpnfsc = 0
        tmpnofgc = 0
        tmpnsgc = 0
        tmpnsgc2 = 0
        tmpnsgc1 = 0
        tmpncrfs = 0
        tmpncrofg = 0
        tmpncrsg3 = 0
        tmpncrsg2 = 0
        tmpncrsg1 = 0
        tmpnspc = 0
  
  ! Open output files pixdat.out.X and fillpix.out.X
  
        mcb = 0
        mub = 0
        ! iii = iproc
        iii=kk1
        tmpDiff=tmpDiff+1
        WRITE(intform(3:3),'(I1)') icount_digits(iii)
  
        filename(1:12) = 'pixdata.out.'
        WRITE(filename(13:),intform) kk1
        OPEN(2,FILE=filename,STATUS='unknown')
        ! WRITE(2,*) lnx,lny,lnz
        WRITE(2,*) lnx,tmpDiff,lnz
        WRITE(2,*) nodeLeft,nodeRight
        ! WRITE(2,*) nodeFront,nodeBack
        WRITE(2,*) tmpNodeFront1,tmpNodeBack1
        WRITE(2,*) nodeUp,nodeDown
        WRITE(2,*) nx, ny, nz
        CLOSE(2)
        
        IF (ibin == 1) THEN
          filename(1:11) = 'testin.out.'
          WRITE(filename(12:),intform) kk1 !iproc
          OPEN(6667,FILE=filename,FORM='formatted', &
               STATUS='unknown')
        
          filename(1:11) = 'pixdat.out.'
          WRITE(filename(12:),intform) kk1 !iproc
          OPEN(7,FILE=filename,FORM='unformatted', &
               STATUS='unknown')
          CLOSE(7, STATUS='delete')
          OPEN(7,FILE=filename,FORM='unformatted', &
               STATUS='new')
        ELSE
        
          filename(1:11) = 'pixdat.out.'
          WRITE(filename(12:),intform) kk1 !iproc
          OPEN(7,FILE=filename,FORM='formatted', &
               STATUS='unknown')
          CLOSE(7, STATUS='delete')
          OPEN(7,FILE=filename,FORM='formatted', &
               STATUS='new')
        ENDIF
  
  
        IF (iproc == 0) THEN
          WRITE(*,101)
          WRITE(*,1) 'Beginning Calculations ... '
          WRITE(*,101)
        ENDIF
  
  ! Parameters for array sizes
  
        nup = 17
        nbp = 2
        nhp = 1
        np = 2
  
  ! Read input files by compound bar
  
        DO  kk = 1, nocb
  
          !IF(iproc == 0)WRITE(*,*) 'Reading compound bar number: ',kk        
  
  ! Get array sizes before allocating arrays
        iii = kk-1
        WRITE(intform(3:3),'(I1)') icount_digits(iii)
         filename(1:11) = 'params.out.'
        WRITE(filename(12:),intform) kk-1
        OPEN(33, FILE=filename,FORM='unformatted',STATUS='unknown')        
         
  !        READ(33,REC=kk)noub,notss,nocs,nocl        
          READ(33)noub,notss,nocs,nocl        
                   
  !        filename(1:9) = 'wpar.out.'
  !        WRITE(filename(10:),intform) kk-1
  !        OPEN(450,FILE=filename,FORM='formatted',STATUS='unknown')
  !         WRITE(450,*) ': ',noub,notss,nocs,nocl
  !        CLOSE(450)
  
  ! Allocate for ubch.out read
          WRITE(*,*) 'The number of unit bars in the file is ', noub
          ALLOCATE(nop(noub), STAT=IERR)
          adum = 'nop'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(a(noub,nup), STAT=IERR)
          adum = 'a'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(b(noub,nup), STAT=IERR)
          adum = 'b'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(c(noub,nup), STAT=IERR)
          adum = 'c'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(d(noub,nup), STAT=IERR)
          adum = 'd'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(ba(noub,nbp), STAT=IERR)
          adum = 'ba'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(bb(noub,nbp), STAT=IERR)
          adum = 'bb'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(bc(noub,nbp), STAT=IERR)
          adum = 'bc'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(bd(noub,nbp), STAT=IERR)
          adum = 'bd'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(oa(noub,np), STAT=IERR)
          adum = 'oa'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(ob(noub,np), STAT=IERR)
          adum = 'ob'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(oc(noub,np), STAT=IERR)
          adum = 'oc'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(od(noub,np), STAT=IERR)
          adum = 'od'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(ha(noub,nhp), STAT=IERR)
          adum = 'ha'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(hb(noub,nhp), STAT=IERR)
          adum = 'hb'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(hc(noub,nhp), STAT=IERR)
          adum = 'hc'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(hd(noub,nhp), STAT=IERR)
          adum = 'hd'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(a1(noub,nup), STAT=IERR)
          adum = 'a1'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(b1(noub,nup), STAT=IERR)
          adum = 'b1'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(c1(noub,nup), STAT=IERR)
          adum = 'c1'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(d1(noub,nup), STAT=IERR)
          adum = 'd1'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(ba1(noub,nbp), STAT=IERR)
          adum = 'ba1'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(bb1(noub,nbp), STAT=IERR)
          adum = 'bb1'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(bc1(noub,nbp), STAT=IERR)
          adum = 'bc1'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(bd1(noub,nbp), STAT=IERR)
          adum = 'bd1'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(oa1(noub,np), STAT=IERR)
          adum = 'oa1'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(ob1(noub,np), STAT=IERR)
          adum = 'ob1'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(oc1(noub,np), STAT=IERR)
          adum = 'oc1'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(od1(noub,np), STAT=IERR)
          adum = 'od1'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(ha1(noub,nhp), STAT=IERR)
          adum = 'ha1'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(hb1(noub,nhp), STAT=IERR)
          adum = 'hb1'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(hc1(noub,nhp), STAT=IERR)
          adum = 'hc1'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(hd1(noub,nhp), STAT=IERR)
          adum = 'hd1'
          CALL CHECK(adum,ierr)
  
  ! Allocate for coset.out read
  
          ALLOCATE(ncoset(noub), STAT=IERR)
          adum = 'ncoset'
          CALL CHECK(adum,ierr)
  
  ! Allocate for line.out read
  
          ALLOCATE(nline(noub,nocs), STAT=IERR)
          adum = 'nline'
          CALL CHECK(adum,ierr)
  
  ! Allocate for location.out read  - this can be deleted??
  
          ALLOCATE(lcnum(noub,nocs,nocl), STAT=IERR) !Naum
          adum = 'lcnum'
          CALL CHECK(adum,ierr)
  
  ! Allocate for trgena.out read
  
          ALLOCATE(xc(noub,na), STAT=IERR)
          adum = 'xc'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(yc(noub,na), STAT=IERR)
          adum = 'yc'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(xinfl(noub,na), STAT=IERR)
          adum = 'xinfl'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(yinfl(noub,na), STAT=IERR)
          adum = 'yinfl'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(r(noub,na), STAT=IERR)
          adum = 'r'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(alpha(noub,na), STAT=IERR)
          adum = 'alpha'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(larc(noub,na), STAT=IERR)
          adum = 'larc'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(lsum(noub,na), STAT=IERR)
          adum = 'lsum'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(iinf(noub,na), STAT=IERR)
          adum = 'iinf'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(xx(noub,na), STAT=IERR)
          adum = 'xx'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(width(noub), STAT=IERR)
          adum = 'width'
          CALL CHECK(adum,ierr)
  
          ALLOCATE(length(noub), STAT=IERR)
          adum = 'length'
          CALL CHECK(adum,ierr)
  
          ALLOCATE(height(noub), STAT=IERR)
          adum = 'height'
          CALL CHECK(adum,ierr)
  
          ALLOCATE(minx(noub), STAT=IERR)
          adum = 'minx'
          CALL CHECK(adum,ierr)
  
          ALLOCATE(miny(noub), STAT=IERR)
          adum = 'miny'
          CALL CHECK(adum,ierr)
  
          ALLOCATE(minz(noub), STAT=IERR)
          adum = 'minz'
          CALL CHECK(adum,ierr)
  
          ALLOCATE(maxx(noub), STAT=IERR)
          adum = 'maxx'
          CALL CHECK(adum,ierr)
  
          ALLOCATE(maxy(noub), STAT=IERR)
          adum = 'maxy'
          CALL CHECK(adum,ierr)
  
          ALLOCATE(maxz(noub), STAT=IERR)
          adum = 'maxz'
          CALL CHECK(adum,ierr)
  
  ! Allocate for indicator.out read
  
          ALLOCATE(ind(notss), STAT=IERR)
          adum = 'ind'
          CALL CHECK(adum,ierr)
  
  ! Allocate for out1.out (trough plane coefficients) read
  
          nps = 16
          ALLOCATE(ta(nps), STAT=IERR)
          adum = 'ta'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(tb(nps), STAT=IERR)
          adum = 'tb'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(tc(nps), STAT=IERR)
          adum = 'tc'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(td(nps), STAT=IERR)
          adum = 'td'
          CALL CHECK(adum,ierr)
  
          npf = 15
          ALLOCATE(pta(npf), STAT=IERR)
          adum = 'pta'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(ptb(npf), STAT=IERR)
          adum = 'ptb'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(ptc(npf), STAT=IERR)
          adum = 'ptc'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(ptd(npf), STAT=IERR)
          adum = 'ptd'
          CALL CHECK(adum,ierr)
  
  ! Allocate for lenstat.out read
  
          ALLOCATE(slope(notss), STAT=IERR)
          adum = 'slope'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(yint(notss), STAT=IERR)
          adum = 'yint'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(xmin(notss), STAT=IERR)
          adum = 'xmin'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(xmax(notss), STAT=IERR)
          adum = 'xmax'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(cymin(notss), STAT=IERR)
          adum = 'cymin'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(czmax(notss), STAT=IERR)
          adum = 'czmax'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(twidth(notss), STAT=IERR)
          adum = 'twidth'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(ppdn(notss), STAT=IERR)
          adum = 'ppdn'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(pslope(notss), STAT=IERR)
          adum = 'pslope'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(pyint(notss), STAT=IERR)
          adum = 'pyint'
          CALL CHECK(adum,ierr)
  
          ALLOCATE(llen(notss), STAT=IERR)
          adum = 'llen'
          CALL CHECK(adum,ierr)
  
          ALLOCATE(delh(notss), STAT=IERR)
          adum = 'delh'
          CALL CHECK(adum,ierr)
  
  ! Allocate for x-y-z loop
  
          ALLOCATE(xpr(noub), STAT=IERR)
          adum = 'xpr'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(ypr(noub), STAT=IERR)
          adum = 'ypr'
          CALL CHECK(adum,ierr)
          
          ALLOCATE(ierru(noub), STAT=IERR)
          adum = 'ierru'
          CALL CHECK(adum,ierr)
  
  ! End allocation; open files for read
  
          ij = 1         
          
          iii = kk-1
          WRITE(intform(3:3),'(I1)') icount_digits(iii)
  
          filename(1:9) = 'ubch.out.'
          WRITE(filename(10:),intform) kk-1
          OPEN(3,FILE=filename,FORM='unformatted',STATUS='unknown')
  
          READ(3)nop(1:noub)
  
          READ(3)a(1:noub,1:nup)
          READ(3)b(1:noub,1:nup)
          READ(3)c(1:noub,1:nup)
          READ(3)d(1:noub,1:nup)
  
          READ(3)a1(1:noub,1:nup)
          READ(3)b1(1:noub,1:nup)
          READ(3)c1(1:noub,1:nup)
          READ(3)d1(1:noub,1:nup)
  
          READ(3)ba(1:noub,1:nbp)
          READ(3)bb(1:noub,1:nbp)
          READ(3)bc(1:noub,1:nbp)
          READ(3)bd(1:noub,1:nbp)
  
          READ(3)ba1(1:noub,1:nbp)
          READ(3)bb1(1:noub,1:nbp)
          READ(3)bc1(1:noub,1:nbp)
          READ(3)bd1(1:noub,1:nbp)
  
          READ(3)ha(1:noub,1:nhp)
          READ(3)hb(1:noub,1:nhp)
          READ(3)hc(1:noub,1:nhp)
          READ(3)hd(1:noub,1:nhp)
  
          READ(3)ha1(1:noub,1:nhp)
          READ(3)hb1(1:noub,1:nhp)
          READ(3)hc1(1:noub,1:nhp)
          READ(3)hd1(1:noub,1:nhp)
  
          READ(3)oa(1:noub,1:np)
          READ(3)ob(1:noub,1:np)
          READ(3)oc(1:noub,1:np)
          READ(3)od(1:noub,1:np)
  
          READ(3)oa1(1:noub,1:np)
          READ(3)ob1(1:noub,1:np)
          READ(3)oc1(1:noub,1:np)
          READ(3)od1(1:noub,1:np)
  
          CLOSE(3)
  
          IF(iproc == 0)WRITE(*,*) 'Read UB planes... '
          
  !        IF(iproc == 0)WRITE(9081,*) 'ncbar,iproc ', ncbar,iproc,kk
          filename(1:9) = 'cbch.out.'
          WRITE(filename(10:),intform) kk-1
          OPEN(34,FILE=filename,FORM='unformatted',STATUS='unknown')
  
          jjcb = 15
          iicb = 15
          mlj = 3
          READ(34)ca
          READ(34)cb
          READ(34)cc
          READ(34)cd
          READ(34)ca1
          READ(34)cb1
          READ(34)cc1
          READ(34)cd1
          READ(34)mla
          READ(34)mlb
          READ(34)mlc
          READ(34)mld
          READ(34)mra
          READ(34)mrb
          READ(34)mrc
          READ(34)mrd
  
          CLOSE(34)
  
          IF (iproc == 0)WRITE(*,*)'Read unit/compound bar data... '
          iii = kk-1
          WRITE(intform(3:3),'(I1)') icount_digits(iii)
          filename(1:11) = 'points.out.'
          WRITE(filename(12:),intform) kk-1
          OPEN(19,FILE=filename,FORM='unformatted',STATUS='unknown')
  
          READ(19)x1,y1,widthc,lengthc,alphacb,indc1_t
          indc1 = indc1_t !Naum
  
          CLOSE(19)
  
          IF (iproc == 0)WRITE(*,*)'Read CB data... '
  
  ! READ IN PLAN VIEW GEOMETRIC DATA
  
          filename(1:11) = 'trgena.out.'
          WRITE(filename(12:),intform) kk-1
          OPEN(4,FILE=filename,FORM='unformatted',STATUS='unknown')
  
  ! Need to check this file read
          narc = 2
          READ(4) xc(1:noub,1:2)
          READ(4) yc(1:noub,1:2)
          READ(4) xinfl(1:noub,1:2)
          READ(4) yinfl(1:noub,1:2)
          READ(4) r(1:noub,1:2)
          READ(4) alpha(1:noub,1:2)
          READ(4) larc(1:noub,1:2)
          READ(4) lsum(1:noub,1:2)
          READ(4) iinf(1:noub,1:2)
          READ(4) xx(1:noub,1:2)
  !        READ(4) width(1,1:noub)
  !        READ(4) length(1,1:noub)
  !        READ(4) height(1,1:noub)
          READ(4) minx(1:noub)
          READ(4) miny(1:noub)
          READ(4) minz(1:noub)
          READ(4) maxx(1:noub)
          READ(4) maxy(1:noub)
          READ(4) maxz(1:noub)
          READ(4) height(1:noub)
          width(1:noub)=maxx(1:noub)
          length(1:noub)=maxy(1:noub)
          CLOSE(4)
  
  
          IF (iproc == 0)WRITE(*,*) 'Read trgena data... '
  
  ! Read in CBC data
  
          filename(1:10) = 'cbarc.out.'
          WRITE(filename(11:),intform) kk-1
          OPEN(919,FILE=filename,FORM='unformatted',STATUS='unknown')
  
          READ(919)xcc
          READ(919)ycc
          READ(919)xinflc
          READ(919)yinflc
          READ(919)rc
          READ(919)alphacbc
          READ(919)lsumc
          READ(919)widthcb
  
          CLOSE(919)
  
          filename(1:10) = 'cbcch.out.'
          WRITE(filename(11:),intform) kk-1
          OPEN(92,FILE=filename,FORM='unformatted',STATUS='unknown')
  
          jjcbc = 5
          READ(92)cba
          READ(92)cbb
          READ(92)cbc
          READ(92)cbd        
          CLOSE(92)
          IF (iproc == 0) WRITE(*,*)'Finished read channel fills planes'
  
  !  First, the # of coset boundaries (ncoset) in each unit bar is read.
  !  Next # of centre lines (of trough sets) in each coset is read.
  !  Finally, the location of trough sets (lcnum) is read.
  
  
            filename(1:10) = 'coset.out.'
            WRITE(filename(11:),intform) kk-1
              OPEN(13,FILE=filename,FORM='unformatted',STATUS='old')
              filename(1:9) = 'line.out.'
              WRITE(filename(10:),intform) kk-1
              OPEN(14,FILE=filename,FORM='unformatted',STATUS='old')
  
             READ(13)ncoset(1:noub)
  !         write(511,*)ncoset(1:noub)
  !            WRITE(34,*)(ncoset(ij,j), j=1, noub)!((nline(ij,i,j),i=1,noub),j=1,nocs)
  !            write(5681,*)nocs,nocl            
              READ(14)nline(1:noub,1:nocs)  !((nline(ij,i,j),i=1,noub),j=1,ncoset(ij,i))
                  
              !DO j = 1, noub
              !DO i = 1, nocs
                 
             
              !enddo
              !enddo
  !            WRITE(34,*)((nline(ij,i,j),i=1,noub),j=1,ncoset(ij,i))
  !            CLOSE(34)
              IF (iproc == 0) WRITE(*,*)'read of coset and line'
  !        ELSE
  
  !          filename(1:10) = 'coset.out.'
  !          WRITE(filename(11:),intform) kk-1
  !          OPEN(13,FILE=filename,FORM='formatted',STATUS='old')
  !  
  !          filename(1:9) = 'line.out.'
  !          WRITE(filename(10:),intform) kk-1
  !          OPEN(14,FILE=filename,FORM='formatted',STATUS='old')
  !
  !          READ(13,*)(ncoset(ij,i),i=1,noub)
  !          READ(14,*)((nline(ij,i,j),i=1,noub),j=1,nocs)
  
  !        ENDIF
  
  !vlf
  !        nq1 = 0
  !        DO  i = 1,noub
  !          DO  j = 1,nocs !ncoset(i)
  !            DO  k = 1,nocl !nline(i,j)
  !              nq1 = nq1 + 1
  !              lcnum(i,j,k) = nq1
  !            if(lcnum(ij,i,j,k)==0) then
  !           write(538,*)lcnum(ij,i,j,k),ij,i,j,k
  !          endif
   
  !            END DO
  !          END DO
  !        END DO        
  !vlf
          CLOSE(13)
          CLOSE(14)
          CLOSE(11)
    
          IF (iproc == 0)WRITE(*,*)'Read location of trough sets... '
  
  !  Reads in indicator values (as texture) assigned for each trough
  
          filename(1:14) = 'indicator.out.'
          WRITE(filename(15:),intform) kk-1
  
  !        IF (ibin == 1) THEN
            OPEN(15,FILE=filename,FORM='unformatted',STATUS='old')
            READ(15)ind(1:notss)
  !         do j=1,notss
  !          if(ind(ij,j)==0) then
  !           write(537,*)ind(ij,j),j
  !          endif
  !         enddo
  !        ELSE
  !          OPEN(15,FILE=filename,FORM='formatted',STATUS='old')
  !          READ(15,*)(ind(ij,jj),jj=1,notss)
  !        ENDIF
       
          CLOSE(15)
  
          IF(iproc == 0)WRITE(*,*)'Read indicator values... '
  
  !  Read in trough plane coefficients
  
  
          filename(1:9) = 'OUT1.out.'
          WRITE(filename(10:),intform) kk-1
  
  !        IF (ibin == 1) THEN
            OPEN(17,FILE=filename,FORM='unformatted',STATUS='old')
  !          READ(17)ta(ij,1:notss,1:16)
  !          READ(17)tb(ij,1:notss,1:16)
  !          READ(17)tc(ij,1:notss,1:16)
  !          READ(17)td(ij,1:notss,1:16)
  !          READ(17)pta(ij,1:notss,1:15)
  !          READ(17)ptb(ij,1:notss,1:15)
  !          READ(17)ptc(ij,1:notss,1:15)
  !          READ(17)ptd(ij,1:notss,1:15)
                         
            READ(17)(ta(ii),ii=1,16)
  !          IF(iproc == 0)WRITE(*,*)'Read ta'
            READ(17)(tb(ii),ii=1,16)
  !          IF(iproc == 0)WRITE(*,*)'Read tb'
            READ(17)(tc(ii),ii=1,16)
  !          IF(iproc == 0)WRITE(*,*)'Read tc'
            READ(17)(td(ii),ii=1,16)
  !          IF(iproc == 0)WRITE(*,*)'Read td'
            READ(17)(pta(ii),ii=1,15)
  !          IF(iproc == 0)WRITE(*,*)'Read pta'
            READ(17)(ptb(ii),ii=1,15)
  !          IF(iproc == 0)WRITE(*,*)'Read ptb'
            READ(17)(ptc(ii),ii=1,15)
  !          IF(iproc == 0)WRITE(*,*)'Read ptc'
            READ(17)(ptd(ii),ii=1,15)
  !          IF(iproc == 0)WRITE(*,*)'Read ptd'
            
  !          write(1700,*)((ta(ij,jj,ii),jj=1,notss),ii=1,16)
  !          write(1700,*)((tb(ij,jj,ii),jj=1,notss),ii=1,16)
  !          write(1700,*)((tc(ij,jj,ii),jj=1,notss),ii=1,16)
  !          write(1700,*)((td(ij,jj,ii),jj=1,notss),ii=1,16)
  !          write(1700,*)((pta(ij,jj,ii),jj=1,notss),ii=1,15)
  !          write(1700,*)((ptb(ij,jj,ii),jj=1,notss),ii=1,15)
  !          write(1700,*)((ptc(ij,jj,ii),jj=1,notss),ii=1,15)
  !          write(1700,*)((ptd(ij,jj,ii),jj=1,notss),ii=1,15)
  
  !        ELSE
  !          OPEN(17,FILE=filename,FORM='formatted',STATUS='old')
  !          READ(17,*)((ta(ij,jj,ii),jj=1,notss),ii=1,16)
  !          READ(17,*)((tb(ij,jj,ii),jj=1,notss),ii=1,16)
  !          READ(17,*)((tc(ij,jj,ii),jj=1,notss),ii=1,16)
  !          READ(17,*)((td(ij,jj,ii),jj=1,notss),ii=1,16)
  !          READ(17,*)((pta(ij,jj,ii),jj=1,notss),ii=1,15)
  !          READ(17,*)((ptb(ij,jj,ii),jj=1,notss),ii=1,15)
  !          READ(17,*)((ptc(ij,jj,ii),jj=1,notss),ii=1,15)
  !          READ(17,*)((ptd(ij,jj,ii),jj=1,notss),ii=1,15)
  !        ENDIF
  
          CLOSE(17)
  
          IF(iproc == 0)WRITE(*,*)'Read trough plane coefficients... '
  
  !  Read in the input for line parameters
  
          filename(1:12) = 'lenstat.out.'
          WRITE(filename(13:),intform) kk-1
  
  !        IF (ibin == 1) THEN
            OPEN(9,FILE=filename,FORM='unformatted',STATUS='unknown')
            READ(9)slope(1:notss)
            READ(9)yint(1:notss)
            READ(9)xmin(1:notss)
            READ(9)xmax(1:notss)
            READ(9)cymin(1:notss)
            READ(9)czmax(1:notss)
            READ(9)twidth(1:notss)
            READ(9)ppdn(1:notss)
            READ(9)pslope(1:notss)
            READ(9)pyint(1:notss)
  !        ELSE
  !          OPEN(9,FILE=filename,FORM='formatted',STATUS='unknown')
  !          READ(9,*)(slope(ij,jj),jj = 1,notss)
  !          READ(9,*)(yint(ij,jj),jj = 1,notss)
  !          READ(9,*)(xmin(ij,jj),jj = 1,notss)
  !          READ(9,*)(xmax(ij,jj),jj = 1,notss)
  !          READ(9,*)(cymin(ij,jj),jj = 1,notss)
  !          READ(9,*)(czmax(ij,jj),jj = 1,notss)
  !          READ(9,*)(twidth(ij,jj),jj = 1,notss)
  !          READ(9,*)(ppdn(ij,jj),jj = 1,notss)
  !          READ(9,*)(pslope(ij,jj),jj = 1,notss)
  !          READ(9,*)(pyint(ij,jj),jj = 1,notss)
  !        ENDIF
           IF(iproc == 0)WRITE(*,*)'Read trough slope coefficients... '
           !WRITE(*,*)'line 1907:'
         nq1 = 0
         DO  i=1,noub
           DO  j=1,ncoset(i)
             DO  k=1,nline(i,j)
               nq1 = nq1+1
               lcnum(i,j,k) = nq1 !Naum
             END DO
           END DO
         END DO
         CLOSE(9)
  
           IF(iproc == 0)WRITE(*,*)'closed slope coeff file... '
  ! Read in the input for line parameters
  
          filename(1:11) = 'rdleng.out.'
          WRITE(filename(12:),intform) kk-1
  
  !        IF (ibin == 1) THEN
            OPEN(43,FILE=filename,FORM='unformatted',STATUS='unknown')
  !        filename(1:9) = 'chc1.out.'
  !        WRITE(filename(10:),intform) kk-1
  !          OPEN(93,FILE=filename,FORM='formatted',STATUS='unknown')
            READ(43)llen(1:notss)
  
  !          write(93,*)llen(ij,1:notss)
            READ(43)delh(1:notss)
  !          write(93,*)llen(ij,1:notss)
            CLOSE(43)
  !          CLOSE(93)
  
          IF (iproc == 0) THEN
              WRITE(*,*)'Read trough centre line parameters... '
              WRITE(*,*)' '
              WRITE(*,*)'Completed compound bar read: ',kk-1
              WRITE(*,*)' '
  !           WRITE(5671,*)'Completed compound bar read: ',kk,iproc
          ENDIF
  
  ! Grid calculations
  
        zmax=zulc+zspdim
  
  ! Will n*pp \= 1 ever be used in large scale simulations?
  
        IF (ans == 1) THEN
          nxpp = 1
          nypp = 1
          nzpp = 1
        ELSE
          WRITE(*,101)
          WRITE(*,1) 'Enter number of pixel subdivisions: '
          WRITE(*,101)
          WRITE(*,1) '      X - direction? '
          READ(*,*) nxpp
          WRITE(*,1) '      Y - direction? '
          READ(*,*) nypp
          WRITE(*,1) '      Z - direction? '
          READ(*,*) nzpp
          WRITE(*,101)
        END IF
  
        dsubx = xdim/FLOAT(nxpp)
        dsuby = ydim/FLOAT(nypp)
        dsubz = zdim/FLOAT(nzpp)
  
        ! nctot = 0
        tmpnctot=0
        tmpnpsc=0
        !write(*,*) iproc, 'X: ',nodeLeft,nodeRight, dsubx, dsuby, dsubz
        !write(*,*) iproc, 'Y: ',nodeFront,nodeBack
        !write(*,*) iproc, 'Z: ',nodeUp,nodeDown,nzpp*nodeDown
          ! DO  j = nodeFront, nodeBack
          ! write(*,*) 'In the calculation the starting node is ', nodeFront1(localIndices)
          ! DO  j = nodeFront1(localIndices),nodeBack1(localIndices)
          DO  j = tmpNodeFront1,tmpNodeBack1
            DO  i = nodeLeft, nodeRight
              x = xulc + FLOAT(i-1)*xdim + dsubx/2.0
              DO  ix = 1, nxpp
                y = yulc + FLOAT(j-1)*ydim + dsuby/2.0
                DO  iy = 1, nypp
                  z = zulc + dsubz/2.0
                  loopIZ:  DO  iz = nodeUp, nzpp*nodeDown
                  tmpnctot = tmpnctot+1
                  IF (numCB(i,j,iz) > 0) THEN
                    z = z + dsubz
                    CYCLE loopIZ
                  ENDIF
  
  !  Checking the point is inside the CB
  !  Axis Transformation
  !  Checks to see if inside plane equation for that compound bar
  
                    ierrc = 0
                    CALL cotranm(x,y,xprc,yprc,x1,y1,alphacb, &
                                 widthc,lengthc,indc1,ierrc)
              
  !
  ! ierrc=2 is outside compound bar
  !
                    IF(ierrc == 2)THEN
                      IF(kk==nocb) THEN
                           icval(i,j,iz) =iocb
                           tmpnosf = tmpnosf + 1
                           GOTO 171
                      ELSE
                       z = z + dsubz
                       CYCLE loopIZ
  !                    GO TO 170
                      END IF
                     ENDIF
  ! Loop over number of planes
  !
                    DO  jj = 1, jjcb
                      CALL above(xprc,yprc,z,ca(jj), &
                                 cb(jj),cc(jj),cd(jj),iflag)
  !
  ! not found; cycle, can be in the zone, but not within the compound bar
  !
                      IF (iflag == 1) THEN
                       IF(kk==nocb) THEN
                          icval(i,j,iz) =iocb
                          tmpnosf = tmpnosf + 1
                          GOTO 171
                       ENDIF
                        z = z + dsubz
                        CYCLE loopIZ
  !                     GO TO 170
                      END IF
                    END DO
                    numCB(i,j,iz) = kk
                    tmpnnacb = tmpnnacb+1
  !
  ! Bottom plane; boundary plane of compound bar; loop over number of 
  ! planes which defined bottom boundary
  !
                    DO  jj = 1, iicb
                      CALL above(xprc,yprc,z,ca1(jj), &
                                 cb1(jj),cc1(jj),cd1(jj),iflag)
  !                     
                      IF (iflag == 1) THEN
                        ibcval (i,j,iz) = 1
                        tmpnbbcb = tmpnbbcb+1
                        icval(i,j,iz)=indcb
                        GO TO 235
                        z = z + dsubz
                        CYCLE loopIZ             
                      END IF
                    END DO
  
                        icval(i,j,iz)=indcb
  ! Major Channel Fills(Left)
              
  ! Loop over number of planes that define major channel fills
  ! Two major channel fills, left and right - two loop nests
  
                    DO  jj=1,mlj
                      CALL above(xprc,yprc,z,mla(jj), &
                                 mlb(jj),mlc(jj),mld(jj),iflag)
                      IF(iflag == 0) THEN
                        icval(i,j,iz) = indlcf
  !                      z = z + dsubz
  !                      CYCLE loopIZ
                       GO TO 170
                      END IF
                    END DO
                    DO  jj=1,mlj
                      CALL above(xprc,yprc,z,mra(jj), &
                                 mrb(jj),mrc(jj),mrd(jj),iflag)
                      IF(iflag == 0) THEN
                        icval(i,j,iz) = indrcf
  !                      z = z + dsubz
  !                      CYCLE loopIZ
                       GO TO 170
                      END IF
                    END DO
                   
  ! Cross bar channel; narcc(ij) = number of cross bar channels; same as compound bar
  ! number because each compound bar has one cross bar channel.
  ! Inflection points of arcs; is it in cross bar channel for compound bar?
  
                    narcc1 = 2
                    ierrcb = 0.0
                    CALL cotran(xprc,yprc,xprcb,yprcb, &
                            narcc1,xcc,ycc,xinflc,yinflc,rc, &
                            alphacbc,lsumc,widthcb,ierrcb)
                    IF (ierrcb == 2) THEN
                      GO TO 235
                    END IF
  
  ! plane equation
  
                    DO  jj = 1, jjcbc
                      CALL above(xprcb,yprcb,z,cba(jj), &
                                 cbb(jj),cbc(jj),cbd(jj),iflag)
                      IF (iflag == 1) THEN
                        GO TO 235
                      END IF
                    END DO
                    icval(i,j,iz) = indcf
  !                  z = z + dsubz
  !                  CYCLE loopIZ
                   GO TO 170
    235             CONTINUE
                    tmpnncb = tmpnncb + 1            
                  
  ! Perform coordinate transform of (x,y)
  ! [If outside limits, set ivalue to 0]
  
                    
  ! Unit bar loop
                    narc1 = 1
                    na1=2
                    loop111:  DO  m = 1,noub
                      ierru(m) = 0
                      CALL cotran(xprc,yprc,xpr(m),ypr(m), &
                                  narc1,xc(m,1:na1),yc(m,1:na1), &
                                  xinfl(m,1:na1),yinfl(m,1:na1), &
                                  r(m,1:na1),alpha(m,1:na1), &
                                  lsum(m,1:na1),width(m),ierru(m))
                      IF (ierru(m) == 2) THEN
                        IF(m == noub)THEN
                          icval(i,j,iz) = ioub
                          tmpnopub = tmpnopub + 1
                          GOTO 171
  !                      ENDIF
  !                      if(iproc==0) write(*,*)'ierru'
  !                        z = z + dsubz
  !                        CYCLE loopIZ
  !                       GO TO 171
                        ELSE
  !                        if(iproc==0) write(*,*)'ierru!!!',m
                          CYCLE loop111
                        END IF
                      END IF
                      
  ! Checking the point is inside the  UB.  Point has to clear all four loops to be inside
  ! unit bar.
  !                    if (iproc==0) write (*,*) 'UB planes'
                      DO  jj = 1, nup
                        CALL above(xpr(m),ypr(m),z,a(m,jj), &
                                   b(m,jj),c(m,jj),d(m,jj),iflag)
                        IF (iflag == 0) THEN
  
  ! Read maximum number of unit bars?  If outside, will go to trough set
  
                          IF(m == noub)THEN
                            icval(i,j,iz) = ioub
                            mcb=kk
            !                lcval(i,j,iz)=icval(i,j,iz)
                            GO TO 171
                          ELSE
                            CYCLE loop111
                          END IF
                        END IF
                      END DO
  !                   DO  jj =1,nbp
  !                     CALL above(xpr(ij,m),ypr(ij,m),z,ba(m,jj), &
  !                                bb(m,jj),bc(m,jj),bd(m,jj), &
  !                                iflag)
  !                     IF(iflag == 1) THEN
  !                       IF(m == noub)THEN
  !                         icval(i,j,iz) = 1
  !                         lcval(i,j,iz)=icval(i,j,iz)
  !                         GO TO 289
  !                       ELSE
  !                         CYCLE loop111
  !                       END IF
  !                     END IF
  !                   END DO
                      DO  jj =1,nop(m)
                        CALL above(xpr(m),ypr(m),z,oa(m,jj), &
                                   ob(m,jj),oc(m,jj),od(m,jj), &
                                   iflag)
                        IF (iflag == 0) THEN
                          IF(m == noub)THEN
                            icval(i,j,iz) = ioub
                            tmpnopub = tmpnopub + 1
                            mcb = kk
  !                         lcval(i,j,iz)=icval(i,j,iz)
                            GO TO 171
                          ELSE
                            CYCLE loop111
                          END IF
                        END IF
                      END DO
                      DO  jj =1,nhp
                        CALL above(xpr(m),ypr(m),z,ha(m,jj), &
                                   hb(m,jj),hc(m,jj),hd(m,jj), &
                                   iflag)
                        IF (iflag == 0) THEN
                          IF(m == noub)THEN
                            icval(i,j,iz) = ioub
                            tmpnopub = tmpnopub + 1
                            mcb = kk
  !                         lcval(i,j,iz)=icval(i,j,iz)
                            GO TO 171
                          ELSE
                            CYCLE loop111
                          END IF
                        END IF
                      END DO
                      tmpnnub = tmpnnub + 1
                      mcb = kk
  !                    if(iproc==0) write(*,*)'cleared UB planes'
  ! checking the bottom planes
                      IF(ibcval(i,j,iz) == 1) THEN
                        GOTO 288
                      ENDIF
                      DO  jj = 1, nup
                        CALL above(xpr(m),ypr(m),z,a1(m,jj), &
                                   b1(m,jj),c1(m,jj),d1(m,jj), &
                                   iflag)
  
                        IF (iflag == 0) THEN
                          ibuval(i,j,iz) = 1
                          icval(i,j,iz) = indub
  !                        z = z + dsubz
  !                        CYCLE loopIZ
                        
                         GO TO 288
                        END IF
                      END DO
  
  !                    DO  jj =1,nbp
  !                      CALL above(xpr(ij,m),ypr(ij,m),z,ba1(m,jj), &
  !                                 bb1(m,jj),bc1(m,jj),bd1(m,jj),iflag)
  !                      IF(iflag == 1) THEN
  !                        icval(i,j,iz) = 11
  !                        z = z + dsubz
  !                        CYCLE loopIZ
  !                       GO TO 170
  !                      END IF
  !                    END DO
                      DO  jj =1,nop(m)
                        CALL above(xpr(m),ypr(m),z,oa1(m,jj), &
                                   ob1(m,jj),oc1(m,jj),od1(m,jj),iflag)
                        IF (iflag == 0) THEN
                          ibuval(i,j,iz) = 1
                          icval(i,j,iz) = indub
  !                        z = z + dsubz
  !                        CYCLE loopIZ
                         GO TO 288
                        END IF
                      END DO
                      DO  jj =1,nhp
                        CALL above(xpr(m),ypr(m),z,ha1(m,jj), &
                                   hb1(m,jj),hc1(m,jj), hd1(m,jj),iflag)
                        IF (iflag == 0) THEN
                          ibuval(i,j,iz) = 1
                          icval(i,j,iz) = indub
  !                        z = z + dsubz
  !                        CYCLE loopIZ
                         GO TO 288
                        END IF
                      END DO
                      tmpnspc = tmpnspc + 1
  !                    if(iproc==0) write(*,*)'cleared bottom UB pl'
    288               IF (itrough == 0) THEN
                       IF (ibcval(i,j,iz) == 1) THEN
                         icval(i,j,iz) = indcb
                         GOTO 170
                       ELSE IF (ibuval(i,j,iz) == 1) THEN
                         icval(i,j,iz) = indub
                         GOTO 170
                       ENDIF
                       icval(i,j,iz) = iots
                       GOTO 170
                      ENDIF                  
   289                mub=m
                       numUB(i,j,iz) = mub
  !                   if(iproc==0) write(*,*)'UBNUM',  numUB(i,j,iz),mub
                      mcb=numCB(i,j,iz)
                      axpr = xpr(m)
                      aypr = ypr(m)
                      azpr = z
                      icval(i,j,iz)=indub
                      GO TO 270
    171              CONTINUE
                     IF (ifill == 0) THEN
                      icval( i,j,iz) = ifind
                     
                    !  icval( i,j,iz) = 5!ifind
                      GOTO 170
                     ENDIF
                     mnx = minx(ifub)
                     mny = miny(ifub)
                     mnz = minz(ifub)
                     mxx = maxx(ifub)
                     mxy = maxy(ifub)
                     mxz = maxz(ifub)
                     mxlub = length(ifub)
                     wxub = width(ifub)
                     hxub = height(ifub)
                     GOTO 173                  
                    END DO loop111
                    z = z + dsubz
                    CYCLE loopIZ
  !                 GO TO 170
  173               CONTINUE
  !                  IF (ifill == 0) THEN
  !                   icval(i,j,iz) = ifind
  !                   GOTO 170
  !                  ENDIF
                    CALL dtran (x, xd, mnx, mxx, wxub, icerr)
                    CALL dtran (y, yd, mny, mxy, mxlub, icerr)
                    CALL dtran (z, zd, mnz, mxz, hxub, icerr)
                    mub = ifub
                    axpr = xd
                    aypr = yd
                    azpr = zd                   
    270           CONTINUE
                  IF (itrough == 0 ) THEN
                    GOTO 170
                  ENDIF
  !          
  ! Identifying Trough sets
  !       
                  nl=0
                  nneq=0
  
  !                    IF(icval(i,j,iz)==0) THEN
  !                       write(8900,*)i,j,iz,mub,mcb
  !                    ENDIF
  ! Loops over number of cosets; inner loop is number of lines inside coset
  !               if(iproc==0) write(*,*) 'entering trough'
                  loop781:  DO  jct=1,ncoset(mub)
                    loop782:  DO  knl = 1, nline(mub,jct)
                      iterr=0
                      nneq=jct
                      nl=knl
  
  ! Determine if point is within zone of trough set
                      !IF(icval(i,j,iz==0) THEN
                      ! write(8900,*)i,j,iz,mub,mcb
                      !ENDIF
                      xtran = 0.0
                      ytran = 0.0
                      ztran = 0.0
                      
                      CALL tranf(x,y,axpr,aypr,azpr,ytran,ztran, &
                                 xtran,nneq,mub,mloc,nl,mcb,iterr,notss)
  
  ! Iterr=3 outside trough set
                      !IF(icval(i,j,iz)==indub) THEN
                      ! icval(i,j,iz)=5
                      !ENDIF
  ! This step is done to make sure all unfilled space is filled with sg.
                      IF (iterr == 3) THEN
                        
                       ! IF(nl == nline(mcb,mub,nneq)) THEN
                        ! IF(nneq == ncoset(mcb,mub)) THEN
                          IF(nl == nline(mub,nneq)) THEN
                          IF(nneq == ncoset(mub)) THEN
                         
  !                    IF (iterr == 3) THEN
                        IF(iz==1) THEN
                          icval(i,j,iz)=5
                        ELSE IF(iz/=1) THEN
                             IF(mub/=tmub) THEN
                                icval(i,j,iz)=5
                             ELSE IF(icval(i,j,iz-1)==indlcf.or. &
                      icval(i,j,iz-1)==indrcf.or.icval(i,j,iz-1)==indcf) THEN
                               icval(i,j,iz)=5
                             ELSE IF(icval(i,j,iz-1)==indcs.or. &
                      icval(i,j,iz-1)==indus.or.icval(i,j,iz-1)==indts) THEN
                               icval(i,j,iz)=3
                              ELSE IF(icval(i,j,iz-1)==indcofg.or. &
                   icval(i,j,iz-1)==induofg.or.icval(i,j,iz-1)==indtofg) THEN
                               icval(i,j,iz)=4
                               ELSE IF(icval(i,j,iz-1)==indcsg3.or. &
                   icval(i,j,iz-1)==indusg3.or.icval(i,j,iz-1)==indtsg3) THEN
                               icval(i,j,iz)=5
                               ELSE IF(icval(i,j,iz-1)==indcsg2.or. &
                   icval(i,j,iz-1)==indusg2.or.icval(i,j,iz-1)==indtsg2) THEN
                               icval(i,j,iz)=21
                               ELSE IF(icval(i,j,iz-1)==indcsg1.or. &
                   icval(i,j,iz-1)==indusg1.or.icval(i,j,iz-1)==indtsg1) THEN
                               icval(i,j,iz)=20
                               ELSE
                                icval(i,j,iz)=mtind
                             ENDIF
                       ENDIF
                        tmub=mub
  
                                 !  icval(i,j,iz)=5
                            IF (ibcval(i,j,iz) == 1.OR.ibuval(i,j,iz) &
                                == 1) THEN
  
                             GOTO 169
                            ENDIF 
                             nsp=nsp+1
                            EXIT loop781
                          ELSE
                            CYCLE loop781
                          END IF
                        ELSE
  
                          CYCLE loop782
                        END IF
  
                      END IF
                     
                      icval(i,j,iz) = ind(mloc)
                      mtind=icval(i,j,iz)                    
  !                 IF(iproc==0)write(4443,*)icval(i,j,iz),ibcval(i,j,iz),&
  !                          ibuval(i,j,iz)
  !                 IF(iproc==1)write(4444,*)icval(i,j,iz),ibcval(i,j,iz),&
  !                          ibuval(i,j,iz)
  
  ! Found within the trough set; 16 planes per trough; 16th plane only considered for 
  ! open framework gravel
  
    139               DO  k = 1,15
                        CALL above(xtran,ytran,ztran,ta(k), &
                                   tb(k),tc(k), &
                                   td(k),iflag)
                        IF (iflag == 1) THEN
                        
                            IF(nl == nline(mub,nneq)) THEN
                            IF(nneq == ncoset(mub)) THEN
  
                        IF(iz==1) THEN
                          icval(i,j,iz)=5
                        ELSE IF(iz/=1) THEN
                          IF(mub/=tmub) THEN
                           icval(i,j,iz)=5
                          ELSE IF(icval(i,j,iz-1)==indlcf.or. &
                     icval(i,j,iz-1)==indrcf.or.icval(i,j,iz-1)==indcf) THEN
                               icval(i,j,iz)=5
                          ELSE IF(icval(i,j,iz-1)==indcs.or. &
                      icval(i,j,iz-1)==indus.or.icval(i,j,iz-1)==indts) THEN
                               icval(i,j,iz)=3
                         ELSE IF(icval(i,j,iz-1)==indcofg.or. &
                   icval(i,j,iz-1)==induofg.or.icval(i,j,iz-1)==indtofg) THEN
                               icval(i,j,iz)=4
                         ELSE IF(icval(i,j,iz-1)==indcsg3.or. &
                   icval(i,j,iz-1)==indusg3.or.icval(i,j,iz-1)==indtsg3) THEN
                               icval(i,j,iz)=5
                         ELSE IF(icval(i,j,iz-1)==indcsg2.or. &
                   icval(i,j,iz-1)==indusg2.or.icval(i,j,iz-1)==indtsg2) THEN
                               icval(i,j,iz)=21
                         ELSE IF(icval(i,j,iz-1)==indcsg1.or. &
                   icval(i,j,iz-1)==indusg1.or.icval(i,j,iz-1)==indtsg1) THEN
                               icval(i,j,iz)=20
  
                         ELSE
                           icval(i,j,iz)=icval(i,j,iz-1)
                         ENDIF
                        ENDIF
                        tmub=mub
  
  
                              IF (ibcval(i,j,iz) == 1.OR.ibuval(i,j,iz) &
                                      == 1) THEN
  
                               GOTO 169
                              ENDIF 
  
                              GOTO 170
                              EXIT loop781
                            ELSE
                              CYCLE loop781
                            END IF
                          ELSE
                            CYCLE loop782
                          END IF
                        END IF
                      END DO
  
                   icval(i,j,iz)=ind(mloc)  
  !                 IF(iproc==0)write(4441,*)icval(i,j,iz),ibcval(i,j,iz),&
  !                          ibuval(i,j,iz)
  !                 IF(iproc==1)write(4442,*)icval(i,j,iz),ibcval(i,j,iz),&
  !                          ibuval(i,j,iz)
                           
  
  !  Inside a trough set. If below the bottom bed plane, then assign a different indicator (ind=7)
                
                      DO  k=1,15  
                        CALL above(xtran,ytran,ztran,pta(k), &
                                   ptb(k),ptc(k), &
                                   ptd(k),iflag)
  
                        IF (iflag == 1) THEN
  ! Bottom flag indicator
                          IF (ibcval(i,j,iz)==1.OR.ibuval(i,j,iz)==1) THEN
                           GOTO 169
                          ENDIF
                          IF (icval(i,j,iz)==3) THEN
                           icval(i,j,iz)=indts
                           tmpncrfs = ncrfst + 1
                           ncrfst = tmpncrfs
                          ELSE IF (icval(i,j,iz)== 4) THEN
                           icval(i,j,iz) = indtofg
  !vf                         ncrofg = ncrofgt + 1
  !vf                         ncrofgt = ncrofg 
                           tmpncrofg = tmpncrofg+1
                          ELSE IF (icval(i,j,iz)==5) THEN
                           icval(i,j,iz) = indtsg3
                            tmpncrsg3 = tmpncrsg3+1
                          ELSE IF (icval(i,j,iz)==21) THEN
                           icval(i,j,iz) = indtsg2
                            tmpncrsg2 = tmpncrsg2+1
                          ELSE IF (icval(i,j,iz)==20) THEN
                           icval(i,j,iz) = indtsg1
                            tmpncrsg1 = tmpncrsg1+1
                          ENDIF
                          GOTO 170 
                          EXIT loop781
                        END IF
                      END DO
                   !  icval(i,j,iz)=ind(ij,lcnum(ij,mub,nneq,nl))
                   icval(i,j,iz)=ind(mloc)  
  !                    icval(i,j,iz) = ind(mcb,mub,nneq,nl)
                    !  lcval(i,j,iz)=icval(i,j,iz)
  !                    if( icval(i,j,iz)==0) then
  !                write(711,*) icval(i,j,iz),ind(ij,mloc)
  !                write(711,*) lcnum(ij,mub,nneq,nl),mloc,ij,mub,nneq,nl
  !                    endif
                      IF (ibcval(i,j,iz)==1.OR.ibuval(i,j,iz)==1) THEN
                           GOTO 169
                      ENDIF
                      GOTO 170
                      EXIT loop781
                    END DO loop782
                  END DO loop781
  169            CONTINUE
                      mtind=icval(i,j,iz)                    
             IF(ibcval(i,j,iz)==1) THEN
              IF(icval(i,j,iz)==3) THEN
               icval(i,j,iz)=indcs
  !             nfsc=nfsct+1
  !             nfsct=nfsc
               tmpncrfs=ncrfst+1
               ncrfst=tmpncrfs
              ELSE IF (icval(i,j,iz)==4) THEN
               icval(i,j,iz)=indcofg
  !             nofgc=nofgct+1
  !             nofgct=nofgc
  !vf              ncrofg=ncrofgt+1
  !vf              ncrofgt=ncrofg
                 tmpncrofg = tmpncrofg+1
              ELSE IF(icval(i,j,iz)==5) THEN
               icval(i,j,iz)=indcsg3
                tmpncrsg3=ncrsgt3+1
                ncrsgt3=tmpncrsg3
              ELSE IF(icval(i,j,iz)==21) THEN
               icval(i,j,iz)=indcsg2
                tmpncrsg2=ncrsgt2+1
                ncrsgt2=tmpncrsg2
              ELSE IF(icval(i,j,iz)==20) THEN
               icval(i,j,iz)=indcsg1
                tmpncrsg1=ncrsgt1+1
                ncrsgt1=tmpncrsg1
              ENDIF
             ELSE IF (ibuval(i,j,iz)==1) THEN
              IF(icval(i,j,iz)==3) THEN
               icval(i,j,iz)=indus
  !             nfsc=nfsct+1
  !             nfsct=nfsc
                tmpncrfs=ncrfst+1
                ncrfst=tmpncrfs
              ELSE IF(icval(i,j,iz)==4) THEN
               icval(i,j,iz)=induofg
  !             nofgc=nofgct+1
  !             nofgct=nofgc
  !vf              ncrofg=ncrofgt+1
  !vf              ncrofgt=ncrofg
               tmpncrofg = tmpncrofg+1
              ELSE IF(icval(i,j,iz)==5) THEN
               icval(i,j,iz)=indusg3
                tmpncrsg3=ncrsgt3+1
                ncrsgt3=tmpncrsg3
              ELSE IF(icval(i,j,iz)==21) THEN
               icval(i,j,iz)=indusg2
                tmpncrsg2=ncrsgt2+1
                ncrsgt2=tmpncrsg2
              ELSE IF(icval(i,j,iz)==20) THEN
               icval(i,j,iz)=indusg1
                tmpncrsg1=ncrsgt1+1
                ncrsgt1=tmpncrsg1
              ENDIF
             ENDIF
  170            CONTINUE
                
  !                 IF(iproc==0)write(4451,*)icval(i,j,iz),i,j,iz
  !                 IF(iproc==1)write(4452,*)icval(i,j,iz),i,j,iz
  !      IF(icval(i,j,iz)==0) THEN
  !        write(8001,*)i,j,iz,xtran,ytran,ztran
  !      ENDIF
  !               ibcval = 0
  !               ibuval = 0
                      mtind=icval(i,j,iz)                    
  !                IF (ibin == 1) THEN
  !                  WRITE(7) x,y,ztr,ival
  !                ELSE
  !!                 WRITE(7,*) x,y,ztr,icval (i,j,iz)
  !                  WRITE(7,'(6(I6))') i,j,iz,icval(i,j,iz),numCB(i,j,iz), &
  !                             numUB(i,j,iz)
  !                ENDIF
                   
  ! Calculate volume proportion of each texture.
  !  3 = fine sand, 4 = open framework gravel, 5= sandy gravel
            
  !                nctot=nctemp+1
  !                nctemp=nctot
                  IF(icval (i,j,iz) == 3) THEN
                    ! nfsc=nfsct+1
                    tmpnfsc=tmpnfsc+1
                    nfsct=tmpnfsc
  
                  ELSE IF (icval (i,j,iz) == 4) THEN
  !vf                  nofgc=nofgt+1
  !vf                  nofgt=nofgc
                    ! nofgc = nofgc+1
                    tmpnofgc=tmpnofgc+1
                  ELSE IF (icval (i,j,iz) == 5) THEN
                    ! nsgc3=nsgct3+1
                    tmpnsgc3=nsgct3+1
                    nsgct3=tmpnsgc3
                  ELSE IF (icval (i,j,iz) == 21) THEN
                    tmpnsgc2=nsgct2+1
                    nsgct2=tmpnsgc2
                  ELSE IF (icval (i,j,iz) == 20) THEN
                    tmpnsgc1=nsgct1+1
                    nsgct1=tmpnsgc1
                  END IF
                  IF(icval (i,j,iz) == indts) THEN
                   ! IF(lcval == 3) THEN
                      tmpncrfs=ncrfst+1
                      ncrfst=tmpncrfs
                    ELSE IF (icval (i,j,iz) == indtofg) THEN
  !vf                    ncrofg=ncrofgt+1
  !vf                    ncrofgt=ncrofg
                      tmpncrofg = tmpncrofg+1
                    ELSE IF(icval (i,j,iz) == indtsg3) THEN
                      tmpncrsg3=ncrsgt3+1
                      ncrsgt3=tmpncrsg3
                    ELSE IF(icval (i,j,iz) == indtsg2) THEN
                      tmpncrsg2=ncrsgt2+1
                      ncrsgt2=tmpncrsg2
                    ELSE IF(icval (i,j,iz) == indtsg1) THEN
                      tmpncrsg1=ncrsgt1+1
                      ncrsgt1=tmpncrsg1
                    !END IF
                  END IF
                    z = z + dsubz
                  END DO loopIZ
                  y = y + dsuby
                END DO
                x = x + dsubx
              END DO
            END DO
          END DO
          
          ! lcval(nodeLeft:nodeRight,nodeFront:nodeBack, &
          !  nodeUp:nodeDown)=icval(nodeLeft:nodeRight,nodeFront:nodeBack, &
          !        nodeUp:nodeDown)
          ! lcval(nodeLeft:nodeRight,&
          ! nodeFront1(localIndices):nodeBack1(localIndices), &
          ! nodeUp:nodeDown)=icval(nodeLeft:nodeRight,&
          !   nodeFront1(localIndices):nodeBack1(localIndices), &
          !   nodeUp:nodeDown)
          lcval(nodeLeft:nodeRight,&
          tmpNodeFront1:tmpNodeBack1, &
          nodeUp:nodeDown)=icval(nodeLeft:nodeRight,&
          tmpNodeFront1:tmpNodeBack1, &
            nodeUp:nodeDown)
           
  !         write(4889,*)iproc
  !         write(4889,*)lcval
  !        write(5782,*)'compound bar, iproc', kk,iproc
  ! Deallocate arrays
  
          DEALLOCATE(nop, STAT=IERR)
          adum = 'nop-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(a, STAT=IERR)
          adum = 'a-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(b, STAT=IERR)
          adum = 'b-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(c, STAT=IERR)
          adum = 'c-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(d, STAT=IERR)
          adum = 'd-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(ba, STAT=IERR)
          adum = 'ba-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(bb, STAT=IERR)
          adum = 'bb-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(bc, STAT=IERR)
          adum = 'bc-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(bd, STAT=IERR)
          adum = 'bd-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(oa, STAT=IERR)
          adum = 'oa-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(ob, STAT=IERR)
          adum = 'ob-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(oc, STAT=IERR)
          adum = 'oc-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(od, STAT=IERR)
          adum = 'od-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(ha, STAT=IERR)
  
          adum = 'ha-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(hb, STAT=IERR)
          adum = 'hb-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(hc, STAT=IERR)
          adum = 'hc-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(hd, STAT=IERR)
          adum = 'hd-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(a1, STAT=IERR)
          adum = 'a1-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(b1, STAT=IERR)
          adum = 'b1-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(c1, STAT=IERR)
          adum = 'c1-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(d1, STAT=IERR)
          adum = 'd1-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(ba1, STAT=IERR)
          adum = 'ba1-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(bb1, STAT=IERR)
          adum = 'bb1-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(bc1, STAT=IERR)
          adum = 'bc1-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(bd1, STAT=IERR)
          adum = 'bd1-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(oa1, STAT=IERR)
          adum = 'oa1-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(ob1, STAT=IERR)
          adum = 'ob1-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(oc1, STAT=IERR)
          adum = 'oc1-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(od1, STAT=IERR)
          adum = 'od1-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(ha1, STAT=IERR)
          adum = 'ha1-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(hb1, STAT=IERR)
          adum = 'hb1-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(hc1, STAT=IERR)
          adum = 'hc1-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(hd1, STAT=IERR)
          adum = 'hd1-deal'
          CALL CHECK(adum,ierr)
  
  ! Allocate for coset.out read
  
          DEALLOCATE(ncoset, STAT=IERR)
          adum = 'ncoset-deal'
          CALL CHECK(adum,ierr)
  
  ! Allocate for line.out read
  
          DEALLOCATE(nline, STAT=IERR)
          adum = 'nline-deal'
          CALL CHECK(adum,ierr)
  
  ! Allocate for location.out read  - this can be deleted??
  
          DEALLOCATE(lcnum, STAT=IERR)
          adum = 'lcnum-deal'
          CALL CHECK(adum,ierr)
  
  ! Deallocate for trgena.out read
  
          DEALLOCATE(xc, STAT=IERR)
          adum = 'xc-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(yc, STAT=IERR)
          adum = 'yc-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(xinfl, STAT=IERR)
          adum = 'xinfl-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(yinfl, STAT=IERR)
          adum = 'yinfl-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(r, STAT=IERR)
          adum = 'r-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(alpha, STAT=IERR)
          adum = 'alpha-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(larc, STAT=IERR)
          adum = 'larc-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(lsum, STAT=IERR)
          adum = 'lsum-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(iinf, STAT=IERR)
          adum = 'iinf-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(xx, STAT=IERR)
          adum = 'xx-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(width, STAT=IERR)
          adum = 'width-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(length, STAT=IERR)
          adum = 'length-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(height, STAT=IERR)
          adum = 'height-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(minx, STAT=IERR)
          adum = 'minx-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(miny, STAT=IERR)
          adum = 'miny-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(minz, STAT=IERR)
          adum = 'minz-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(maxx, STAT=IERR)
          adum = 'maxx-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(maxy, STAT=IERR)
          adum = 'maxy-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(maxz, STAT=IERR)
          adum = 'maxz-deal'
          CALL CHECK(adum,ierr)
  
  ! Deallocate for indicator.out read
  
          DEALLOCATE(ind, STAT=IERR)
          adum = 'ind-deal'
          CALL CHECK(adum,ierr)
  
  ! Allocate for out1.out (trough plane coefficients) read
  
          DEALLOCATE(ta, STAT=IERR)
          adum = 'ta-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(tb, STAT=IERR)
          adum = 'tb-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(tc, STAT=IERR)
          adum = 'tc-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(td, STAT=IERR)
          adum = 'td-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(pta, STAT=IERR)
          adum = 'pta-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(ptb, STAT=IERR)
          adum = 'ptb-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(ptc, STAT=IERR)
          adum = 'ptc-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(ptd, STAT=IERR)
          adum = 'ptd-deal'
          CALL CHECK(adum,ierr)
  
  ! Allocate for lenstat.out read
  
          DEALLOCATE(slope, STAT=IERR)
          adum = 'slope-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(yint, STAT=IERR)
          adum = 'yint-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(xmin, STAT=IERR)
          adum = 'xmin-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(xmax, STAT=IERR)
          adum = 'xmax-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(cymin, STAT=IERR)
          adum = 'cymin-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(czmax, STAT=IERR)
          adum = 'czmax-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(twidth, STAT=IERR)
          adum = 'twidth-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(ppdn, STAT=IERR)
          adum = 'ppdn-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(pslope, STAT=IERR)
          adum = 'pslope-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(pyint, STAT=IERR)
          adum = 'pyint-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(llen, STAT=IERR)
          adum = 'llen-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(delh, STAT=IERR)
          adum = 'delh-deal'
          CALL CHECK(adum,ierr)
  
  ! Deallocate for x-y-z loop
  
          DEALLOCATE(xpr, STAT=IERR)
          adum = 'xpr-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(ypr, STAT=IERR)
          adum = 'ypr-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(ierru, STAT=IERR)
          adum = 'ierru-deal'
          CALL CHECK(adum,ierr)
  
        END DO  !End compound bar loop read
        CLOSE(33)
  !      icval(nodeLeft:nodeBack,nodeFront:nodeBack,nodeUp:nodeDown)= &
  !          icval(nodeLeft:nodeBack,nodeFront:nodeBack,nodeDown:nodeUp:-1)
        IF (ibin == 1) THEN
  !        filename(1:11) = 'testin.out.'
  !        WRITE(filename(12:),intform) iproc
  !        OPEN(6667,FILE=filename,FORM='unformatted', &
  !             STATUS='unknown')
  !      WRITE(*,*) 'nodeFront, nodeBack,nodeLeft, nodeRight',nodeFront, nodeBack,nodeLeft, nodeRight
      !    DO j = nodeFront, nodeBack
          ! DO j = nodeFront1(localIndices), nodeBack1(localIndices)
          DO j = tmpNodeFront1,tmpNodeBack1
           DO i = nodeLeft, nodeRight
             DO  iz = nodeUp, nzpp*nodeDown
                WRITE(7)icval(i,j,iz), i,j,iz
                IF(icval(i,j,iz)==0.OR.icval(i,j,iz)==1.OR.&
                       icval(i,j,iz)==2.OR.icval(i,j,iz)==15) THEN
                 WRITE(6667,*)icval(i,j,iz),i,j,iz
  !               IF(iproc==1)WRITE(6668,*)icval(i,j,iz),i,j,iz
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ELSE
          ! DO j = nodeFront, nodeBack
          ! DO j = nodeFront1(localIndices), nodeBack1(localIndices)
          DO j = tmpNodeFront1,tmpNodeBack1
            DO i = nodeLeft, nodeRight
              DO  iz = nodeUp, nzpp*nodeDown
  !vf              IF(iproc==0) THEN
  !vf                 WRITE(796,*) i,j,iz,icval(i,j,iz)
  !vf              ELSE IF(iproc==1) THEN
  !vf                 WRITE(797,*) i,j,iz,icval(i,j,iz)
  !vf              ELSE IF(iproc==2) THEN
  !vf                 WRITE(798,*) i,j,iz,icval(i,j,iz)
  !vf              ENDIF
                 WRITE(7,'(6(I6))') i,j,iz,icval(i,j,iz),numUB(i,j,iz), &
                               numCB(i,j,iz)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        CLOSE(7)
          
          DEALLOCATE(icval, STAT=IERR)
          adum = 'icval-deal'
          CALL CHECK(adum,ierr)
          
          DEALLOCATE(lcval, STAT=IERR)
          adum = 'icval-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(ibcval, STAT=IERR)
          adum = 'ibcval-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(ibuval, STAT=IERR)
          adum = 'ibuval-deal'
          CALL CHECK(adum,ierr)
  
          !Deallocate arrays used for each component bar
          DEALLOCATE(cba, STAT=IERR)
          adum = 'cba-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(cbb, STAT=IERR)
          adum = 'cbb-deal'
          CALL CHECK(adum,ierr)        
  
          DEALLOCATE(cbc, STAT=IERR)
          adum = 'cbc-deal'
          CALL CHECK(adum,ierr) 
  
          DEALLOCATE(cbd, STAT=IERR)
          adum = 'cbd-deal'
          CALL CHECK(adum,ierr) 
  
          DEALLOCATE(xcc, STAT=IERR)
          adum = 'xcc-deal'
          CALL CHECK(adum,ierr) 
  
          DEALLOCATE(ycc, STAT=IERR)
          adum = 'ycc-deal'
          CALL CHECK(adum,ierr) 
  
          DEALLOCATE(xinflc, STAT=IERR)
          adum = 'xinflc-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(yinflc, STAT=IERR)
          adum = 'yinflc-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(rc, STAT=IERR)
          adum = 'rc-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(alphacbc, STAT=IERR)
          adum = 'alphacbc-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(lsumc, STAT=IERR)
          adum = 'lsumc-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(numCB, STAT=IERR)
          adum = 'numCB-deal'
          CALL CHECK(adum,ierr)
  
          DEALLOCATE(numUB, STAT=IERR)
          adum = 'numUB-deal'
          CALL CHECK(adum,ierr)
  
          tmpNodeFront1 = tmpNodeBack1+1
          ! IF (tmpnfsc>nfsc) THEN
          nfsc=nfsc+tmpnfsc
          ! ENDIF
          nofgc=nofgc+tmpnofgc
          nsgc3=nsgc3+tmpnsgc3
          nsgc2=nsgc2+tmpnsgc2 
          nsgc1=nsgc1+tmpnsgc1
          nctot=nctot+tmpnctot
          ncrofg=ncrofg+tmpncrofg
          ncrfs=ncrfs+tmpncrfs
          ncrsg3=ncrsg3+tmpncrsg3
          ncrsg2=ncrsg2+tmpncrsg2
          ncrsg1=ncrsg1+tmpncrsg1
          nspc=nspc+tmpnspc
          nnub=nnub+tmpnnub
          nosf=nosf+tmpnosf 
          nncb=nncb+tmpnncb
          nbbcb=nbbcb+tmpnbbcb
          nopub=nopub+tmpnopub
          nnacb=nnacb+tmpnnacb
        ENDDO
  !  Get global sums
        WRITE(*,*) 'The number of cells is ', nctot
        WRITE(*,*) 'nspc number is ', nspc
        CALL MPI_AllReduce(nctot,nctotG,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
        adum = 'nctotG'
        CALL CHECK(adum,ierr)
  
        CALL MPI_AllReduce(nspc,nspcG,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
        adum = 'nspc'
        CALL CHECK(adum,ierr)
  
        !CALL MPI_AllReduce(ncitx,ncitxG,1,MPI_INTEGER,MPI_SUM,comm3d,ierr) !Naum never used
        !adum = 'ncitx'
        !CALL CHECK(adum,ierr)
  
        CALL MPI_AllReduce(nfsc,nfscG,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
        adum = 'nfsc'
        CALL CHECK(adum,ierr)
  
        CALL MPI_AllReduce(nofgc,nofgcG,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
        adum = 'nofgc'
        CALL CHECK(adum,ierr)
  
        CALL MPI_AllReduce(nsgc3,nsgc3G,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
        adum = 'nsgc3'
        CALL CHECK(adum,ierr)
  
        CALL MPI_AllReduce(nsgc2,nsgc2G,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
        adum = 'nsgc2'
        CALL CHECK(adum,ierr)
  
        CALL MPI_AllReduce(nsgc1,nsgc1G,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
        adum = 'nsgc1'
        CALL CHECK(adum,ierr)
  
        CALL MPI_AllReduce(ncrfs,ncrfsG,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
        adum = 'ncrfs'
        CALL CHECK(adum,ierr)
  
        CALL MPI_AllReduce(ncrofg,ncrofgG,1,MPI_INTEGER,MPI_SUM,comm3d, &
                           ierr)
        adum = 'ncrofg'
        CALL CHECK(adum,ierr)
  
        CALL MPI_AllReduce(ncrsg3,ncrsg3G,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
        adum = 'ncrsg3'
        CALL CHECK(adum,ierr)
  
        CALL MPI_AllReduce(ncrsg2,ncrsg2G,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
        adum = 'ncrsg2'
        CALL CHECK(adum,ierr)
  
        CALL MPI_AllReduce(ncrsg1,ncrsg1G,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
        adum = 'ncrsg1'
        CALL CHECK(adum,ierr)
  
        CALL MPI_AllReduce(nosf,nosfG,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
        adum = 'nosf'
        CALL CHECK(adum,ierr)
        
        CALL MPI_AllReduce(nnacb,nnacbG,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
        adum = 'nnacb'
        CALL CHECK(adum,ierr)
        
        CALL MPI_AllReduce(nbbcb,nbbcbG,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
        adum = 'nbbcb'
        CALL CHECK(adum,ierr)
  
        CALL MPI_AllReduce(nncb,nncbG,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
        adum = 'nncb'
        CALL CHECK(adum,ierr)
       
        CALL MPI_AllReduce(nopub,nopubG,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
        adum = 'nopub'
        CALL CHECK(adum,ierr)
  
        CALL MPI_AllReduce(nnub,nnubG,1,MPI_INTEGER,MPI_SUM,comm3d,ierr)
        adum = 'nnub'
        CALL CHECK(adum,ierr)
  
  !  Writeout results of column proportion to results.out
  
        IF (iproc == 0) THEN
          WRITE(16,59)'total # cells', nctotG !OK
          WRITE(16,59)'total # os cb cells',nosfG !OK
          WRITE(16,59)'to textured cbar cells',nnacbG !OK
          WRITE(16,59)'totalcb boundary cells',nbbcbG !OK
          WRITE(16,59)'exc cf cb cells',nncbG !OK
          WRITE(16,59)'outside ub cells',nopubG !OK
          WRITE(16,59)'inside ub cells',nnubG !OK
          WRITE(16,59)'total # fs cells', nfscG !OK
          WRITE(16,59)'total # ofg cells', nofgcG !OK
          WRITE(16,59)'total # pth sg cells', nsgc3G !OK
          WRITE(16,59)'total # m sg cells', nsgc2G !OK
          WRITE(16,59)'total # bm sg cells', nsgc1G !OK
          WRITE(16,59)'excl bbed # fs cells', ncrfsG !OK
          WRITE(16,59)'excl bbed # ofg cells', ncrofgG !OK
          WRITE(16,59)'excl bbed # pth sg cells', ncrsg3G !OK
          WRITE(16,59)'excl bbed # m sg cells', ncrsg2G !OK
          WRITE(16,59)'excl bbed # bm sg cells', ncrsg1G !OK
          
          WRITE(16,59)'nspc', nspcG ! Naum What is it?
          !WRITE(16,59)'unfilled space', nspeG
          !write(16,59)'outside fill ub',nopuboG
          !write(16,59)'outside fill ub',nocubeG
          close(16) !Naum
          WRITE(*,59)'total # cells', nctotG
          WRITE(*,59)'total # os cb cells',nosfG
          WRITE(*,59)'to textured cbar cells',nnacbG
          WRITE(*,59)'totalcb boundary cells',nbbcbG
          WRITE(*,59)'exc cf cb cells',nncbG
          WRITE(*,59)'outside ub cells',nopubG
          WRITE(*,59)'inside ub cells',nnubG
          WRITE(*,59)'total # fs cells', nfscG
          WRITE(*,59)'total # ofg cells', nofgcG
          WRITE(*,59)'total # pth sg cells', nsgc3G
          WRITE(*,59)'total # bm sg cells', nsgc1G
          WRITE(*,59)'total # m sg cells', nsgc2G
          WRITE(*,59)'excl bbed # fs cells', ncrfsG
          WRITE(*,59)'excl bbed # ofg cells', ncrofgG
          WRITE(*,59)'excl bbed # pth sg cells', ncrsg3G
          WRITE(*,59)'excl bbed # bm sg cells', ncrsg2G
          WRITE(*,59)'excl bbed # m sg cells', ncrsg1G
          WRITE(*,59)'nspc', nspcG ! Naum What is it?
          !    WRITE(*,59)'unfilled space', nspeG
          !write(*,59)'outside fill ub',nopuboG
          !write(*,59)'outside fill ub',nocubeG
        ENDIF
  !      CALL crfgamt(nx,ny,nz)
  
        999 CONTINUE
  !     WRITE(*,*)mm
  !      STOP
        END SUBROUTINE mergeAll
  !***********************************************************************
  
        SUBROUTINE READ_GRID(ans,xulc,yulc,zulc,xdim,ydim,zdim,xspdim, &
                             yspdim,zspdim,nx,ny,nz)
  
        INTEGER*4                                   :: nx
        INTEGER*4                                   :: ny
        INTEGER*4                                   :: nz
        INTEGER*4                                   :: ans
  
        REAL*4                                      :: xulc
        REAL*4                                      :: yulc
        REAL*4                                      :: zulc
        REAL*4                                      :: xdim
        REAL*4                                      :: ydim
        REAL*4                                      :: zdim
        REAL*4                                      :: xspdim
        REAL*4                                      :: yspdim
        REAL*4                                      :: zspdim
  
  
  ! Screen read of control data: Coordinate limits
  ! of universe space and pixel size
  
  !    WRITE(*,101)
  !    WRITE(*,1) 'Output data can be written in two formats: '
  !    WRITE(*,1) '  1) Integer (takes less space and time)'
  !    WRITE(*,1) '  2) Real (allows averaging)'
  !    WRITE(*,101)
  !    WRITE(*,1) 'Which would you like (1 or 2)? '
  !    READ(*,*) ans
  !    WRITE(*,1) 'X-coord. of upper left corner of space? '
  !    READ(*,*) xulc
  !    WRITE(*,1) 'Y-coord. of upper left corner of space? '
  !    READ(*,*) yulc
  !   WRITE(*,1) 'Z-coord. of upper left corner of space? '
  !   READ(*,*) zulc
  !   WRITE(*,1) 'X-dimension of pixel? '
  !   READ(*,*) xdim
  !   WRITE(*,1) 'Y-dimension of pixel? '
  !   READ(*,*) ydim
  !   WRITE(*,1) 'Z-dimension of pixel? '
  !   READ(*,*) zdim
  !   WRITE(*,1) 'X-dimension of space? '
  !   READ(*,*) xspdim
  !   WRITE(*,1) 'Y-dimension of space? '
  !   READ(*,*) yspdim
  !   WRITE(*,1) 'Z-dimension of space? '
  !   READ(*,*) zspdim
  
  !     ans = 1
  !     xulc = 50.0
  !     yulc = 20.0
  !     zulc = 0.0
  !     xdim = 1.0
  !     ydim = 1.0
  !     zdim = 0.1
  !     xspdim = 300.0
  !     yspdim = 600.0
  !     zspdim = 3.0
  
          OPEN(21,FILE='sample.dat',STATUS='old') ! Naum 06.07.13
          READ(21,*)xulc     
          READ(21,*)yulc    
          READ(21,*)zulc    
          READ(21,*)xdim    
          READ(21,*)ydim    
          READ(21,*)zdim    
          READ(21,*)xspdim  
          READ(21,*)yspdim  
          READ(21,*)zspdim  ! Naum 06.07.13
          
       ans = 1
       
        nx = INT((xspdim+xdim*1.e-02)/xdim)
        ny = INT((yspdim+ydim*1.e-02)/ydim)
        nz = INT((zspdim+zdim*1.e-02)/zdim)
  
        WRITE(*,101)
        WRITE(*,*) 'nx = ',nx
        WRITE(*,*) 'ny = ',ny
        WRITE(*,*) 'nz = ',nz
        WRITE(*,101)
  
      1 FORMAT(1X,A)
    101 FORMAT(1X,/)
  
        RETURN
  
        END SUBROUTINE READ_GRID
  
  !***********************************************************************
  !     Subroutine ABOVE checks the location of a point relative to a
  !        plane and returns the appropriate flag:
  
  !        IFLAG = 0  ....   Point is "above" plane
  !        IFLAG = 1  ....   Point is "below" plane
  
  !     Where "above" and "below" are defined in terms of an inverted
  !        z axis (i.e., if z(point) > z(plane), the point is below the
  !        plane)
  
  !***********************************************************************
  
        SUBROUTINE above(x,y,z,a,b,c,d,iflag)
       
        REAL*4, INTENT(IN)                         :: x
        REAL*4, INTENT(IN)                         :: y
        REAL*4, INTENT(IN OUT)                     :: z
        REAL*4, INTENT(IN OUT)                     :: a
        REAL*4, INTENT(IN)                         :: b
        REAL*4, INTENT(IN)                         :: c
        REAL*4, INTENT(IN OUT)                     :: d
        REAL*4                                     :: zplane
        INTEGER*4, INTENT(OUT)                     :: iflag
        
        zplane = -(a*x + b*y + d)/c
        
        IF (z >= zplane) THEN
          iflag = 1
        ELSE
          iflag = 0
        END IF
        
        RETURN
        END SUBROUTINE above
  
  !***********************************************************************
  
  !     Subroutine COTRAN transforms an (x,y) coordinate point in the
  !        real space to a point (x',y') in the straight-centerline
  !        geometric space (without symmetry).
  
  !***********************************************************************
  
        SUBROUTINE cotran(x,y,xpr,ypr,narc,xc,yc,xinfl,yinfl, r,alpha, &
                          lsum,width,ierr)
  
        REAL*4, INTENT(IN)                         :: x
        REAL*4, INTENT(IN)                         :: y
        REAL*4, INTENT(OUT)                        :: xpr
        REAL*4, INTENT(OUT)                        :: ypr
        REAL*4, INTENT(IN)                         :: xc(3)
        REAL*4, INTENT(IN)                         :: yc(3)
        REAL*4, INTENT(IN)                         :: xinfl(3)
        REAL*4, INTENT(IN)                         :: yinfl(3)
        REAL*4, INTENT(IN)                         :: r(3)
        REAL*4, INTENT(IN OUT)                     :: alpha(3)
        REAL*4, INTENT(IN)                         :: lsum(3)
        REAL*4, INTENT(IN)                         :: width
        REAL*4                                     :: rtest
        REAL*4                                     :: arg
        REAL*4                                     :: a
        REAL*4                                     :: b
        REAL*4                                     :: c
        REAL*4                                     :: cang
        REAL*4                                     :: blow
        REAL*4                                     :: ahigh
        REAL*4                                     :: bhigh
        REAL*4                                     :: rplus
        REAL*4                                     :: rminus
        REAL*4                                     :: dx
        INTEGER*4, INTENT(OUT)                     :: ierr
        INTEGER*4, INTENT(IN)                      :: narc
        INTEGER*4                                  :: jj
        INTEGER*4                                  :: jerr1
        INTEGER*4                                  :: jerr2
  
      1 FORMAT(1X,a)
    101 FORMAT(1X,/)
  
  ! Find the arc associated with the specified point
  
        DO  i = 1, narc
          rtest = xydist(x,y,xc(i),yc(i))
          rplus = r(i)+width/2
          rminus = r(i)-width/2
          IF (rtest <= rplus) THEN
            IF (rtest >= rminus) THEN
              CALL line(xc(i),yc(i),xinfl(i),yinfl(i),alow,blow,jerr1)
              CALL line(xc(i),yc(i),xinfl(i+1),yinfl(i+1),ahigh, &
                        bhigh, jerr2)
              IF (jerr1 == 1 .or. jerr2 == 1) THEN
                WRITE(*,101)
                WRITE(*,1) '*** Error code from LINE = 1 *** '
                WRITE(*,101)
                ierr = 1
              END IF
              ytest = alow*x + blow
              IF (y > ytest) THEN
                ytest = ahigh*x + bhigh
                IF (y < ytest) THEN
                  jj = i
                  GO TO 110
                END IF
              END IF
            END IF
          END IF
        END DO
        ierr = 2
        GO TO 999
    110 CONTINUE
  
  !     Calculate xpr, ypr
  
        a = xydist(x,y,xinfl(jj),yinfl(jj))
        b = xydist(x,y,xc(jj),yc(jj))
        c = xydist(xinfl(jj),yinfl(jj),xc(jj),yc(jj))
        arg = (b**2.0 + c**2.0 - a**2.0)/(2.0*b*c)
  !      if(iproc==0) write(*,*)a,b,c,arg
  !      if(iproc==0) write(*,*)x,y
        cang = ACOS(arg)
        arcl = cang*r(jj)
        ypr = lsum(jj) + arcl
        dx = r(jj)-b
        IF (xc(jj) < x) THEN
          xpr = width/2.0 - dx
        ELSE
          xpr = width/2.0 + dx
        END IF
  
    999 CONTINUE
        RETURN
        END SUBROUTINE cotran
  
  
  !***********************************************************************
  
  !     Subroutine LINE calculates the slope and y-intercept parameters
  !     of a line passing through the two points given (in 2D space).
  !     This program also checks if the line is vertical, in which
  !     case the index is set to 1.
  
  !***********************************************************************
  
        SUBROUTINE line(x1,y1,x2,y2,slope,yint,indx)
       
       
        REAL*4, INTENT(IN)                         :: x1
        REAL*4, INTENT(IN)                         :: y1
        REAL*4, INTENT(IN OUT)                     :: x2
        REAL*4, INTENT(IN OUT)                     :: y2
        REAL*4, INTENT(OUT)                        :: slope
        REAL*4, INTENT(OUT)                        :: yint
        INTEGER*4, INTENT(OUT)                     :: indx
        
        indx = 0
        diff = ABS(x2-x1)
        IF (diff < 0.00001) THEN
          indx = 1
        ELSE
          slope = (y2-y1)/(x2-x1)
          yint = y1 - slope*x1
        END IF
        
        RETURN
        END SUBROUTINE line
  !***********************************************************************
  
  !     Function XYDIST calculates the distance in cartesian coordinates
  !     between two points  (x1,y1) and (x2,y2).
  
  !***********************************************************************
  
        FUNCTION xydist(x1,y1,x2,y2)
       
        REAL*4, INTENT(IN OUT)                     :: x1
        REAL*4, INTENT(IN OUT)                     :: y1
        REAL*4, INTENT(IN OUT)                     :: x2
        REAL*4, INTENT(IN OUT)                     :: y2
        REAL*4                                     :: xydist
        REAL*4                                     :: arg
        
        arg = (ABS(x1-x2))**2.0 + (ABS(y1-y2))**2.0
        xydist = SQRT(arg)
        
        RETURN
        END FUNCTION xydist
  !***********************************************************************
  
  !     Subroutine TRANF transforms an (x,y,z) coordinate point in the
  !     real space to a point (xtran,ytran,ztran) in the straight-centerline
  !     geometric space. This subroutine identifies the location of a 
  !     trough sets for each point. First, perpendicular distance of the 
  !     point from the centre line and line parallel to it are determined.
  !     Next, checked to see if it meets some criteria specified (here it 
  !     is the thickness of trough set). If yes, then calculate xtran,ztran
  !     then calculate the point of intersection (ytran) of the centre line
  !     and a line perpendicular to the centre line that passes through (x,y,z).
  !     Return xtran, ytan, ztran to the main program along with the info 
  !     about the location of trough set. If outisde all trough sets, then 
  !     iterr=3. Input file: requires an include file tpix.inc
  !***********************************************************************
  
        SUBROUTINE tranf(x,y,axpr,aypr,azpr,ytran,ztran,xtran,nneq,  &
                         mub,mloc,nl,mcb,iterr,notss)
        USE tpix
  
        REAL*4, INTENT(IN OUT)                     :: x
        REAL*4, INTENT(IN OUT)                     :: y
        REAL*4, INTENT(IN OUT)                     :: axpr
        REAL*4, INTENT(IN OUT)                         :: aypr
        REAL*4, INTENT(IN OUT)                         :: azpr
        REAL*4, INTENT(OUT)                        :: ytran
        REAL*4, INTENT(OUT)                        :: ztran
        REAL*4, INTENT(OUT)                        :: xtran
        REAL*4                                     :: nr
        REAL*4                                     :: dr
        REAL*4                                     :: mdist
        REAL*4                                     :: xpr
        REAL*4                                     :: ypr
        REAL*4                                     :: rdist1
        REAL*4                                     :: rdist2
        REAL*4                                     :: rdist
        REAL*4                                     :: pdtest
        REAL*4                                     :: aslope
        REAL*4                                     :: ayint
        REAL*4                                     :: yintr
        REAL*4                                     :: zintr
        REAL*4                                     :: llength
        REAL*4                                     :: xc
        INTEGER*4, INTENT(IN OUT)                  :: nneq
        INTEGER*4, INTENT(IN)                      :: mub
        INTEGER*4, INTENT(OUT)                     :: mloc
        INTEGER*4, INTENT(IN OUT)                  :: nl
        INTEGER*4, INTENT(IN)                      :: mcb
        INTEGER*4, INTENT(OUT)                     :: iterr
        INTEGER*4                                  :: jj
        INTEGER*4                                  :: i
        INTEGER*4                                  :: jerr1
        INTEGER*4                                  :: jerr2
        INTEGER*4                                  :: kj
        INTEGER*4                                  :: kn
        INTEGER*4                                  :: jn
        INTEGER*4                                  :: ik
        INTEGER*4                                  :: nind
        INTEGER*4                                  :: nq1
        INTEGER*4                                  :: j
        INTEGER*4                                  :: k
  
      1 FORMAT(1X,a)
     97 FORMAT(31X,i2)
     98 FORMAT(1X,f10.4,2X,f10.4)
    101 FORMAT(1X,/)
    141 FORMAT(4(1X,f7.2), 3(1X, i6))
   
  
        i=mub
        !ik=mcb Naum
        k=nneq
        j=nl
        
   550  nq1=lcnum(i,k,j)
  !     write(*,*) '3539 : ',i,k,j, nq1,notss
  !  550     nq1=lcnum(1,mub,nneq,nl)
  !   Calculate perpendicular distance of the point from the line and its parallel line (rdist, rdist2).
  !   If combined perpendicular distance (rdist+rdist2) does not meet the specification,
  !   then the point is not in this trough region.
  !   Proceed to next region and repeat the same steps.
  
        rdist=ABS((aypr-azpr*slope(nq1)-yint(nq1))/ &
              (SQRT(1+(slope(nq1))**2)))
        a=(SQRT(1+(pslope(nq1))**2))
        b=(aypr-azpr*pslope(nq1)-pyint(nq1))
        rdist2=ABS(b/a)
        pdtest = rdist2 + rdist
        diff = ABS(pdtest - ppdn(nq1))
  
        IF (diff > 0.1) THEN
          GO TO 940
        END IF
  
  ! To determine the trough region based on x limits
  
        IF(axpr < xmin(nq1).OR.axpr > xmax(nq1)) THEN
          GO TO 940
        END IF
  
  ! Calculate xtran
  
        xc = (ABS(xmax(nq1)+xmin(nq1)))/2.0
        IF (xc >= axpr) THEN
          xtran = twidth(nq1)/2.0 - ABS(axpr-xc)
        ELSE
          xtran = twidth(nq1)/2.0 + ABS(axpr-xc)
        END IF
  
  !  Calculate ztran
  
    131 ztran = rdist
        !IF(ztran>delh(mcb,nq1) THEN
        ! GO TO 940
        !ENDIF
  
  ! Calculate slope of a line perpendicular to the central line passing through the point (aypr, z)
  
        aslope=-1/slope(nq1)
        ayint=aypr-aslope*azpr
  
  ! To find the intersection of these two lines
  
        yintr=(yint(nq1)*aslope-ayint*slope(nq1)) &
             /(aslope-slope(nq1))
        zintr=(ayint-yint(nq1))/(slope(nq1)-aslope)
  
  !  Calculate ytran
  
        llength=SQRT((yintr-cymin(nq1))**2+(zintr-czmax(nq1))**2)
  !	  IF(llen(1,nq1)<=0.67*llenm) THEN
  !	   IF(llength>0.67*llen(1,nq1))THEN
  !	    ytran=0.67*llenm+(abs(0.67*llen(1,nq1)-llength))
  !		GOTO 499
  !	   ELSE
  !	     GOTO 399
  !	   ENDIF
  !	  ELSE
  !	     GOTO 399
  !	  ENDIF
  ! 399  CONTINUE
  !      ytran=llength !+(abs(llen(mcb,nq1)-llenm))
  ! 499  CONTINUE
       ytran=llength
        nneq=k
        nl=j
        mloc=nq1
  !      IF(mloc==0) THEN
  !       write(9872,*)lcnum(mcb,mub,nneq,nl)
  !      ENDIF
  !vf      if(ind(ik,nq1)==4) then
  !vf       write(534,*)xtran,ytran,ztran,mub,k,j
  !vf      endif
        GO TO 999
  
    940 CONTINUE
  
        iterr=3
  
    999 CONTINUE
  
        RETURN
        END SUBROUTINE tranf
  
  !***********************************************************************
  
  !     Subroutine COTRANM transforms an (x,y) coordinate point in the
  !     real space to a point (x',y') in the straight-centerline
  !     geometric space.
  
  !***********************************************************************
  
        SUBROUTINE cotranm(x,y,xpr,ypr,x1,y1,alpha,width,length,indc,ierr)
  
        REAL*4, INTENT(IN)                         :: x
        REAL*4, INTENT(IN)                         :: y
        REAL*4, INTENT(OUT)                        :: xpr
        REAL*4, INTENT(OUT)                        :: ypr
        REAL*4, INTENT(IN)                         :: x1
        REAL*4, INTENT(IN)                         :: y1
        REAL*4, INTENT(IN OUT)                     :: alpha
        REAL*4, INTENT(IN)                         :: width
        REAL*4, INTENT(IN)                         :: length
        REAL*4                                     :: rtest
        REAL*4                                     :: xt
        REAL*4                                     :: yt
        REAL*4                                     :: xint
        REAL*4                                     :: yint
        REAL*4                                     :: acl
        REAL*4                                     :: bcl
        REAL*4                                     :: ap
        REAL*4                                     :: bp
        REAL*4                                     :: xl
        REAL*4                                     :: xr
        REAL*4                                     :: xtdist
        REAL*4                                     :: ytdist
        REAL*4                                     :: xtest
        REAL*4                                     :: ytest
        REAL*4                                     :: yy
        INTEGER*4, INTENT(IN OUT)                  :: indc
        INTEGER*4, INTENT(OUT)                     :: ierr
        INTEGER*4                                  :: narc
        INTEGER*4                                  :: jj
        INTEGER*4                                  :: jerr1
        INTEGER*4                                  :: jerr2
        INTEGER*4                                  :: indx
        INTEGER*4                                  :: indflag
   
      1 FORMAT(1X,a)
    101 FORMAT(1X,/)
   
  ! Slope and intercept of transformed center line
  
        xtdist =length*SIN(alpha)
        IF (indc == 1)THEN
          xt = x1-xtdist
        ELSE
          xt = x1+xtdist
        END IF
        yt = y1+length*COS(alpha)
        
        CALL line(x1,y1,xt,yt,acl,bcl,indx)
  
  !  Slope/intercept of a line (Perpendicular to trasformed center line
  !  and passing through the point (X,Y)
  
        ap = -1.0/acl
        bp = y-ap*x
  
  !      Check the point is inside CB
  
        xint = (bp-bcl)/(acl-ap)
        yint = ap*xint+bp
        xtest = xydist(x,y,xint,yint)
        ytest = xydist(xint,yint,x1,y1)
        IF (xtest < width/2.0)THEN
          IF(yint > y1)THEN
            IF (ytest < length)THEN
              ierr =1
              GO TO 10
            END IF
          END IF
        END IF
        ierr = 2
        GO TO 999
  
  ! Xprime and Yprime calculation
  
     10 CONTINUE
        yy = x*acl+bcl
        IF (indc == 1)THEN
          IF (yy < y)THEN
            xpr = width/2.0 + xtest
          ELSE
            xpr = width/2.0 - xtest
          END IF
        ELSE
          IF (yy > y)THEN
            xpr = width/2.0 + xtest
          ELSE
            xpr = width/2.0 - xtest
          END IF
        END IF
        ypr = ytest
  
    999 CONTINUE
  
        RETURN
        END SUBROUTINE cotranm
  
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
  
        END FUNCTION icount_digits
  
  !***********************************************************************
  !  Function check checks the status of error flags
  !***********************************************************************
  
        SUBROUTINE check(adum,ierr)
  
        INTEGER*4 :: ierr
        CHARACTER (LEN=10) :: adum
  
        IF (ierr /= 0) THEN
          WRITE(*,*) 'Error in allocating array or MPI call:  ',trim(adum)
          CALL MPI_Finalize(ierr)
          STOP
        ENDIF
  
        END SUBROUTINE check
  
  !***********************************************************************
  !***********************************************************************
  !
  !     Subroutine DTRAN transforms a point in the global coordinate to a point
  !                      in a local coordinate for the donor meduim scale sets.
  !
  !***********************************************************************
  
  SUBROUTINE dtran(cord, tcord, min, max, len)
  
  
  
  REAL*4, INTENT(IN OUT)                     :: cord
  REAL*4, INTENT(OUT)                        :: tcord
  REAL*4, INTENT(IN OUT)                     :: min
  REAL*4, INTENT(IN OUT)                     :: max
  REAL*4, INTENT(IN OUT)                     :: len
  !INTEGER*2, INTENT(OUT)                     :: icerr
  INTEGER*2                                  :: icord, iratio
  
  call chcrt (cord, min, max, icord)
  !
   if(icord==0) then
     goto 10
   else if (icord==1) then
         cord1=abs(cord)+min
         if(cord1<= max) then
               goto 10
             else
               goto 11
         endif
   else if( icord==2) then
        goto 11
   endif
  11 call scale(cord, len, iratio)
     tcord=cord-len*iratio
    ! write(541,*)tcord, cord, len, iratio
     goto 888
  10 tcord=cord
  888 continue
  RETURN
  END SUBROUTINE dtran
  !***********************************************************************
  !
  !     Subroutine CHECK determines criteria passed by a point
  !
  !***********************************************************************
  
  SUBROUTINE chcrt (cord2, min1, max1,lflag)
  
  
  
  REAL*4, INTENT(IN OUT)                     :: cord2
  !REAL*4, INTENT(OUT)                        :: tcord
  REAL*4, INTENT(IN OUT)                     :: min1
  REAL*4, INTENT(IN OUT)                     :: max1
  !REAL*4, INTENT(IN OUT)                     :: len
  INTEGER*2                                  :: lflag
  
  if ((cord2>=min1).and.(cord2<=max1)) then
  lflag=0
  else if(cord2<min1) then
  lflag=1
  else if (cord2> max1) then
  lflag=2
  endif
  887 continue
  RETURN
  END SUBROUTINE chcrt
  !***********************************************************************
  !
  !     Subroutine SCALE determines the location of a global coordinate
  !                      point in a local coordinate system
  !
  !***********************************************************************
  
  SUBROUTINE scale(cord3, len1, lratio1)
  
  
  
  REAL*4, INTENT(IN OUT)                     :: cord3
  !REAL*4, INTENT(OUT)                        :: ratio1
  REAL*4                                     :: rratio
  REAL*4, INTENT(IN OUT)                     :: len1
  INTEGER*2                                  :: ltr
  INTEGER*2, INTENT(OUT)                     ::lratio1
  !
  !
  !
   rratio= abs(cord3/len1)
  ltr=INT(rratio)
   if(ltr>1.and.ltr>rratio) then
   lratio1=ltr-1
   else
   lratio1=ltr
   endif
   !write(601,*)rratio, ltr, lratio1, cord3, len1
  RETURN
  END SUBROUTINE scale
  !***************************************************************************************************************************
  !***********************************************************************
  
  !     Subroutine NORGEN uses a set of subroutines to generate a single
  !        random variable from a normal distribution with mean and
  !        variance as specified in the call statement.  The
  !        subroutine TSEED should be called before the first
  !        call of NORGEN in a main program, in order to set the
  !        initial seed values from the clock.
  
  !***********************************************************************
  
  SUBROUTINE norgen(x,idt,mean,var,lind,is1,is2,is3)
  
  
  REAL*4, INTENT(OUT)                        :: x
  INTEGER*2, INTENT(IN OUT)                  :: idt,lind
  REAL*4, INTENT(IN)                         :: mean
  REAL*4, INTENT(IN OUT)                     :: var
  INTEGER*4, INTENT(IN OUT)                  :: is1
  INTEGER*4, INTENT(IN OUT)                  :: is2
  INTEGER*4, INTENT(IN OUT)                  :: is3
  !INTEGER*4, PARAMETER :: nidt=20
  REAL*8                                     :: xarr(2)
  INTEGER*2                                  :: n
  
  !     Number of uniform random variables to be generated = 2
  
  n = 2
  
  !     First, AS183 is called to generate two U[0,1] numbers
  !        (and to modify the seed values)
  
  CALL as183(is1,is2,is3,n,xarr)
  
  !     Next, BOXMUL is called to perform the Box-Muller transformation,
  !        which transforms the uniform variables into standard normal
  !        variables.
  
  CALL boxmul(n,xarr,xarr)
  
  !     Then, transform the standard normal into the distribution with
  !        specified mean and variance.
  
  x = xarr(1)*SQRT(var) + mean
  RETURN
  END SUBROUTINE norgen
  !***********************************************************************
  
  SUBROUTINE tseed(is1,is2,is3)
  
  !  TSEED constructs three integer seed values for the random number
  !  generator AS183 using the internal clock of the computer (PC/AT).
  
  !  The arguments must be type INTEGER.
  
  !  TSEED requires subroutine GETTIM, which is supplied by PROFORT.LIB.
  
  !  Written 8/85 by Tom Black
  
  
  INTEGER*4, INTENT(OUT)                     :: is1
  INTEGER*4, INTENT(OUT)                     :: is2
  INTEGER*4, INTENT(OUT)                     :: is3
  REAL :: r(3)
  doubleprecision dseed
  DOUBLE PRECISION :: d2p31m,d2p31
  !                                  D2P31M=(2**31) - 1
  !                                  D2P31 =(2**31)
  DATA               d2p31m/2147483647.d0/
  DATA               d2p31/2147483648.d0/
  
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
  !***********************************************************************
  
  SUBROUTINE boxmul(n,unif,norml)
  
  !  BOXMUL performs the standard Box-Muller normal transfomation on the
  !  rectangular [0,1] random numbers in UNIF and places them in NORML.
  !  N is the number of generated numbers in array UNIF.  UNIF and NORML
  !  have assumed dimensions corresponding to those in the calling program.
  
  !  THE TEMP STRUCTURE PERMITS UNIF AND NORML TO BE THE SAME VECTOR.
  
  !  N is type INTEGER.
  
  !  Arrays corresponding to UNIF and NORML must be single-precision REAL
  !  in the calling program.
  
  !  BOXMUL requires no FUNCTIONs or SUBROUTINEs.
  
  !  Written 3/85 by Tom Black.  REVISED 4/86.
  
  
  
  INTEGER*2, INTENT(IN)                      :: n
  REAL*8, INTENT(IN)                         :: unif(*)
  REAL*8, INTENT(OUT)                        :: norml(*)
  
  
  100 DO  i=1,n-1,2
    temp1=unif(i)
    temp2=unif(i+1)
    norml(i)=SQRT(0.0-2.0*LOG(temp1))*COS(6.28318*temp2)
    norml(i+1)=SQRT(0.0-2.0*LOG(temp1))*SIN(6.28318*temp2)
  END DO
  RETURN
  END SUBROUTINE boxmul
  !***********************************************************************
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
  INTEGER*2, INTENT(IN)                      :: n
  REAL*8, INTENT(OUT)                        :: unif(*)
  real*4:: num1, num2
  
  
  100 DO  i=1,n
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
  !       write(315,*)unif(i), num1
    endif
  END DO
  RETURN
  END SUBROUTINE as183
  !***********************************************************************
  
  !     Subroutine expgen uses a set of subroutines to generate a single
  !        random variable from a exponential distribution with mean and
  !        variance as specified in the call statement.  The
  !        subroutine TSEED should be called before the first
  !        call of expgen in a main program, in order to set the
  !        initial seed values from the clock.
  
  !***********************************************************************
  
  
  SUBROUTINE expgen(x,idt,mean,var,lind,is1,is2,is3)
  
  
  REAL*4                        :: x
  INTEGER*2, INTENT(IN OUT)                  :: idt,lind
  REAL*4, INTENT(IN)                         :: mean
  REAL*4, INTENT(IN OUT)                     :: var
  INTEGER*4, INTENT(IN OUT)                  :: is1
  INTEGER*4, INTENT(IN OUT)                  :: is2
  INTEGER*4, INTENT(IN OUT)                  :: is3
  INTEGER*4, PARAMETER :: nidt=20
  REAL*8 :: xarr(2)
  INTEGER*2 :: n
  !integer*2                        :: lind
  
  !     Number of uniform random variables to be generated = 2
  
  n = 2
  
  !     First, AS183 is called to generate two U[0,1] numbers
  !        (and to modify the seed values)
  
  CALL as183(is1,is2,is3,n,xarr)
  
  !     Next, exptransf is called to perform the Exponential transformation,
  !        which transforms the uniform variables into standard exponential
  !        variables.
  
  CALL exptransf(n,xarr,xarr,var)
  
  !     Then, transform the standard exponential into the distribution with
  !        specified mean and variance.
  
  x = xarr(1)*SQRT(var) + mean
  !x = xarr(1)
  
  RETURN
  END SUBROUTINE expgen
  !***********************************************************************
  
  SUBROUTINE exptransf(n,unif,expdev, vari)
  
  !  BOXMUL performs the standard Box-Muller normal transfomation on the
  !  rectangular [0,1] random numbers in UNIF and places them in NORML.
  !  N is the number of generated numbers in array UNIF.  UNIF and NORML
  !  have assumed dimensions corresponding to those in the calling program.
  
  !  THE TEMP STRUCTURE PERMITS UNIF AND NORML TO BE THE SAME VECTOR.
  
  !  N is type INTEGER.
  
  !  Arrays corresponding to UNIF and NORML must be single-precision REAL
  !  in the calling program.
  
  !  BOXMUL requires no FUNCTIONs or SUBROUTINEs.
  
  !  Written 3/85 by Tom Black.  REVISED 4/86.
  
  
  
  INTEGER*2, INTENT(IN)                      :: n
  REAL*8, INTENT(IN)                         :: unif(*)
  REAL*8, INTENT(OUT)                        :: expdev(*)
  Real*4                                     :: temp1, temp2
  Real*4, INTENT(IN OUT)                          ::vari
  
  
  104 DO  i=1,n-1,2
    temp1=unif(i)
    temp2=unif(i+1)
    expdev(i)=-(LOG(temp1))
    expdev(i+1)=-(LOG(temp2))
  !  write(164,*)expdev(i),expdev(i+1)
  END DO
  RETURN
  END SUBROUTINE exptransf
  !***********************************************************************
  !********************************************************************************************************
        SUBROUTINE chknam(str,LEN)
  !*************************************************************************************************
  !                   Check for a Valid File Name
  !                   ***************************
  !
  ! This subroutine takes the character string "str" of length "len" and
  ! removes all leading blanks and blanks out all characters after the
  ! first blank found in the string (leading blanks are removed first).
  !**********************************************************************************************************
  
            INTEGER, PARAMETER                       :: maxlen=132
            CHARACTER (LEN=41), INTENT(OUT)          :: str(maxlen)
            INTEGER, INTENT(IN)                      :: LEN
  
  ! Remove leading blanks:
  
             DO i=1,LEN-1
              IF(str(i) /= ' ') THEN
               IF(i == 1) GO TO 1
               DO j=1,LEN-i+1
                 k = j + i - 1
                 str(j) = str(k)
               END DO
               DO j=LEN,LEN-i+2,-1
                 str(j) = ' '
               END DO
               GO TO 1
              END IF
            END DO
            1 CONTINUE
  
  ! Find first blank and blank out the remaining characters:
  
            DO i=1,LEN-1
             IF(str(i) == ' ') THEN
               DO j=i+1,LEN
                str(j) = ' '
               END DO
               GO TO 2
             END IF
           END DO
           2 CONTINUE
  
  ! Return with modified file name:
  
           RETURN
           END SUBROUTINE chknam
  !**********************************************************************************************************************
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
  