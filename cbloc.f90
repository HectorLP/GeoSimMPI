!********************************************************************************
! Date: 2007-07-27  Time: 15:05:28
 
! CBLOC.f90, the code generates points (locates level III unit types) randomly in a domain.
!    Output files:  X,Y,Z coordinates for each point
!                   CBLWH.out (All length,width,thickness of each level III polyhedron)
!                    CBC.out (Cross bar channel fill info)

!Last modified 06/16/2011 RWR
!*********************************************************************************
REAL*4          :: x,y,z,xi,yi,zi,xmx,ymx,zmx,wcbmn,wcbvar,  &
    lcbmn,lcbvar,hcbmn,hcbvar,xs,ys,  &
    xxmn,xxvar,yymn,yyvar,alpha,xcbc1,ycbc1,xcbc2,ycbc2, xcbc3,ycbc3, &
    wcb,lcb,hcb,yy,xx,zz,zzmn,zzvar,alphamn,alphavar,a,b,xxc,wcbc,zcbc,lcbc
INTEGER*2       :: iflag ,m,l,nflag,ind,ninfl,iinf
integer*4       :: is1,is2,is3
!CHARACTER (LEN=12) :: filein, fileot, fileot1
!DATA      lin/1/
CHARACTER (LEN=30) :: filein, fileot, fileout1, fileout2, fileout3
DATA      lin/1/
CHARACTER (LEN=4)   :: intform !Naum
DATA intform / '(I )' / !Naum

!10-29-08

INTEGER*2             :: ians
INTEGER*4             :: icountl, icountw, icounth
INTEGER*4             :: icountcl,icountcw,icountch
CHARACTER(LEN=40)     :: strl, strw, strh
CHARACTER(LEN=40)     :: strcl, strcw, strch

! 07-20-10 RR
REAL*4                :: xtemp
REAL*4                :: ytemp
REAL*4                :: ztemp

1     FORMAT(1X,a)
11   FORMAT(a12)
101 FORMAT(1X,/)
32  FORMAT(2(i4,2X),4(f13.2,2X),i4)
31  FORMAT(2(i4,2X),2X,3(f13.2),2X)
33  FORMAT(6(f13.2),i4,4(f13.2))
41  FORMAT(i4,2X,2(f15.5,2X))

!     FILENAME READS AND OPEN STATEMENTS

WRITE(*,101)
WRITE(*,101)
WRITE(*,101)
WRITE(*,1) 'Program CBLOC'
WRITE(*,1) '**************'
WRITE(*,101)
WRITE(*,1) 'Name of input file?'
READ(*,11) filein

      fileot='cbpts.out'

WRITE(*,*)'Do you want to store the IGL that are generated? Enter 1 (yes) or 0 (no)'
WRITE(*,*)' (Answer usually is 0 (no) )'

  READ(*,*)ians
  
  IF(ians==1) THEN
     WRITE(*,*)'Enter the file for length'
     READ(*,*)strl
 
     WRITE(*,*)'Enter the file for width'
     READ(*,*)strw

     WRITE(*,*)'Enter the file for height'
     READ(*,*)strh
     
     WRITE(*,*)'Enter the file for cbc_length'
     READ(*,*)strcl
 
     WRITE(*,*)'Enter the file for cbc_width'
     READ(*,*)strcw

     WRITE(*,*)'Enter the file for cbc_height'
     READ(*,*)strch
  ENDIF
!

OPEN(UNIT=4,FILE=filein,STATUS='old')
OPEN(UNIT=5,FILE=fileot,FORM='unformatted',STATUS='unknown')
OPEN(UNIT=6,FILE='cblwh.out',FORM='unformatted',STATUS='unknown')
!OPEN(UNIT=8,FILE='cbc.out',FORM='unformatted',STATUS='unknown')
OPEN(UNIT=9,FILE='cbclw.out',FORM='unformatted',STATUS='unknown')
OPEN(UNIT=99,FILE='maxv.out',STATUS='unknown')

!10-29-08
 IF(ians==1) THEN
   OPEN(UNIT=10, FILE=strl, STATUS='unknown')
   OPEN(UNIT=11, FILE=strw, STATUS='unknown')
   OPEN(UNIT=12, FILE=strh, STATUS='unknown')
   OPEN(UNIT=15, FILE=strcl, STATUS='unknown')
   OPEN(UNIT=16, FILE=strcw, STATUS='unknown')
   OPEN(UNIT=17, FILE=strch, STATUS='unknown')
   OPEN(UNIT=18, FILE='cbccnt.out', STATUS='unknown')
   OPEN(UNIT=14, FILE='cbcnt.out', STATUS='unknown')
 ENDIF
! INITIALIZE
  icountl=0
  icountw=0
  icounth=0
  icountcl=0
  icountcw=0
  icountch=0
!
READ(4,*)xi,yi,zi !starting location of first CB
READ(4,*)wcbmn,wcbvar,lcbmn,lcbvar,hcbmn,hcbvar
READ(4,*)yymn,yyvar,xxmn,xxvar,zzmn,zzvar
READ(4,*)xmx, ymx, zmx
READ(4,*)alphamn, alphavar

      ztemp=zi
!     Set initial values of random number generator seeds, and
!        print to output:

WRITE(*,101)
WRITE(*,1) 'Are seed values to be specified (Enter 1) or randomly generated (Enter 0)? '
READ(*,*) ans
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


!     Loop with each layer


nflag=0
!1. What does this 'l=1' mean?
l=1
m=0
z = zi
!2. Why does xi and yi are shifted -100 and -200? (Is there any standard to choose these values?)
xs=xi-100                                                            !WSU 4/4/08
ys=yi-200                                                           !WSU 4/4/08
ytemp=ys      ! 07-20-10 RR
xtemp=xs       ! 07-20-10 RR


10   IF (z < (zmx+2.0))THEN                                         !WSU 4/4/08
  iflag=1
!          Z = Zs
            IF(xs>xtemp/2.0) THEN           ! 07-20-10 RR
             xs=xtemp                       ! 07-20-10 RR
            ENDIF                           ! 07-20-10 RR
            IF(ys>ytemp/2.0) THEN           ! 07-20-10 RR
             ys=ytemp                       ! 07-20-10 RR
            ENDIF                           ! 07-20-10 RR

  x = xs
  y = ys
  20     IF(x < xmx)THEN
    m=m+1
    WRITE(intform(3:3),'(I1)') icount_digits(m-1)
    fileout2(1:8) = 'cbc.out.' !Naum
    WRITE(fileout2(9:),intform) m-1 !Naum
    OPEN(UNIT=8,FILE=fileout2,FORM='unformatted',STATUS='unknown') !Naum
    100     CONTINUE
    CALL norgen(wcb,wcbmn,wcbvar,is1,is2,is3)
    IF (wcb <= 0.0) GO TO 100
    IF(ians==1) THEN
     write(11,*)wcb
     icountw=icountw+1
    ENDIF
    110     CONTINUE
    CALL norgen(lcb,lcbmn,lcbvar,is1,is2,is3)
    IF (lcb <= 0.0) GO TO 110

     
    IF(ians==1) THEN
     write(10,*)lcb
     icountl=icountl+1
    ENDIF
    120     CONTINUE
    CALL norgen(hcb,hcbmn,hcbvar,is1,is2,is3)
    IF (hcb <= 0.0) GO TO 120

     
    IF(ians==1) THEN
     write(12,*)hcb
     icounth=icounth+1
    ENDIF
    CALL norgen(yy,yymn,yyvar,is1,is2,is3)
    CALL norgen(xx,xxmn,xxvar,is1,is2,is3)
    CALL norgen(zz,zzmn,zzvar,is1,is2,is3)
    150     CONTINUE
    CALL norgen(alpha,alphamn,alphavar,is1,is2,is3)
    IF (alpha <= 0.0) GO TO 150

    
    CALL norgen1(ind,is1,is2,is3)

    
    
    
    WRITE(6)m,l,wcb,lcb,hcb,alpha,ind !cblwh.out 
    !WRITE(66,*)wcb,lcb,hcb 
    WRITE(5)m,l,x,y,z !cbpts.out
	!WRITE(67,*)m,l,x,y,z,wcb,lcb,hcb,alpha,ind
    !3.Does ind = 1 or 2 mean the boundary between cross bar channel and compound bar? 
    !4.Where can we find the details about the calculation below? Although the paper (Ramya et.al., 2010)
    ! mentioned 'The modules CBLOC and XBGEN locate a cross-bar channel fill polyhedron created by XBPLANE into
    !a compound bar deposit. These have archetypal polyhedra with a simple convex-up geometry formed by five planes.
    !locator points and series of arcs are created in a procedure similar to the way UBLOC and UBGEN locate the arcs
    !of unit bars.' For example, how is the value of a (-0 & 0.22) chosen?
    IF (ind==1)THEN    
     xcbc1=xi+wcb/2.0
     ycbc1=yi
     xcbc2=xi
     ycbc2=ycbc1+lcb/2.0
     xcbc3=xi-wcb                                                    !WSU 4/4/08
     ycbc3=ycbc1+lcb+100
     a =-0
    ELSE
     xcbc1=xi-wcb/2.0
     ycbc1=yi
     xcbc2=xi
     ycbc2=ycbc1+lcb/2.0
     xcbc3=xi+wcb                                                    !WSU 4/4/08
     ycbc3=ycbc1+lcb+100
     a = 0.22
    ENDIF
    b = ycbc1-a*xcbc1                                              !WSU 4/4/08
    ninfl=3
    iinf=0
    xxc = 0.0
    wcbc = wcb/15.0                                                !WSU 4/4/08
    zcbc = z-hcb 
    lcbc=lcb+50                                                  !WSU 4/4/08
    WRITE(9)wcbc,lcbc,zcbc !cbclw.out
    !WRITE(91,*)wcbc,lcbc,zcbc !cbclw.out
    IF(ians==1) THEN
       WRITE(16,*)wcbc
       WRITE(15,*)lcbc
       WRITE(17,*)zcbc
       icountcl=icountcl+1
       icountcw=icountcw+1
       icountch=icountch+1
    ENDIF
    !WRITE(8)m,ninfl !cbc.out
    !WRITE(8)xcbc1,ycbc1,xcbc2,ycbc2,xcbc3,ycbc3,iinf,xx,a,b,wcbc
    !This is for parallel program
    WRITE(8)xcbc1,xcbc2,xcbc3,ycbc1,ycbc2,ycbc3,a,b,wcbc  !Naum 
    !WRITE(68,*)xcbc1,xcbc2,xcbc3,ycbc1,ycbc2,ycbc3,a,b,wcbc  !Naum  
    y = y+yy
   ! WRITE(33,*)zz
    IF (ind==1)THEN                                                !WSU 4/4/08                                               
    z = z-zz                                                       !WSU 4/4/08
    ELSE                                                           !WSU 4/4/08
    z = z+zz
    ENDIF                                                          !WSU 4/4/08

    IF (l==1)then                                                  !WSU 4/4/08           
      x = x+wcb                                                      !WSU 4/4/08
    ELSE                                                         !WSU 4/4/08
      x = x+(wcb/3.00)                                               !WSU 4/4/08
    ENDIF                                                        !WSU 4/4/08

    GO TO 20
  ELSE
    IF(y < (ymx-lcbmn/2.0))THEN                                      !WSU 4/4/08

      x = xs - xx                                                  !WSU 4/4/08

      IF (l==1)then                                                !WSU 4/4/08
        y = y + lcb                                                  !WSU 4/4/08
      ELSE                                                         !WSU 4/4/08
        y = y+(lcb/3.00)
      ENDIF                                                        !WSU 4/4/08
        GO TO 20

    ELSE
      xs = xs+xx
      IF (l==1)THEN                                                 !WSU 4/4/08
        zi = zi+hcbmn/4.00                                            !WSU 4/4/08
      ELSE                                                          !WSU 4/4/08
        zi = zi+hcbmn/2.00                                            !WSU 4/4/08
      ENDIF                                                         !WSU 4/4/08
      z = zi                                                        !WSU 4/4/08
      ys=ys+yy
      l=l+1
	  !write(9998,*)xs,ys,z,m,l   ! 07-20-10 RR

      GO TO 10
    END IF
  END IF
ELSE
  GO TO 30
END IF
30   CONTINUE
  WRITE(99,*)m
IF(ians==1) THEN
  WRITE(14,*)icountl
  WRITE(14,*)icountw
  WRITE(14,*)icounth
  
  WRITE(18,*)icountcl
  WRITE(18,*)icountcw
  WRITE(18,*)icountch
ENDIF
STOP
END
!***********************************************************************

SUBROUTINE tseed(is1,is2,is3)

!  TSEED constructs three integer seed values for the random number
!  generator AS183 using the internal clock of the computer.

!  The arguments must be type INTEGER.

!  TSEED requires subroutine GETTIM, which is supplied by PROFORT.LIB.

!  Written 8/85 by Tom Black


INTEGER*4, INTENT(OUT)                     :: is1
INTEGER*4, INTENT(OUT)                     :: is2
INTEGER*4, INTENT(OUT)                     :: is3
REAL :: r(3)
doubleprecision dseed
DOUBLE PRECISION :: d2p31m,d2p31

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
REAL*4, INTENT(OUT)                        :: unif(*)


100 DO  i=1,n
  ix=171*MOD(ix,177)-2*(ix/177)
  iy=172*MOD(iy,176)-35*(iy/176)
  iz=170*MOD(iz,178)-63*(iz/178)
  IF(ix < 0) ix=ix+30269
  IF(iy < 0) iy=iy+30307
  IF(iz < 0) iz=iz+30323
  unif(i)=AMOD(FLOAT(ix)/30269.0+FLOAT(iy)/30307.0+ FLOAT(iz)/30323.0,1.0)
END DO
RETURN
END SUBROUTINE as183
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
REAL*4, INTENT(IN)                         :: unif(*)
REAL*4, INTENT(OUT)                        :: norml(*)


100 DO  i=1,n-1,2
  temp1=unif(i)
  temp2=unif(i+1)
  norml(i)=SQRT(0.0-2.0*LOG(temp1))*COS(6.28318*temp2)
  norml(i+1)=SQRT(0.0-2.0*LOG(temp1))*SIN(6.28318*temp2)
END DO
RETURN
END SUBROUTINE boxmul
!***********************************************************************

!     Subroutine NORGEN uses a set of subroutines to generate a single
!        random variable from a normal distribution with mean and
!        variance as specified in the call statement.  The
!        subroutine TSEED should be called before the first
!        call of NORGEN in a main program, in order to set the
!        initial seed values from the clock.

!***********************************************************************

SUBROUTINE norgen(x,mean,var,is1,is2,is3)


REAL*4, INTENT(OUT)                        :: x
REAL*4, INTENT(IN)                         :: mean
REAL*4, INTENT(IN OUT)                     :: var
INTEGER*4, INTENT(IN OUT)                  :: is1
INTEGER*4, INTENT(IN OUT)                  :: is2
INTEGER*4, INTENT(IN OUT)                  :: is3
REAL*4 :: xarr(2)
INTEGER*2 :: n

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

!     Subroutine NORGEN1 uses a set of subroutines to generate a single
!        random variable from a Uniform distribution

!***********************************************************************

SUBROUTINE norgen1(ind,is1,is2,is3)


INTEGER*2, INTENT(OUT)                     :: ind
INTEGER*4, INTENT(IN OUT)                  :: is1
INTEGER*4, INTENT(IN OUT)                  :: is2
INTEGER*4, INTENT(IN OUT)                  :: is3
REAL*4 :: x,xarr(2)
INTEGER*2 :: n

!     Number of uniform random variables to be generated = 2

n = 2

!     First, AS183 is called to generate two U[0,1] numbers
!        (and to modify the seed values)

CALL as183(is1,is2,is3,n,xarr)


x = xarr(1)
IF (x < 0.3333)THEN
!     if (x.gt.0.66)then
  ind = 1
!     else
!      ind = 2
!     endif
ELSE
  ind = 2
END IF

RETURN
END SUBROUTINE norgen1
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
