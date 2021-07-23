! compute the isotopic mass distribution for
! given chemical formula (nat atoms with OZ it())
! returns nsig masses in mass() with probability
! mint()
! this is a quick and dirty MC algorithm and its run-time depends
! critically on the number of trials nrnd
! TK currently 10 elements covered:
! TK 1H, 6C, 7N, 8O, 9F, 14Si, 15N, 16S, 17Cl, 26Fe
! TK increased mass acuracy values, currently unit masses
! TK See: ATOMIC WEIGHTS OF THE ELEMENTS: REVIEW 2000 (IUPAC Technical Report)
! TK Pure Appl. Chem., Vol. 75, No. 6, pp. 683�800, 2003.
! SW calculate accurate mass
     subroutine isotope(nat,it,ndim,nrnd,rnd,mass,mint,nsig, &
       dict,mdict,key_list,formula)
     use dictionary_m
     implicit none
     character(len=16), allocatable :: key_list(:)
     character(len=16) charstrk, charstrv
     type(dictionary_t):: dict,mdict
     real*16 intensity
     integer nat,it(*),nsig,nrnd,ndim, s, dictsize
     integer mass(*)
     real*8  mint(*)
     real*4  rnd(ndim,nrnd)
     character*80 formula
     real*8  prob(200,4),massiso(200,4),p1,p2,x,xmass,examass
     integer niso(200)
     integer nmass(10000)
     integer n,i,j,iso,imass,isum,k,iti
     real*8 r,xm

     niso=0
     prob=0
     massiso=0
! SW accuratemass
! TK 1H
     niso(1)=2
     prob(1,1)=0.0115
     prob(1,2)=99.9885
     massiso(1,1)=2.014102
     massiso(1,2)=1.007825

! TK 6C
     niso(6)=2
     prob(6,1)=1.07
     prob(6,2)=98.93
     massiso(6,1)=13.003355
     massiso(6,2)=12.000000

! TK 7N
     niso(7)=2
     prob(7,1)=0.368
     prob(7,2)=99.632
     massiso(7,1)=15.000109
     massiso(7,2)=14.003074

! TK 8O (Oxygen)
     niso(8)=3
     prob(8,1)=0.038
     prob(8,2)=0.205
     prob(8,3)=99.757
     massiso(8,1)=16.999132
     massiso(8,2)=17.999160
     massiso(8,3)=15.994915

! TK 9F
     niso(9)=1
     prob(9,1)=100.
     massiso(9,1)=18.998403

! TK 14Si
     niso(14)=3
     prob(14,1)=92.223
     prob(14,2)=4.685
     prob(14,3)=3.092
     massiso(14,1)=27.976927
     massiso(14,2)=28.976495
     massiso(14,3)=29.973770

! SW 15P
     niso(15)=1
     prob(15,1)=100.
     massiso(15,1)=30.973762

! TK 16S
     niso(16)=4
     prob(16,1)=0.02
     prob(16,2)=0.76
     prob(16,3)=4.29
     prob(16,4)=94.93
     massiso(16,1)=35.967081
     massiso(16,2)=32.971458
     massiso(16,3)=33.967867
     massiso(16,4)=31.972071


! TK 17Cl - changed to isotopic masses and abundances OK
     niso(17)=2
     prob(17,1)=75.76
     prob(17,2)=24.24
     massiso(17,1)=34.968853
     massiso(17,2)=36.965903

! TK 26Fe
     niso(26)=4
     prob(26,1)=5.845
     prob(26,2)=91.754
     prob(26,3)=2.119
     prob(26,4)=0.2
     massiso(26,1)=53.939
     massiso(26,2)=55.934
     massiso(26,3)=56.935
     massiso(26,4)=57.933

! TK added 35Br
     niso(35)=2
     prob(35,1)=50.69
     prob(35,2)=49.31
     massiso(35,1)=78.91833
     massiso(35,2)=80.91629

! SW monoisotopic noble gas
! TK Helium (assume to be monoisotopic)
     niso(2)=1
     prob(2,1)=100.
     massiso(2,1)=4.0026
! TK Neon (assume to be monoisotopic)
     niso(10)=1
     prob(10,1)=100.
     massiso(10,1)=20.180
! TK Argon (assume to be monoisotopic)
     niso(18)=1
     prob(18,1)=100.
     massiso(18,1)=39.948
! TK Krypton (assume to be monoisotopic)
     niso(36)=1
     prob(36,1)=100.
     massiso(36,1)=83.798
! TK Xenon (assume to be monoisotopic)
     niso(54)=1
     prob(54,1)=100.
     massiso(54,1)=131.29
! SW new added
     !  28 Ni (Nickel)
            niso(28)=5
            prob(28,1)=68.08
            prob(28,2)=26.22
            prob(28,3)=1.14
            prob(28,4)=3.63
            prob(28,5)=0.93
            massiso(28,1)=57.935343
            massiso(28,2)=59.930786
            massiso(28,3)=60.931056
            massiso(28,4)=61.928345
            massiso(28,5)=63.927966

       ! 29 Cu (Copper)
            niso(29)=2
            prob(29,1)=69.15
            prob(29,2)=30.85
            massiso(29,1)=62.929597
            massiso(29,2)=64.927789

       ! 30 Zn (Zinc)
            niso(30)=5
            prob(30,1)=48.6
            prob(30,2)=27.9
            prob(30,3)=4.1
            prob(30,4)=18.8
            prob(30,5)=0.6
            massiso(30,1)=63.929142
            massiso(30,2)=65.926033
            massiso(30,3)=66.927127
            massiso(30,4)=67.924884
            massiso(30,5)=69.925319

     ! 27 Co (Cobalt)
          niso(27)=1
          prob(27,1)=100.000
          massiso(27,1)=58.933195

    ! 74 W(Tungsten)
          niso(74)=5
          prob(74,1)=0.12
          massiso(74,1)=179.946706
          prob(74,2)=26.50
          massiso(74,2)=181.948205
          prob(74,3)=14.31
          massiso(74,3)=182.9502242
          prob(74,4)=30.64
          massiso(74,4)=183.9509323
          prob(74,5)=28.43
          massiso(74,5)=185.954362

     prob = prob * 0.01
     dictsize=1024
     call dict%init(dictsize)




! TK mass currently only loops to element 36  (Krypton)
! check isotope data saved in this program
     do i=1,74
        xm=0
        do j=1,niso(i)
           xm=xm+prob(i,j)
        enddo
        if(niso(i).gt.0.and.abs(xm-1.).gt.0.01) &
        stop 'internal isotope error 1'
     enddo

     do i=1,nat
        if(it(i).gt.100)then
           niso   (it(i))  =1
           prob   (it(i),1)=1.
           massiso(it(i),1)=it(i)-100.
        endif
     enddo

     do i=1,nat
        if(niso(it(i)).eq.0) stop 'internal isotope error 2'
     enddo

     nmass=0
!      nrnd=200
     do n=1,nrnd
     xmass=0
     do i=1,nat
        iti=it(i)
        r=rnd(i,n)
        p1=0.0
        p2=prob(iti,1)
        do iso=1,niso(iti)
           if(r.ge.p1.and.r.le.p2)then
              x=massiso(iti,iso)
              exit
           endif
           p1=p2
           p2=p2+prob(iti,iso+1)
        enddo
        xmass=xmass+x

     enddo
     xmass = aint( xmass * 10000.0) / 10000.0

     !SW mass and intensity
     imass=int(xmass)
     nmass(imass)=nmass(imass)+1
     !SW exact mass and intensity
     write(charstrk,'(F16.4)')xmass
     call mdict%set(charstrk, formula)

   !  write(*,*) charstrk,'<key vaule>',dict%get(charstrk)
     !WRITE (STRING,*) variable.
     if (dict%get(charstrk) .eq. ' ') THEN
       charstrv='0.00002'
       call dict%set(charstrk, charstrv)
     else
       charstrv = dict%get(charstrk)
       read (charstrv,*)intensity
       !write (*,*)'intensity',intensity
       intensity = (intensity+0.00002) !normalization
       write(charstrv,'(F16.6)')intensity
     !  write (*,*)'charstrv',charstrv
       call dict%set(charstrk, charstrv)



     endif

     !SW write formula

     enddo


     isum=sum(nmass)
     k=0
     do i=1,2000
        if(nmass(i).gt.0)then
           k=k+1
           mass(k)=i
           mint(k)=float(nmass(i))/float(isum)

        endif
     enddo


     nsig=k
     end


! TK sub function as explained below
     REAL FUNCTION snorm()

!C**********************************************************************C
!C                                                                      C
!C                                                                      C
!C     (STANDARD-)  N O R M A L  DISTRIBUTION                           C
!C                                                                      C
!C                                                                      C
!C**********************************************************************C
!C**********************************************************************C
!C                                                                      C
!C     FOR DETAILS SEE:                                                 C
!C                                                                      C
!C               AHRENS, J.H. AND DIETER, U.                            C
!C               EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM             C
!C               SAMPLING FROM THE NORMAL DISTRIBUTION.                 C
!C               MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.          C
!C                                                                      C
!C     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  C
!C     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  C
!C                                                                      C
!C     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   C
!C     SUNIF.  The argument IR thus goes away.                          C
!C                                                                      C
!C**********************************************************************C
!C
     DIMENSION a(32),d(31),t(31),h(31)
!C
!C     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
!C     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
!C

     DATA a/0.0,.3917609E-1,.7841241E-1,.1177699,.1573107,.1970991, &
          .2372021,.2776904,.3186394,.3601299,.4022501,.4450965, &
          .4887764,.5334097,.5791322,.6260990,.6744898,.7245144, &
          .7764218,.8305109,.8871466,.9467818,1.009990,1.077516, &
          1.150349,1.229859,1.318011,1.417797,1.534121,1.675940, &
          1.862732,2.153875/
     DATA d/5*0.0,.2636843,.2425085,.2255674,.2116342,.1999243, &
          .1899108,.1812252,.1736014,.1668419,.1607967,.1553497, &
          .1504094,.1459026,.1417700,.1379632,.1344418,.1311722, &
          .1281260,.1252791,.1226109,.1201036,.1177417,.1155119, &
          .1134023,.1114027,.1095039/
     DATA t/.7673828E-3,.2306870E-2,.3860618E-2,.5438454E-2, &
          .7050699E-2,.8708396E-2,.1042357E-1,.1220953E-1,.1408125E-1, &
          .1605579E-1,.1815290E-1,.2039573E-1,.2281177E-1,.2543407E-1, &
          .2830296E-1,.3146822E-1,.3499233E-1,.3895483E-1,.4345878E-1, &
          .4864035E-1,.5468334E-1,.6184222E-1,.7047983E-1,.8113195E-1, &
          .9462444E-1,.1123001,.1364980,.1716886,.2276241,.3304980, &
          .5847031/
     DATA h/.3920617E-1,.3932705E-1,.3950999E-1,.3975703E-1, &
          .4007093E-1,.4045533E-1,.4091481E-1,.4145507E-1,.4208311E-1, &
          .4280748E-1,.4363863E-1,.4458932E-1,.4567523E-1,.4691571E-1, &
          .4833487E-1,.4996298E-1,.5183859E-1,.5401138E-1,.5654656E-1, &
          .5953130E-1,.6308489E-1,.6737503E-1,.7264544E-1,.7926471E-1, &
          .8781922E-1,.9930398E-1,.1155599,.1404344,.1836142,.2790016, &
          .7010474/

!C
  10 call random_number(u)
     s = 0.0
     IF (u.GT.0.5) s = 1.0
     u = u + u - s
  20 u = 32.0*u
     i = int(u)
     IF (i.EQ.32) i = 31
     IF (i.EQ.0) GO TO 100
!C
!C                                START CENTER
!C
  30 ustar = u - float(i)
     aa = a(i)
  40 IF (ustar.LE.t(i)) GO TO 60
     w = (ustar-t(i))*h(i)
!C
!C                                EXIT   (BOTH CASES)
!C
  50 y = aa + w
     snorm = y
     IF (s.EQ.1.0) snorm = -y
     RETURN
!C
!C                                CENTER CONTINUED
!C
  60 call random_number(u)
     w = u* (a(i+1)-aa)
     tt = (0.5*w+aa)*w
     GO TO 80

  70 tt = u
     call random_number(ustar)
  80 IF (ustar.GT.tt) GO TO 50
  90 call random_number(u)
     IF (ustar.GE.u) GO TO 70
     call random_number(ustar)
     GO TO 40
!C
!C                                START TAIL
!C
 100 i = 6
     aa = a(32)
     GO TO 120

 110 aa = aa + d(i)
     i = i + 1
 120 u = u + u
     IF (u.LT.1.0) GO TO 110
 130 u = u - 1.0
 140 w = u*d(i)
     tt = (0.5*w+aa)*w
     GO TO 160

 150 tt = u
 160 call random_number(ustar)
     IF (ustar.GT.tt) GO TO 50
 170 call random_number(u)
     IF (ustar.GE.u) GO TO 150
     call random_number(u)
     GO TO 140

     END

!     *****************************************************************
! TK  This routine handles
! TK  Called by: Main program
! TK  USES:      READAA

     SUBROUTINE READL(A1,X,N)
     IMPLICIT REAL*8 (A-H,O-Z)
     CHARACTER*(*) A1
     DIMENSION X(*)
     I=0
     IS=1
 10  I=I+1
     X(I)=READAA(A1,IS,IB,IE)
     IF(IB.GT.0 .AND. IE.GT.0) THEN
                               IS=IE
                               GOTO 10
     ENDIF
     N=I-1
     RETURN
     END

! *******
! SW function to  computes order of character strings
! SW Return -1 if C1 before C2, 0 if C1 = C2, and 1 if C1 after C2.
     integer(2) function cmp_function(a1, a2)
     character*16 a1, a2
     Real*16 b1, b2
     read (a1,*) b1
     read (a2,*) b2
     if (b1 < b2) THEN
       cmp_function = -1
     else
       cmp_function = 1
     end if
     end function



! *******
! SW function to get value as float
     FUNCTION floatv(charstrk, d)
     use dictionary_m
     implicit NONE
     CHARACTER*16 charstrk, charstrv
     Real*16 floatv,intensity
     type(dictionary_t) :: d
     character(len=16), allocatable :: key_list(:)
     charstrv = d%get(charstrk)
!      write(*,*)'charstrk',charstrk
!      write(*,*)'charstrv',charstrv
     read (charstrv,*) intensity
     floatv = intensity
     return
   end function floatv

!     *****************************************************************
! TK  This function handles
! TK  Called by: Main program
! TK  USES:      no other sub

     FUNCTION READAA(A,ISTART,IEND,IEND2)
     IMPLICIT REAL*8 (A-H,O-Z)
     REAL*8 READAA
     CHARACTER*(*) A
     NINE=ICHAR('9')
     IZERO=ICHAR('0')
     MINUS=ICHAR('-')
     IDOT=ICHAR('.')
     ND=ICHAR('D')
     NE=ICHAR('E')
     IBL=ICHAR(' ')
     IEND=0
     IEND2=0
     IDIG=0
     C1=0
     C2=0
     ONE=1.D0
     X = 1.D0
     NL=LEN(A)
     DO 10 J=ISTART,NL-1
        N=ICHAR(A(J:J))
        M=ICHAR(A(J+1:J+1))
        IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20
        IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO .OR. M.EQ.IDOT)) GOTO 20
  10 CONTINUE
     READAA=0.D0
     RETURN
  20 CONTINUE
     IEND=J
     DO 30 I=J,NL
        N=ICHAR(A(I:I))
        IF(N.LE.NINE.AND.N.GE.IZERO) THEN
           IDIG=IDIG+1
           IF (IDIG.GT.10) GOTO 60
           C1=C1*10+N-IZERO
        ELSEIF(N.EQ.MINUS.AND.I.EQ.J) THEN
           ONE=-1.D0
        ELSEIF(N.EQ.IDOT) THEN
           GOTO 40
        ELSE
           GOTO 60
        ENDIF
  30 CONTINUE
  40 CONTINUE
     IDIG=0
     DO 50 II=I+1,NL
        N=ICHAR(A(II:II))
        IF(N.LE.NINE.AND.N.GE.IZERO) THEN
           IDIG=IDIG+1
           IF (IDIG.GT.10) GOTO 60
           C2=C2*10+N-IZERO
           X = X /10
        ELSEIF(N.EQ.MINUS.AND.II.EQ.I) THEN
           X=-X
        ELSE
           GOTO 60
        ENDIF
  50 CONTINUE

!C
!C PUT THE PIECES TOGETHER
!C

  60 CONTINUE
     READAA= ONE * ( C1 + C2 * X)
     DO 55 J=IEND,NL
        N=ICHAR(A(J:J))
        IEND2=J
        IF(N.EQ.IBL)RETURN
  55 IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57
     RETURN

  57 C1=0.0D0
     ONE=1.0D0
     DO 31 I=J+1,NL
        N=ICHAR(A(I:I))
        IEND2=I
        IF(N.EQ.IBL)GOTO 70
        IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO
        IF(N.EQ.MINUS)ONE=-1.0D0
  31 CONTINUE
  61 CONTINUE
  70 READAA=READAA*10**(ONE*C1)
     RETURN
     END
