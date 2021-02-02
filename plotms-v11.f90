! ============================================================================
! Name        : plotms.f90
! Author      : S. Grimme, modifications T. Kind (FiehnLab 2013)
! Version     : 2.4  (Feb 19 2014)
! Copyright   :
! Description : plotms from QCEIMS
! ============================================================================

!====================================================================
!     Original PlotMS Program for extraction of mass spectra
!     from quantum mechanical simulations used within QCEIMS
!
!                      *********************************************
!                      *                S. Grimme                  *
!                      * Mulliken Center for Theoretical Chemistry *
!                      *             Universitaet Bonn             *
!                      *                  2008-13                  *
!                      *********************************************
!
!     Please cite as:
!     S.Grimme, Angew.Chem.Int.Ed. 52 (2013) 6306-6312.
!
!     Adapted version:
!     Tobias Kind - FiehnLab 2013
!
!     Changes:
!     1) Output MSP or JCAMP formated unit mass mass spectra (requires exp.dat and mass.agr)
!     2) Elements Br and W added
!     3) Chlorine masses changed to accurate masses
!     4) Output of exact masses
!
!     Adapted version:
!     Shunyang Wang - FiehnLab 2020
!     Shunyang Wang - Dec 2020
!
!====================================================================




      program plotms
      use IFPORT ! To get QSORT
      use dictionary_m ! get hash table
      use, intrinsic :: iso_c_binding, only: c_size_t
      implicit none
      integer(2), external :: cmp_function
      integer (C_SIZE_T) array_len, array_size
!     treat i,j,k,l,m,n as integers and all other cariables as real

!    declare parameters in dictionary
      type(dictionary_t) :: dict, tdict, mdict
      character(len=16), allocatable :: key_list(:), tkey_list(:),mkey_list(:)
      integer dictsize,s,si,pos
      character*16 charstrk, charstrv, ctmp
      real*16 ttmp, tmp, floatv, kmax, ftmp

      integer n,i,j,k,kk,kkk,kkkk,nn,natot,nsig
      integer nrnd,ndim,imax,imin,nagrfile,length
! maximum mass = 10000
      real*8  iexp (2,10000)
      integer iat  (10000)
      integer nat  (10000)
      integer idum (10000)
      integer mass (1000)
      integer isec,jsec,ial,jal(0:10),nf,irun,intint
      real*8  mint (1000)
! TK  tmass contains abundances (real) for each unit mass (i)
      real*8  tmass(10000)
      real*8  checksum2(50000)
      real*4, allocatable :: rnd(:,:)
! TK  cthr and cthr contain ions (counts) with different charges
!     tmax the maximum abundance factor over the whole spectrum
      real*8 xx(100),tmax,r,rms,norm,chrg,cthr,cthr2
      real*8 chrg1,chrg2,dum,checksum,cw,intensity
      real   snorm
      logical ex,sel,echo,exdat,mpop

! TK  fname=<qceims.res> or result file, xname contains the mass.agr plot file
      character*80 arg(10),line,fname,xname,formula,aaa,formula_former
      character*80 formula2
      character*2 a2
      character*5 aa,adum
      character*2 symbol(200)
      character*255 atmp

      ! TK number of peaks in in-silico spectra, needed for JCAMP-DX export
      integer numspec
      ! initialize the hash TABLE
      dictsize=1024
      array_size = 16
      call dict%init(dictsize)
      call mdict%init(dictsize)
      iexp =0
      mpop =.false.
      echo =.false.
      isec =0
      cthr =1.d-3
      cw   =0
      norm =1.0
      cthr2=1.01
      nagrfile=410
!SW symbol
      symbol(1)='H'
      symbol(6)='C'
      symbol(7)='N'
      symbol(8)='O'
      symbol(9)='F'
      symbol(14)='Si'
      symbol(15)='P'
      symbol(16)='S'
      symbol(17)='Cl'
      symbol(35)='Br'
      symbol(2)='He'
      symbol(10)='Ne'
      symbol(18)='Ar'
      symbol(36)='Kr'
      symbol(54)='Xe'
! edit this path name to some standard xmgrace plot file
! TK changed to direct working path

      xname='mass.agr'
      fname=''

      ! TK start loop reading arguments
      do i=1,9
         arg(i)=' '
         call getarg(i,arg(i))
      enddo
      ! TK end loop reading arguments


      ! TK start mega loop processing aruments
      do i=1,9

! TK comand line parameters
! TK -a no idea (related to ion charge count)
! TK -p print spectra "mass % intensity  counts   Int. exptl" to stdout;
!       with "Int. exptl" (experimental) taken from exp.dat but not all exp peaks are exported
!       if no theoretical counterpart exists
! TK -f filename or  -f <name_of_res_file>
! TK -t no idea
! TK -w no idea
! TK -s no idea

         if(index(arg(i),'-a').ne.0)  cthr=-1000.
         if(index(arg(i),'-p').ne.0)  echo=.true.
         if(index(arg(i),'-f').ne.0)  fname=arg(i+1)
         if(index(arg(i),'-t').ne.0) then

            ! TK call subroutine READL
            call readl(arg(i+1),xx,nn)
            cthr=xx(1)
         endif
         if(index(arg(i),'-w').ne.0) then
            call readl(arg(i+1),xx,nn)
            cw=xx(1)
         endif
         if(index(arg(i),'-s').ne.0) then
            call readl(arg(i+1),xx,nn)
            isec=int(xx(1))
         endif
      ! TK end mega loop processing aruments
      enddo

! fname contains the results from each calculation or the temporary result tmpqceims.res
! xname contains the xmgrace plot file

      if(fname.eq.'')then
        fname='qceims.res'
        inquire(file=fname,exist=ex)
        if(.not.ex) fname='tmpqceims.res'
      endif

      inquire(file=fname,exist=ex)
      if(.not.ex) stop 'res file does not exist'

      write(*,*) 'QCEIMS output reader PLOTMS'
      write(*,*) 'V 2.2, Nov 2013'
      write(*,*)
      write(*,*) 'xmgrace file body ',trim(xname)
      write(*,*) 'Reading ... ', trim(fname)

      if(cthr.ge.0)then
      write(*,'( &
      '' couting ions with charge from '',f4.1,'' to '',f4.1)')   cthr,cthr2
      else
      write(*,'( &
      '' counting all fragments with unity charge (frag. overview)'')')
      endif
      if(cw.gt.0)then
      write(*,'( &
     '' broadening the charges by an SD, wdth :'',f4.1)')cw*100.
      endif

      if(isec.ne.0)then
      write(*,*) &
      'Taking only secondary, tertiary ... fragmentations ',isec
      endif

      tmass = 0
      call tdict%init(dictsize)
      call mdict%init(dictsize)
      call random_seed()

! read file once
      i=1
      ndim=0

      ! TK contains qceims.res as standard option
      ! TK example output of qceims.res with 13 variables (in this case)
      ! 0.1986979   21 1 1    4     1  5   6  6   7  5   8  1

      !  other example from qceims.out
      !  trajectory          100           1
      !  mass                 formula              q pop   spin    q IPB  diss time (ps)
      !  M=121.12               H5C7O2   100   1   0.361   0.985   0.000    0.426
      !  M=73.19               H9C3SI1   100   1   0.639   0.015   1.000    0.426 ~

      !  is coded in qceims.res as
      !             trj # sub  #elem H  5   C  7   O  2
      !  0.0000088  100 1 1    3     1  5   6  7   8  2

      !             trj # sub  #elem H  9   C  3  Si  1
      !  0.9999912  100 2 1    3     1  9   6  3  14  1


      ! TK irun is optimized out during compile time
      ! example qceims.res
      ! 0.0000000   13 1 1    2     6  1   8  2

      ! chrg2 = 0
      ! irun = 13
      ! jsec = 1
      ! nf = 1
      ! k = 2 // kk = 3 (?)
      ! iat(kk) = 0
      ! nat(kk) = 0
      ! k = 2

      open(unit=1,file=fname)
      ! iat atomic number; nat number of atoms
 10   read(1,*,end=100)chrg2,irun,jsec,nf,k,(iat(kk),nat(kk),kk=1,k)
      !            0.9999912  100 2 1    3     1  9   6  3  14  1
      ! TK just for debug purposes
      write(*,*) chrg2,irun,jsec,nf,k,(iat(kk),nat(kk),kk=1,k)
      ! TK end debug purposes

      ! TK normal program
      if(isec.gt.0.and.isec.ne.jsec) goto 10
      if(chrg2.gt.cthr) then
         natot=sum(nat(1:k))
         if(natot.gt.ndim)ndim=natot
      endif
      i=i+1
      goto 10


100   continue
      n=i-1
      ! TK n contains number of fragments and ndim contains number of atoms(?) max
      write(*,*) n,' fragments with ',ndim,' atoms max.'
      close(1)

! initialize the random number array (efficiency)
      nrnd=50000
      allocate(rnd(ndim,nrnd))
      do i=1,nrnd
         do j=1,ndim
            call random_number(r)
            rnd(j,i)=r
         enddo
      enddo

! read it again
      write(*,*) 'Computing ...'

! TK contains qceims.res as standard option
      open(unit=1,file=fname)
      imin=100000
      i=1
      ial=0
      jal=0
      checksum =0
      checksum2=0
 11   read(1,*,end=101)chrg2,irun,jsec,nf,k,(iat(kk),nat(kk),kk=1,k)
      sel=.false.
      chrg=chrg2
      checksum2(irun)=checksum2(irun)+chrg2
      if(cthr.lt.0)chrg=1.0
      if(chrg.gt.cthr)then
         sel=.true.
         ial=ial+1
         jal(jsec)=jal(jsec)+1
      endif
      if(isec.gt.0.and.isec.ne.jsec) sel=.false.
      if(sel)then
        natot=sum(nat(1:k))
! all types of atoms
        kkkk=0
        length = 0

        do kk=1,k
! SW formula
          a2 = symbol(iat(kk))
          if(nat(kk).lt.100)then
          if(a2(2:2).eq.' ')then
          write(aa,'(a1,i2)')a2,nat(kk)
          nn = 3
          else
          write(aa,'(a2,i2)')a2,nat(kk)
          nn = 4
          endif
        endif

          if(nat(kk).lt.10)then
          if(a2(2:2).eq.' ')then
          write(aa,'(a1,i1)')a2,nat(kk)
          nn = 2
          else
          write(aa,'(a2,i1)')a2,nat(kk)
          nn = 3
          endif
        endif
          aaa(length+1:length+nn)=aa(1:nn)
          length = length + nn
! all atoms of this type, generate atomic number list, 3 C = [6 6 6]
          do kkk=1,nat(kk)
             kkkk=kkkk+1
             idum(kkkk)=iat(kk)
          enddo
        enddo
        checksum=checksum+chrg
        formula = trim(aaa)
        aaa = ''


! compute pattern, nsig signals at masses mass with int mint
        write(*,*)'calculating traj',irun,'step',jsec,'frag',nf
! natot total atom number
        call isotope(natot,idum,ndim,nrnd,rnd,mass,mint,nsig,dict, &
        mdict,key_list,formula)
        if(cw.gt.1.d-6)chrg=chrg+cw*chrg*snorm()
        call dict%show(key_list)
        s = SIZE(key_list)
! SW debug print
        !write(*,*)'check peaks'
        !do i = 1, s
        !  write(*,*) key_list(i)
        !end do
        !if (s .ne. nsig) write(*,*)s,'wrong with round up',nsig
        do k=1,s !s not nsig!!!!!!!!!!!
         charstrk=key_list(k)
         tmass(mass(k))=tmass(mass(k))+mint(k)*chrg

! SW accurate mass
        !write(*,*)'key value used in floatv function',charstrk
         tmp = floatv(charstrk,dict)
         !SW get frequency
         if (tdict%get(charstrk) .eq. ' ') THEN
           write(charstrv,'(F16.6)')tmp
           !write(*,*)'charstrv',charstrv
           call tdict%set(charstrk, charstrv)
         else
           ttmp = floatv(charstrk,tdict) ! total peak intensity
          ! write(*,*)'ttmp*chrg',ttmp
           ttmp = ttmp +tmp
           WRITE (charstrv,'(F16.6)')ttmp
           call tdict%set(charstrk, charstrv)
         endif
        enddo
        deallocate(key_list)
        i=i+1
      endif
      goto 11
101   continue
      write(*,*) n,' (charged) fragments done.'
      write(*,*)
      write(*,*) 'checksum of charge :',checksum
      write(*,*)

      k=0
      do i=1,50000
         if(abs(checksum2(i)).gt.1.d-6)k=k+1

         if(abs(checksum2(i)).gt.1.d-6.and. &
            abs(checksum2(i)-1.).gt.1.d-3)then
            write(*,*) 'checksum error for trj', i,' chrg=',checksum2(i)
         endif
      enddo
      write(*,*) k,' successfull runs.'

      write(*,*)
      write(*,*) 'ion sources:'
      write(*,'(''% inital run'',f6.1)') &
                100.*float(jal(1))/float(ial)
      do j=2,8
      write(*,'(''%        '',i1,''th'',f6.1)')j, &
                100.*float(jal(j))/float(ial)
      enddo
      write(*,*)

! read exp.
      imax=0
      imin=0
      inquire(file='exp.dat',exist=exdat)
      if(exdat) then
         write(*,*) 'Reading exp.dat ...'
         open(unit=4,file='exp.dat')

! TK (a) is edit descriptor for character strings
20       read(4,'(a)',end=200)line
         if(index(line,'##MAXY=').ne.0)then
            line(7:7)=' '
            call readl(line,xx,nn)
            norm=xx(nn)/100.0d0
         endif
         if(index(line,'##PEAK TABLE').ne.0)then
            kk=0
30          read(4,'(a)',end=200)line

            ! TK JCAMP DX for MS data has "##END="
            if(index(line,'##END').ne.0)goto 200
            do k=1,80
               if(line(k:k).eq.',')line(k:k)=' '
            enddo
            call readl(line,xx,nn)
            do k=1,nn/2
               kk=kk+1
               iexp(1,kk)=xx(2*k-1)
               iexp(2,kk)=xx(2*k)
               if(iexp(1,kk).gt.imax) imax=iexp(1,kk)
               if(iexp(1,kk).lt.imin) imin=iexp(1,kk)
            enddo
            goto 30
         endif
         goto 20
200      continue
         close(4)
      endif

      imin=max(imin,10)
      j=0
      do i=10,10000
         if(tmass(i).ne.0)j=i
      enddo
      imax=max(j,imax)

      tmax=maxval(tmass(10:10000))
      ! TK number of peaks computed in theoretical spectrum
      ! idint(tmax) is not related to the number of spectral peaks
      write(*,*)'Theoretical counts in 100 % signal:',idint(tmax)

      if(echo)then
      write(*,*)'mass % intensity  counts   Int. exptl'
      do i=10,10000
         if(tmass(i).ne.0)then
            dum=0
            do j=1,kk
               if(int(iexp(1,j)).eq.i)dum=iexp(2,j)
            enddo
            write(*,'(i8,F8.2,i8,F8.2)') &
            i,100.*tmass(i)/tmax,idint(tmass(i)),dum/norm
         endif
      enddo
      endif

! when the template exists
!TK Fortran runtime error: File already opened in another unit
!TK At line 285 of file src/plotms.f90 (unit = 3, file = '')

      inquire(file=xname,exist=ex)
      if(ex)then
      open(unit=2,file=xname)

! TK original mass.agr a preconfigured file
      open(unit=7,file='mass-result.agr')

! my xmgrace file mass.agr has nagrfile lines
      write(*,*) 'Writing mass-result.agr ...'
      do i=1,nagrfile
         read (2,'(a)')line
         if(index(line,'world 10').ne.0)then
            write(line,'(''@    world'',i3,'', -105,'' &
                           ,i3,'', 100'')')imin,imax+5
         endif
         if(index(line,'xaxis  tick major').ne.0)then
            line='@    xaxis  tick major 20'
            if(imax.gt.200) line='@    xaxis  tick major 50'
         endif
         if(index(line,'@    s1 symbol size').ne.0)then
            line='@    s1 symbol size 0.160000'
            if(imax.gt.200) line='@    s1 symbol size 0.100000'
         endif
         if(index(line,'@    s2 symbol size').ne.0)then
            line='@    s2 symbol size 0.160000'
            if(imax.gt.200) line='@    s2 symbol size 0.100000'
         endif
         write(7,'(a)')line
      enddo

! write only masses > 10
! TK tmass contains theoretical masses, unit 7 = mass-result.agr

      do i=10,10000
         if(tmass(i).ne.0)then
            write(7,*) i,100.*tmass(i)/tmax
         endif
      enddo
      write(7,*)'&'

! TK establish again if exdat (experimental JCAMPDX) exists
! TK @target G0.S1 means upper theory spectrum in current implementation
! TK @target G0.S2 means lower experimental spectrum
! TK iexp(1,i) - contains the exp masses and  iexp(2,i) contains the abundance
! TK kk contains number of experimental spectra

      if(exdat)then
      write(7,*)'@target G0.S2'
      write(7,*)'@type bar'
      do i=1,kk
         write(7,*) iexp(1,i),-iexp(2,i)/norm
      enddo
      write(7,*)'&'
      endif
      endif

! TK here we export the spectrum to JCAMP-DX, subroutine would be better for
! TK allowing later code merges or inclusion of new and missing elements such as
! TK currently unit mass only
    ! SW also generate accurate mass
    ! TK open file as JCAMP-DX MS result file, open new or replace
    write(*,*)'open file result.jdx'
    open(unit=11, file='result.jdx', STATUS="REPLACE")
    write(11,"(A)")'##TITLE=Theoretical in-silico spectrum (QCEIMS)'
    write(11,"(A)")'##JCAMP-DX=4.24'
    write(11,"(A)")'##DATA TYPE=MASS SPECTRUM'
    write(11,"(A)")'##XUNITS=M/Z'
    write(11,"(A)")'##YUNITS=RELATIVE INTENSITY'
! TK calculate number of in-silico spectra
! TK new: numspec = idint(tmax) ** not related to number ofr spectral peaks
     numspec = 0
     do i=10,10000
         if(tmass(i).ne.0)then
           intint = int(1000.*tmass(i)/tmax)
           if (intint >=1) then
            numspec = numspec + 1
         endif
       endif
     enddo

    write(11,'(A, I0)')'##NPOINTS=' ,numspec
    write(11,"(A)")'##PEAK TABLE=(XY..XY) 1'
    write(*,*)'write result.jdx....'
    s = 0
    do i=10,10000
         if(tmass(i).ne.0)then
           ! unit mass to exact mass
           ! write(11,"(I4, I8)") i, int(100.*tmass(i)/tmax)
           intint = int(1000.*tmass(i)/tmax)
           if (intint >=1) then
             s = s + 1
           write(11,"(I4, I8)") i, intint
         endif
         endif
    enddo
    write(11,"(A)")'##END='
    close(11)

    write(*,*)'open file accuratemass.jdx'
    open(unit=111, file='accuratemass.jdx',STATUS="REPLACE")
    write(111,"(A)")'##TITLE=Theoretical in-silico spectrum (QCEIMS)'
    write(111,"(A)")'##JCAMP-DX=4.24'
    write(111,"(A)")'##DATA TYPE=MASS SPECTRUM'
    write(111,"(A)")'##XUNITS=M/Z'
    write(111,"(A)")'##YUNITS=RELATIVE INTENSITY'
    call tdict%show(tkey_list)
    write(*,*)'***********'
    call mdict%show(mkey_list)
    s = SIZE(tkey_list)
    array_len = s


    write(*,*)'write accurate mass...'
!    write(*,*)'before sorting',tkey_list
    write(*,*)'sort key list'
    CALL qsort(tkey_list,array_len,array_size,cmp_function) ! sort mass in float order
    write(*,*)tkey_list
    kmax = 0.
! normalize m/z
    write(*,*)'normalizing m/z'
    do i = 1, s
      ctmp = tdict%get(tkey_list(i))
!      write(*,*)'ctmp',ctmp
      read (ctmp,*) ftmp
!      write(*,*)'ftmp',ftmp
      if (ftmp > kmax) THEN
        kmax = ftmp
      end if
    enddo
    si = 0
    do i = 1, s
      ctmp = tdict%get(tkey_list(i))
      read (ctmp,*) ftmp
      formula = mdict%get(tkey_list(i))
      intensity = 1000.*ftmp/kmax
      if (intensity >=1) then

        si = si + 1
      end if
    end do
    write(111,'(A, I0)')'##NPOINTS=' ,si
    write(111,"(A)")'##PEAK TABLE=(XY..XY) 1'


formula_former='' !formula from the last loop run
    do i = 1, s
      adum='' !isotpic peak sign
      ctmp = tdict%get(tkey_list(i))
      read (ctmp,*) ftmp
      formula = mdict%get(tkey_list(i))
      intensity = 1000.*ftmp/kmax
      if (intensity >=1) then
        formula2 = trim(formula)
        n =  len_trim( formula2)
        formula = formula2(1:n)
! SW try to remove null value in first entry
        pos = index( formula, achar(0) )
   if ( pos > 0 ) then
       formula = formula(1:pos)
     endif

if (formula==formula_former)adum='(iso)'
formula2  = '"'//trim(formula)//trim(adum)//'"'
write(111,9) ADJUSTL(trim(tkey_list(i))),1000.*ftmp/kmax,formula2
    formula_former = formula
    end if

    end do

    write(111,"(A)")'##END='

    close(111)

    write(*,*)'writing file formula.txt'
    open(unit=1111, file='formula.txt',STATUS="REPLACE", encoding='UTF-8')

    write(1111,'(A, I0)')'##NPOINTS=' ,s
    !write(*,*)'maximum m/z', kmax
    write(1111,'(a,2x,a,2x,a)')'formula','m/z','intensity'
    formula_former=''
    do i = 1, s
      adum=''
      ctmp = tdict%get(tkey_list(i))
      read (ctmp,*) ftmp
      formula = ADJUSTR(trim(mdict%get(tkey_list(i))))
      formula2 = trim(formula)
      n =  len_trim( formula2)
      formula = formula2(1:n)
      if (formula==formula_former)adum='(iso)'
      write(1111,99) trim(formula)//adum,tkey_list(i), 1000.*ftmp/kmax
      formula_former = formula
    end do
    write(1111,"(A)")'##END='
    close(1111)


9   format(a,2x, F7.2, 2x, a20)
8   format(a20,2x,a,2x, I8)
99   format(a20,2x, a, 2x, F7.2)



! compute deviation exp-theor.
! TK here we potentially have to use Mass Spec related terms, such as dot product
! From Stein and Scott

      if(exdat)then
      rms=0
      nn=0
      kkk=0
      kkkk=0
      do i=1,kk
         k=iexp(1,i)
         if(iexp(2,i)/norm.gt.5.0)then
            r=100.*tmass(k)/tmax-iexp(2,i)/norm
            kkk=kkk+1
            if(100.*tmass(k)/tmax.gt.2.5)kkkk=kkkk+1
            rms=rms+abs(r)
         endif
      enddo
      write(*,*)'MAD(exptl./theor.) = ',rms/kkk
      write(*,*)'# exptl. > 5 %     = ',kkk
      write(*,*)'% correctly found  = ',100*float(kkkk)/float(kkk)
      endif

      close(1)
      close(2)
      close(3)
      close(7)
      end
!