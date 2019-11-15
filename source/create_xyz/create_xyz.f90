#define PRC       8

#define ZERO      0.d0
#define ONE       1.d0
#define TWO       2.d0
#define THREE     3.d0
#define FOUR      4.d0
#define FIVE      5.d0
#define SIX       6.d0
#define EIGHT     8.d0
#define NINE      9.d0
#define TEN      10.d0
#define ELEVEN   11.d0
#define TWELVE   12.d0
#define SIXTEEN  16.d0
#define EIGHTEEN  18.d0
#define TWENTYFOUR 24.d0
#define TWENTYSEVEN 27.d0
#define THIRTY   30.d0
#define THIRTYSIX 36.d0
#define FIFTY     50.d0
#define HALF      0.5d0
#define FOURTH    0.25d0

#define MINDENS   1.d-8

#define CHECKQUAT

 program create_xyz
 
  implicit none
  
  integer :: nx=64
  integer :: ny=64
  integer :: nz=64
  
  integer :: ntot,i,j,k,natms,rotate,ongrid,inode,nmax,nxy,ii,jj,kk,l,ll
  integer, parameter :: imcon=2
  integer, parameter :: idrank=0
  
  integer, parameter :: ioxyz=16
  character(len=*), parameter :: filexyz='input.xyz'
  
  real(kind=8), parameter :: & 
   Pi=3.141592653589793238462643383279502884d0
   
  integer, parameter :: maxlen=85
  real(kind=8) :: rcut,dtemp_1,dtemp_2,dtemp_3,cell(9),sqrcut,rdist,sigm
  real(kind=8), allocatable, dimension(:) :: xxx,yyy,zzz,vxx,vyy,vzz, &
   q0,q1,q2,q3
  
  character(len=maxlen) :: arg,redstring,directive
  logical :: safe,lcheck,ltest
  integer :: idum,iseed,inumchar
  
  real(kind=PRC), parameter :: toll=real(1.d-4,kind=PRC)
  logical :: ltestrot=.false.
  real(kind=PRC), dimension(9) :: rot,newrot
  real(kind=PRC) :: matrixmio(3,3),dmio(3)
  real(kind=PRC), dimension(3) :: xsubm,ysubm,zsubm
  
  integer, allocatable, dimension(:) :: idnodes
  
  real(kind=PRC) :: qs(0:3),phi,theta,psi
  
  
   
  if(iargc()/=9)then
    write(6,*) 'error!'
    write(6,*) 'the command line should be'
    write(6,*) '[executable] [natms] [nx] [ny] [nz] [seed] [rcut] [rotate]'
    write(6,*) 'natms = integer indicating the number of particle'
    write(6,*) 'nx    = integer indicating the box lenght along y'
    write(6,*) 'ny    = integer indicating the box lenght along y'
    write(6,*) 'nz    = integer indicating the box lenght along z'
    write(6,*) 'seed  = seed of the random generator             '
    write(6,*) 'rcut  = float indicating the smallest distance   '
    write(6,*) 'sigm  = float indicating the velocity dispersion '
    write(6,*) 'ongrid= 0 for no, 1 for yes                      '
    write(6,*) 'rotate= 0 for no, 1 for yes                      '
    write(6,*) 'STOP!'
    stop
  endif
  
  
  
  do i = 1, iargc()
    call getarg(i, arg)
    if(i==1)then
      call copystring(arg,directive,maxlen)
      natms=intstr(directive,maxlen,inumchar)
      write(6,*) 'natms = ',natms
    elseif(i==2)then
      call copystring(arg,directive,maxlen)
      nx=intstr(directive,maxlen,inumchar)
      write(6,*) 'nx    = ',nx
    elseif(i==3)then
      call copystring(arg,directive,maxlen)
      ny=intstr(directive,maxlen,inumchar)
      write(6,*) 'ny    = ',ny
    elseif(i==4)then
      call copystring(arg,directive,maxlen)
      nz=intstr(directive,maxlen,inumchar)
      write(6,*) 'nz    = ',nz
    elseif(i==5)then
      call copystring(arg,directive,maxlen)
      iseed=intstr(directive,maxlen,inumchar)
      write(6,*) 'seed    = ',iseed
    elseif(i==6)then
      call copystring(arg,directive,maxlen)
      rcut=dblstr(directive,maxlen,inumchar)
      write(6,*) 'rcut  = ',rcut
    elseif(i==7)then
      call copystring(arg,directive,maxlen)
      sigm=dblstr(directive,maxlen,inumchar)
      write(6,*) 'sigm  = ',sigm
    elseif(i==8)then
      call copystring(arg,directive,maxlen)
      ongrid=intstr(directive,maxlen,inumchar)
      write(6,*) 'ongrid= ',ongrid
      if(ongrid/=0 .and. ongrid/=1)then
        write(6,*)'ERROR: ongrid should be 0 or 1'
        stop
      endif
    elseif(i==9)then
      call copystring(arg,directive,maxlen)
      rotate=intstr(directive,maxlen,inumchar)
      write(6,*) 'rotate= ',rotate
      if(rotate/=0 .and. rotate/=1)then
        write(6,*)'ERROR: rotate should be 0 or 1'
        stop
      endif
    endif
  enddo
  
  call init_random_seed(iseed)
  
  allocate(xxx(natms),yyy(natms),zzz(natms),vxx(natms),vyy(natms), &
   vzz(natms))
   
  if(rotate==1)then
    allocate(q0(natms),q1(natms),q2(natms),q3(natms))
  endif
   
  cell(:)=0.d0
  cell(1)=dble(nx)+1.d0
  cell(5)=dble(ny)+1.d0
  cell(9)=dble(nz)+1.d0
  

        
  
  if(ongrid==1)then
    
    sqrcut=rcut**2.d0
    
    nmax=nx*ny*nz
    allocate(idnodes(nmax))
    forall(i=1:nmax)
      idnodes(i)=i
    end forall
    
    do i=1,natms
      lcheck=.true.
      do while(lcheck)
        call random_number(dtemp_1)
        
        ll=floor(dtemp_1*dble(nmax))+1
        
        l=idnodes(ll)
        idnodes(ll)=idnodes(nmax)
        nmax=nmax-1
        if(nmax<1)then
          write(6,*)'numeri da estrarre finiti'
          stop
        endif

        kk=(l-1)/(nx*ny)+1
        jj=(l-1-(nx*ny)*(kk-1))/nx+1
        ii=(l-(nx*ny)*(kk-1)-nx*(jj-1))
        
        dtemp_1=dble(ii)
        dtemp_2=dble(jj)
        dtemp_3=dble(kk)
        
        if(dtemp_1<1 .or. dtemp_1>nx)then
          write(6,*)'cazzi i ',dtemp_1,inode
        endif
        
        if(dtemp_2<1 .or. dtemp_2>ny)then
          write(6,*)'cazzi j ',dtemp_2,inode
        endif
        
        if(dtemp_3<1 .or. dtemp_3>nz)then
          write(6,*)'cazzi k ',dtemp_3,inode
        endif
        
        do j=1,i-1
          vxx(j)=(xxx(j)-dtemp_1)
          vyy(j)=(yyy(j)-dtemp_2)
          vzz(j)=(zzz(j)-dtemp_3)
        enddo
        if(i-1>0)call images_mio(imcon,i-1,cell,vxx,vyy,vzz)
        ltest=.true.
        do j=1,i-1
          rdist=vxx(j)**2.d0+vyy(j)**2.d0+vzz(j)**2.d0
          if(rdist<=sqrcut)then
            ltest=.false.
            exit
          endif
        enddo
        
        if(ltest)then
          lcheck=.false.
          xxx(i)=dtemp_1
          yyy(i)=dtemp_2
          zzz(i)=dtemp_3
          if(mod(i,10)==0)write(6,*)'iatm = ',i
        endif
      enddo
    enddo
  
  else
  
    sqrcut=rcut**2.d0
    
    do i=1,natms
      lcheck=.true.
      do while(lcheck)
        call random_number(dtemp_1)
        call random_number(dtemp_2)
        call random_number(dtemp_3)
        dtemp_1=dtemp_1*(dble(nx)+0.5d0)+0.5d0
        dtemp_2=dtemp_2*(dble(ny)+0.5d0)+0.5d0
        dtemp_3=dtemp_3*(dble(nz)+0.5d0)+0.5d0
        
        do j=1,i-1
          vxx(j)=(xxx(j)-dtemp_1)
          vyy(j)=(yyy(j)-dtemp_2)
          vzz(j)=(zzz(j)-dtemp_3)
        enddo
        if(i-1>0)call images_mio(imcon,i-1,cell,vxx,vyy,vzz)
        ltest=.true.
        do j=1,i-1
          rdist=vxx(j)**2.d0+vyy(j)**2.d0+vzz(j)**2.d0
          if(rdist<=sqrcut)then
            ltest=.false.
            exit
          endif
        enddo
        if(ltest)then
          lcheck=.false.
          xxx(i)=dtemp_1
          yyy(i)=dtemp_2
          zzz(i)=dtemp_3
          if(mod(i,10)==0)write(6,*)'iatm = ',i
        endif
      enddo
    enddo
  
  endif
  
  do i=1,natms
    vxx(i)=sigm*gauss()
    vyy(i)=sigm*gauss()
    vzz(i)=sigm*gauss()
  enddo
  
  if(rotate==1)then
  !transform the rotational matrix in quaternions using an uniform
    !distribution
    do i=1,natms
      xsubm(1:3)=ZERO
      ysubm(1:3)=ZERO
      zsubm(1:3)=ZERO
      xsubm(1)=ONE
      ysubm(2)=ONE
      zsubm(3)=ONE
      call random_number(dmio(1))
      call random_number(dmio(2))
      call random_number(dmio(3))
      dmio(1:3)=dmio(1:3)*TWO*Pi
      do j=1,3
        call matrix_rotaxis(j,dmio(j),matrixmio)
        call rotate_vect(1,3,xsubm,ysubm,zsubm,(/ZERO,ZERO,ZERO/), &
         matrixmio)
      enddo
      rot(1)=xsubm(1)
      rot(2)=xsubm(2)
      rot(3)=xsubm(3)
      rot(4)=ysubm(1)
      rot(5)=ysubm(2)
      rot(6)=ysubm(3)
      rot(7)=zsubm(1)
      rot(8)=zsubm(2)
      rot(9)=zsubm(3)
      call rotmat_2_quat(rot,q0(i),q1(i),q2(i),q3(i))
#ifdef CHECKQUAT
      call quat_2_rotmat(q0(i),q1(i),q2(i),q3(i),newrot)
      do j=1,9
        if(abs(newrot(j)-rot(j))>toll)then
          ltestrot=.true.
        endif
      enddo
#endif
    enddo
  endif
  
  open(unit=ioxyz,file=filexyz,status='replace',action='write')
  
  write(ioxyz,'(i8)')natms
  if(rotate==0)then
  
    write(ioxyz,'(a)')' read list vx vy vz '
    
    do i=1,natms
      write(ioxyz,'(a8,6f16.8)')'C       ',xxx(i),yyy(i),zzz(i), &
       vxx(i),vyy(i),vzz(i) 
    enddo
  
  elseif(rotate==1)then
  
    write(ioxyz,'(a)')' read list vx vy vz phi psi theta'
    
    do i=1,natms
      qs(0)=q0(i)
      qs(1)=q1(i)
      qs(2)=q2(i)
      qs(3)=q3(i)
      call q2eul(qs,phi,psi,theta)
      write(ioxyz,'(a8,9f16.8)')'C       ',xxx(i),yyy(i),zzz(i), &
       vxx(i),vyy(i),vzz(i),phi,psi,theta
    enddo
  
  endif
  
  close(ioxyz)
  
  stop
  
 contains
 
 
 subroutine init_random_seed(myseed)
 
!***********************************************************************
!     
!     LBsoft subroutine for initialising the random generator
!     by the seed given in input file or by a random seed
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  integer,intent(in),optional :: myseed
  integer :: i, n, clock
  
  integer, allocatable :: seed(:)
          
  call random_seed(size = n)
  
  allocate(seed(n))
  
  if(present(myseed))then
!   If the seed is given in input
    seed = myseed*(idrank+1) + 37 * (/ (i - 1, i = 1, n) /)
    
  else
!   If the seed is not given in input it is generated by the clock
    call system_clock(count=clock)
         
    seed = clock*(idrank+1) + 37 * (/ (i - 1, i = 1, n) /)
    
  endif
  
  call random_seed(put = seed)
       
  deallocate(seed)
  
  return
 
 end subroutine init_random_seed
 
 function gauss()
 
!***********************************************************************
!     
!     LBsoft subroutine for generating random number normally
!     distributed by the Box-Muller transformation
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  real(kind=8) :: gauss
  real(kind=8) :: dtemp1,dtemp2
  logical :: lredo
  
  call random_number(dtemp1)
  call random_number(dtemp2)
  
  lredo=.true.
  
! the number is extract again if it is nan
  do while(lredo)
    lredo=.false.
!   Box-Muller transformation
    gauss=dsqrt(- 2.d0 *dlog(dtemp1))*dcos(2*pi*dtemp2)
    if(isnan(dcos(gauss)))lredo=.true.
  enddo
  
 end function gauss
 
 subroutine copystring(oldstring,newstring,lenstring)
 
!***********************************************************************
!     
!     JETSPIN subroutine to copy one character string into another
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(in) :: oldstring
  character(len=*), intent(out) :: newstring
  integer, intent(in) :: lenstring
  
  integer :: i
  
  do i=1,lenstring
    newstring(i:i)=oldstring(i:i)
  enddo
  
  return
  
 end subroutine copystring
 
 function intstr(string,lenstring,laststring)
 
!***********************************************************************
!     
!     JETSPIN subroutine for extracting integers from a character 
!     string
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(inout) :: string
  integer, intent(in) :: lenstring
  integer, intent(out) :: laststring
  
  integer :: intstr
  
  integer :: j,isn
  character*1, parameter, dimension(0:9) :: & 
   n=(/'0','1','2','3','4','5','6','7','8','9'/)
  logical :: flag,lcount,final
  character*1 :: ksn
  character*1, dimension(lenstring) :: word
  
  do j=1,lenstring
    word(j)=string(j:j)
  enddo
  
  isn=1
  laststring=0
  ksn='+'
  intstr=0
  flag=.false.
  final=.false.
  lcount=.false.
  
  
  do while(laststring<lenstring.and.(.not.final))
    
    laststring=laststring+1
    flag=.false.
    
    do j=0,9
      
      if(n(j)==word(laststring))then
        
        intstr=10*intstr+j
        lcount=.true.
        flag=.true.
        
      endif
    
    enddo
    
    if(lcount.and.(.not.flag))final=.true.
    if(flag .and. ksn=='-')isn=-1
    ksn=word(laststring)
    
  enddo

  intstr=isn*intstr

  do j=laststring,lenstring
    word(j-laststring+1)=word(j)
  enddo
  do j=lenstring-laststring+2,lenstring
    word(j)=' '
  enddo
  
  do j=1,lenstring
    string(j:j)=word(j)
  enddo
  
  return
  
  end function intstr

  function dblstr(string,lenstring,laststring)
  
!***********************************************************************
!     
!     JETSPIN subroutine for extracting double precisions from a  
!     character string
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(inout) :: string
  integer, intent(in) :: lenstring
  integer, intent(out) :: laststring
  
  double precision :: dblstr
  
  logical :: flag,ldot,start,final
  integer :: iexp,idum,i,j,fail
  double precision :: sn,ten,one

  character*1, parameter, dimension(0:9) :: & 
   n=(/'0','1','2','3','4','5','6','7','8','9'/)
  character*1, parameter :: dot='.'
  character*1, parameter :: d='d'
  character*1, parameter :: e='e'
  
  character*1 :: ksn
  character*1, dimension(lenstring) :: word
  character(len=lenstring) :: work
  
  do j=1,lenstring
    word(j)=string(j:j)
  enddo

  laststring=0
  sn=1.d0
  ksn='+'
  ten=10.d0
  one=1.d0
 
  dblstr=0.d0
  iexp=0
  idum=0
  start=.false.
  ldot=.false.
  final=.false.
  
  do while(laststring<lenstring .and. (.not.final))
    
    laststring=laststring+1
    flag=.false.
    
    do j=0,9
      
      if(n(j)==word(laststring))then
        
        dblstr=ten*dblstr+one*dble(j)
        flag=.true.
        start=.true.
            
      endif
          
    enddo
        
    if(dot==word(laststring))then
          
      flag=.true.
      ten=1.d0
      ldot=.true.
      start=.true.
          
    endif

    if(flag .and. ksn=='-') sn=-1.d0
    if(ldot)one=one/10.d0
    ksn=word(laststring)
    if(ksn=="D")ksn="d"
    if(ksn=="E")ksn="e"
    
    if(start)then
      if(d==ksn .or. e==ksn)then
        do i=1,lenstring-laststring
          work(i:i)=word(i+laststring)
        enddo
        iexp=intstr(work,lenstring-laststring,idum)
        final=.true.
      endif
      if(.not.flag)final=.true.        
    endif
  enddo
  
  dblstr=sn*dblstr*(10.d0**iexp)
  laststring=laststring+idum
  
  do j=laststring,lenstring
    word(j-laststring+1)=word(j)
  enddo
  do j=lenstring-laststring+2,lenstring
    word(j)=' '
  enddo
  
  do j=1,lenstring
    string(j:j)=word(j)
  enddo
  
  return
  
 end function dblstr

 subroutine strip(string,lenstring)

  implicit none
  
  character(len=*) :: string
  integer, intent(in) :: lenstring
  
  integer :: i,j
  character*1, dimension(lenstring) :: word
  
  do j=1,lenstring
    word(j)=string(j:j)
  enddo 
  
  do i=1,lenstring
    
    if(word(1)==' ')then
      
      do j=1,lenstring-1
        
        word(j)=word(j+1)
        
      enddo
      
      word(lenstring)=' '
      
    endif
    
  enddo
  
  do j=1,lenstring
    string(j:j)=word(j)
  enddo
  
  return
  
 end subroutine strip
 
 subroutine images_mio(imconsub,natmssub,cellsub,xxx,yyy,zzz)
      
!***********************************************************************
!     
!     subroutine for calculating the minimum image
!     of atom pairs within a specified MD cell
!     
!     copyright - daresbury laboratory 1992
!     author    - w. smith march 1992.
!     T3D optimised version. t.forester july 1994
!     serial version M. Lauricella 2016
!     
!     for
!     imconsub=0 no boundary conditions apply
!     imconsub=1 standard cubic boundaries apply
!     imconsub=2 orthorhombic boundaries apply
!     imconsub=3 parallelepiped boundaries apply
!     imconsub=4 truncated octahedron boundaries apply
!     imconsub=5 rhombic dodecahedron boundaries apply
!     imconsub=6 x-y parallelogram boundary conditions : no periodicity in z
!     imconsub=7 hexagonal prism boundaries apply
!     
!     note: in all cases the centre of the cell is at (0,0,0)
!     warning - replicated data version: does not re-merge 
!     coordinate arrays
!     
!***********************************************************************

      
  implicit none
  
  integer, intent(in) :: imconsub,natmssub
  double precision, dimension(9), intent(in) :: cellsub
  integer :: i
  double precision, allocatable, dimension(:), intent(inout) ::  xxx,yyy,zzz
  double precision ::aaa,bbb,ccc,det,ssx
  double precision :: ssy,ssz,ddd,xss,yss,zss
  
  double precision :: rcellsub(9)
  double precision, parameter :: rt2=1.41421356623d0
  double precision, parameter :: rt3=1.7320508075d0
  
      
  if(imconsub.eq.1)then

!   standard cubic boundary conditions
        
        
    aaa=1.d0/cellsub(1)
    
    do i=1,natmssub
      xxx(i)=xxx(i)-cellsub(1)*nint(aaa*xxx(i))
      yyy(i)=yyy(i)-cellsub(1)*nint(aaa*yyy(i))
      zzz(i)=zzz(i)-cellsub(1)*nint(aaa*zzz(i))
    enddo
        
  else if(imconsub.eq.2)then

!   rectangular (slab) boundary conditions
        
    aaa=1.d0/cellsub(1)
    bbb=1.d0/cellsub(5)
    ccc=1.d0/cellsub(9)
        
    do i=1,natmssub
      
      xxx(i)=xxx(i)-cellsub(1)*nint(aaa*xxx(i))
      yyy(i)=yyy(i)-cellsub(5)*nint(bbb*yyy(i))
      zzz(i)=zzz(i)-cellsub(9)*nint(ccc*zzz(i))
          
    enddo
      
  else if(imconsub.eq.3)then

!   parallelepiped boundary conditions
        
    call invert(cellsub,rcellsub,det)
        
    do i=1,natmssub
          
      ssx=(rcellsub(1)*xxx(i)+rcellsub(4)*yyy(i)+rcellsub(7)*zzz(i))
      ssy=(rcellsub(2)*xxx(i)+rcellsub(5)*yyy(i)+rcellsub(8)*zzz(i))
      ssz=(rcellsub(3)*xxx(i)+rcellsub(6)*yyy(i)+rcellsub(9)*zzz(i))
          
      xss=ssx-nint(ssx)
      yss=ssy-nint(ssy)
      zss=ssz-nint(ssz)
          
      xxx(i)=(cellsub(1)*xss+cellsub(4)*yss+cellsub(7)*zss)
      yyy(i)=(cellsub(2)*xss+cellsub(5)*yss+cellsub(8)*zss)
      zzz(i)=(cellsub(3)*xss+cellsub(6)*yss+cellsub(9)*zss)
          
    enddo
        
  else if(imconsub.eq.4)then

!   truncated octahedral boundary conditions
        
    if(.not.(abs(cellsub(1)-cellsub(5)).lt.1.d-6.and. &
     abs(cellsub(5)-cellsub(9)).lt.1.d-6))then
      write(6,'(2a)')'.not.(abs(cellsub(1)-cellsub(5)).lt.1.d-6.and.', &
       'abs(cellsub(5)-cellsub(9)).lt.1.d-6)'
      stop
    endif
        
    aaa=1.d0/cellsub(1)
        
    do i=1,natmssub
          
      xxx(i)=xxx(i)-cellsub(1)*nint(aaa*xxx(i))
      yyy(i)=yyy(i)-cellsub(1)*nint(aaa*yyy(i))
      zzz(i)=zzz(i)-cellsub(1)*nint(aaa*zzz(i))
          
      if((abs(xxx(i))+abs(yyy(i))+abs(zzz(i))).ge.(0.75d0*cellsub(1)))then
            
        xxx(i)=xxx(i)-0.5d0*sign(cellsub(1),xxx(i))
        yyy(i)=yyy(i)-0.5d0*sign(cellsub(1),yyy(i))
        zzz(i)=zzz(i)-0.5d0*sign(cellsub(1),zzz(i))
            
      endif
          
    enddo
        
  else if(imconsub.eq.5)then

!   rhombic dodecahedral boundary conditions
        
    if(.not.(abs(cellsub(1)-cellsub(5)).lt.1.d-6.and. &
     abs(cellsub(9)-cellsub(1)*rt2).lt.1.d-6))then
      write(6,'(2a)')'.not.(abs(cellsub(1)-cellsub(5)).lt.1.d-6.and.', &
       'abs(cellsub(9)-cellsub(1)*rt2).lt.1.d-6)'
      stop
    endif
        
    aaa=1.d0/cellsub(1)
    bbb=1.d0/cellsub(9)
        
    do i=1,natmssub
          
      xxx(i)=xxx(i)-cellsub(1)*nint(aaa*xxx(i))
      yyy(i)=yyy(i)-cellsub(1)*nint(aaa*yyy(i))
      zzz(i)=zzz(i)-cellsub(9)*nint(bbb*zzz(i))
          
      if((abs(xxx(i))+abs(yyy(i))+abs(rt2*zzz(i))).ge.cellsub(1))then
            
        xxx(i)=xxx(i)-0.5d0*sign(cellsub(1),xxx(i))
        yyy(i)=yyy(i)-0.5d0*sign(cellsub(1),yyy(i))
        zzz(i)=zzz(i)-0.5d0*sign(cellsub(9),zzz(i))
            
      endif
          
    enddo
        
  else if(imconsub.eq.6) then

!   x-y boundary conditions 

    det = cellsub(1)*cellsub(5) - cellsub(2)*cellsub(4)

    if(abs(det).lt.1.d-6)then
      write(6,'(a)')'abs(det).lt.1.d-6'
      stop
    endif
        
    det = 1.d0/det

    rcellsub(1) =  det*cellsub(5)
    rcellsub(2) = -det*cellsub(2)
    rcellsub(4) = -det*cellsub(4)
    rcellsub(5) =  det*cellsub(1)
        
    do i=1,natmssub

      ssx = rcellsub(1)*xxx(i) + rcellsub(4)*yyy(i)
      ssy = rcellsub(2)*xxx(i) + rcellsub(5)*yyy(i)

      xss = ssx - nint(ssx)
      yss = ssy - nint(ssy)

      xxx(i)=cellsub(1)*xss + cellsub(4)*yss
      yyy(i)=cellsub(2)*xss + cellsub(5)*yss

    enddo

  else if(imconsub.eq.7) then

!   hexagonal prism boundary conditions
        
    if(abs(cellsub(1)-rt3*cellsub(5)).ge.1.d-6)then
      write(6,'(a)')'abs(cellsub(1)-rt3*cellsub(5)).ge.1.d-6'
      stop
    endif
        
    aaa=cellsub(1)/(rt3*2.d0)
    bbb=cellsub(1)/rt3
    ccc=rt3/cellsub(1)
    ddd=1.d0/cellsub(9)
        
    do i=1,natms
          
      yyy(i)=yyy(i)-bbb*nint(ccc*yyy(i))
      zzz(i)=zzz(i)-cellsub(9)*nint(ddd*zzz(i))
          
      if((abs(yyy(i))+abs(rt3*xxx(i))).ge.bbb)then
            
        xxx(i)=xxx(i)-rt3*sign(aaa,xxx(i))
        yyy(i)=yyy(i)-sign(aaa,yyy(i))
            
      endif
          
    enddo
        
  endif
      
  return

 end subroutine images_mio
 
 subroutine invert(a,b,d)

!***********************************************************************
!     
!     dl_poly subroutine to invert a 3 * 3 matrix using cofactors
!     
!     copyright - daresbury laboratory 1992
!     author    - w. smith       april 1992
!     
!***********************************************************************

  implicit none

  double precision :: a(9),b(9),d,r


! calculate adjoint matrix
  b(1)=a(5)*a(9)-a(6)*a(8)
  b(2)=a(3)*a(8)-a(2)*a(9)
  b(3)=a(2)*a(6)-a(3)*a(5)
  b(4)=a(6)*a(7)-a(4)*a(9)
  b(5)=a(1)*a(9)-a(3)*a(7)
  b(6)=a(3)*a(4)-a(1)*a(6)
  b(7)=a(4)*a(8)-a(5)*a(7)
  b(8)=a(2)*a(7)-a(1)*a(8)
  b(9)=a(1)*a(5)-a(2)*a(4)

! calculate determinant
  d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
  r=0.d0
  if(abs(d).gt.0.d0)r=1.d0/d

! complete inverse matrix
  b(1)=r*b(1)
  b(2)=r*b(2)
  b(3)=r*b(3)
  b(4)=r*b(4)
  b(5)=r*b(5)
  b(6)=r*b(6)
  b(7)=r*b(7)
  b(8)=r*b(8)
  b(9)=r*b(9)

  return
  
 end subroutine invert
 
 subroutine rotmat_2_quat(rot,q0s,q1s,q2s,q3s)
 
!***********************************************************************
!     
!     LBsoft subroutine for determining the quaternion associated
!     to a rotational matrix within a given tollerance value
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in), dimension(9) :: rot
  real(kind=PRC), intent(out) :: q0s,q1s,q2s,q3s
  
  real(kind=PRC) :: aq,bq,cq,dq,eq,fq,gq,hq,rnorm
  
  real(kind=PRC), parameter :: toll=real(1.d-4,kind=PRC)
 
 !determine quaternions from rotational matrix
  aq=rot(1)+rot(5)
  bq=rot(2)-rot(4)
  cq=rot(6)-rot(8)
  dq=rot(2)+rot(4)
  eq=rot(3)+rot(7)
  fq=rot(6)+rot(8)
  gq=rot(3)-rot(7)
  hq=rot(1)-rot(5)
          
  q0s=HALF*sqrt(aq+sqrt(aq*aq+bq*bq))
          
  if(q0s.gt.toll)then
    q1s=-FOURTH*cq/q0s
    q2s=FOURTH*gq/q0s
    q3s=-FOURTH*bq/q0s
  else
    q1s=HALF*sqrt(hq+sqrt(hq*hq+dq*dq))
    if(q1s.gt.toll)then
      q2s=FOURTH*dq/q1s
      q3s=FOURTH*eq/q1s
    else
      q2s=HALF*sqrt(-hq+sqrt(hq*hq+dq*dq))
      if(q2s.gt.toll)then
        q3s=FOURTH*fq/q2s
      else
        q3s=ONE
      endif
    endif
  endif
  
 !normalize quaternions
  rnorm=ONE/sqrt(q0s**TWO+q1s**TWO+q2s**TWO+q3s**TWO)
  q0s=rnorm*q0s
  q1s=rnorm*q1s
  q2s=rnorm*q2s
  q3s=rnorm*q3s
  
  return
  
 end subroutine rotmat_2_quat
 
 subroutine quat_2_rotmat(q0,q1,q2,q3,rot)
 
!***********************************************************************
!     
!     LBsoft subroutine for determining the rotational matrix from
!     quaternion using x convention for euler angles
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
      
  real(kind=PRC), intent(in) :: q0,q1,q2,q3
  real(kind=PRC), intent(out), dimension(9) :: rot
      
  rot(1)=q0**TWO+q1**TWO-q2**TWO-q3**TWO
  rot(2)=TWO*(q1*q2-q0*q3)
  rot(3)=TWO*(q1*q3+q0*q2)
  rot(4)=TWO*(q1*q2+q0*q3)
  rot(5)=q0**TWO-q1**TWO+q2**TWO-q3**TWO
  rot(6)=TWO*(q2*q3-q0*q1)
  rot(7)=TWO*(q1*q3-q0*q2)
  rot(8)=TWO*(q2*q3+q0*q1)
  rot(9)=q0**TWO-q1**TWO-q2**TWO+q3**TWO
  
  return
  
 end subroutine quat_2_rotmat
 
 subroutine rotate_vect(isub,lsub,xsubmy,ysubmy,zsubmy,center,matrix)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying a rotation of an array
!     depending on a rotation matrix
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
  
  implicit none
  
  integer, intent(in) :: isub,lsub
  real(kind=PRC), dimension(*), intent(inout) :: xsubmy,ysubmy,zsubmy
  real(kind=PRC), dimension(3), intent(in) :: center
  real(kind=PRC), dimension(3,3), intent(in) :: matrix
  
  integer :: i
  real(kind=PRC), dimension(3) :: dvec1,dvec2
  
  do i=isub,lsub
    dvec1(1)=xsubmy(i)-center(1)
    dvec1(2)=ysubmy(i)-center(2)
    dvec1(3)=zsubmy(i)-center(3)
    call rotation(dvec1,matrix,dvec2)
    xsubmy(i)=dvec2(1)+center(1)
    ysubmy(i)=dvec2(2)+center(2)
    zsubmy(i)=dvec2(3)+center(3)
  enddo
  
  return
  
 end subroutine rotate_vect
 
 subroutine rotation(vec1,matrix,vec2)
 
!***********************************************************************
!     
!     LBsoft subroutine for applying a rotation of a single vector
!     depending on a rotation matrix
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: vec1(3),matrix(3,3)
  real(kind=PRC), dimension(3), intent(out) :: vec2
  
  integer :: i,j
  
  vec2(:)=ZERO
  do i=1,3
    do j=1,3
      vec2(i)=matrix(i,j)*vec1(j)+vec2(i)
    enddo
  enddo
  
  return
  
 end subroutine rotation
 
 subroutine matrix_rotaxis(iaxis,theta,matrix)
 
!***********************************************************************
!     
!     LBsoft subroutine for building the rotation matrix
!     of a given angle and axis
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2018
!     
!***********************************************************************
 
  implicit none
  
  integer, intent(in) :: iaxis
  real(kind=PRC), intent(in) :: theta
  real(kind=PRC), dimension(3,3), intent(out) :: matrix
  
  real(kind=PRC) :: ct,st
  
  ct=cos(theta)
  st=sin(theta)
  
  matrix(:,:)= ZERO
  
  select case(iaxis)
    case(1)
      matrix(1,1)= ONE
      matrix(2,2)=ct
      matrix(3,2)=-st
      matrix(2,3)=st
      matrix(3,3)=ct
    case(2)
      matrix(1,1)=ct
      matrix(3,1)=st
      matrix(2,2)= ONE
      matrix(1,3)=-st
      matrix(3,3)=ct
    case(3)
      matrix(1,1)=ct
      matrix(2,1)=-st
      matrix(1,2)=st
      matrix(2,2)=ct
      matrix(3,3)= ONE
    case default
      matrix(1:3,1:3)= ZERO
  end select
  
  return
  
 end subroutine matrix_rotaxis 
 
 subroutine eul2q(phis,psis,thetas,qs)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the quaternion associated to
!     the Euler angles
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in) :: phis,psis,thetas
  real(kind=PRC), intent(out), dimension(0:3) :: qs
  
  real(kind=PRC) :: cy,sy,cp,sp,cr,sr
  
  cy = cos(phis * HALF)
  sy = sin(phis * HALF)
  cp = cos(psis * HALF)
  sp = sin(psis * HALF)
  cr = cos(thetas * HALF)
  sr = sin(thetas * HALF)
  
  qs(0) = cy * cp * cr + sy * sp * sr
  qs(1) = cy * cp * sr - sy * sp * cr
  qs(2) = sy * cp * sr + cy * sp * cr
  qs(3) = sy * cp * cr - cy * sp * sr
  
  return
  
 end subroutine eul2q
 
 subroutine q2eul(qs,phis,psis,thetas)
 
!***********************************************************************
!     
!     LBsoft subroutine to compute the three Euler angles associated to
!     a quaternion
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification December 2018
!     
!***********************************************************************
 
  implicit none
  
  real(kind=PRC), intent(in), dimension(0:3) :: qs
  real(kind=PRC), intent(out) :: phis,psis,thetas
  
  real(kind=PRC) :: sinr_cosp,cosr_cosp,siny_cosp,cosy_cosp,sinp
  
  sinr_cosp = TWO * (qs(0) * qs(1) + qs(2) * qs(3))
  cosr_cosp = ONE- TWO * (qs(1)**TWO + qs(2)**TWO)
  thetas = atan2(sinr_cosp, cosr_cosp)

! psis (y-axis rotation)
  sinp = TWO * (qs(0) * qs(2) - qs(3) * qs(1))
  if(abs(sinp)>=ONE)then
    psis = sign(Pi*HALF, sinp) 
  else
    psis = asin(sinp)
  endif

! phis (z-axis rotation)
  siny_cosp = TWO * (qs(0) * qs(3) + qs(1) * qs(2))
  cosy_cosp = ONE - TWO * (qs(2)**TWO + qs(3)**TWO)
  phis = atan2(siny_cosp, cosy_cosp)
  
  return
  
 end subroutine q2eul
  
 end program create_xyz
 
 
  
