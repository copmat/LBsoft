
 program extract_vti
 
  use parse_module
  
  implicit none
  
  integer, parameter :: totmaxfinal=30000
  
  real, parameter :: & 
  Pi=real(3.141592653589793238462643383279502884d0,kind=4)
  
  logical :: lorthogonal=.true.
  logical :: ltriclinic=.false.
  
  integer, parameter :: maxlen=120
  integer, parameter :: maxfilevtk=5
  real, parameter :: zero=real(0.d0)
  character(len=*), parameter :: rootvtk='out_'
  character(len=*), parameter :: filenamevtk='out'
  character(len=8),dimension(maxfilevtk) :: filenamevar
  character(len=maxlen) :: tempstring,myfile
  character(len=maxlen) :: arg,redstring,directive,myarg
  character(len=maxlen),dimension(maxfilevtk) :: extentvtk,namefilein
  logical :: safe
  integer :: idum,natms,i,ifile
  
  logical :: lexist
  
  real, allocatable, dimension(:,:,:) :: rhoR,rhoB,phase,u,v,w
  
  logical :: lendmio
  integer :: im,iw,tramutare,nworig,nmorig,j,k,ii,jj,kk,istep,nx,ny,nz,inumchar


  logical, parameter :: lverbose=.false.
  integer :: istart,iend,ievery,nfl_tot_max,nfl_tot_temp
  
  integer, parameter :: ifileio=222
  integer :: myIPRC,myPRC,myPRCWR,numprocs,numprocstemp
  integer :: inx,enx,iny,eny,inz,enz,enn
  
  integer,dimension(maxfilevtk) :: ndimvtk,vtkcase
  character(len=4) :: string4
  logical :: ltest,lelittle
  character(len=500) :: headervtk
  
  integer :: indent,myoffset,new_myoffset,nfilevtk,nn
  integer(kind=4) :: nbyte
  real(kind=4), allocatable, dimension(:) :: myvar1d,myvar2d,myvar3d
  
  nx=0 
  ny=0
  nz=0
  nn=0
  if(iargc()<5)then
    write(6,*) '#error!'
    write(6,*) '#the command line should be'
    write(6,*) '#[executable] [start] [end] [every] [nfile] [file]'
    write(6,*) '#start = integer indicating the first frame to analyze'
    write(6,*) '#end   = integer indicating the last frame to analyze'
    write(6,*) '#every = integer indicating the interval number of frame to analyze'
    write(6,*) '#nfile = integer indicating number of arguments'
    write(6,*) '#file  = string indicating the input file'
    write(6,*) '#STOP!'
    goto 20
  endif
  
  do i = 1, iargc()
    call getarg(i,myarg)
    if(i==1)then
      call copystring(myarg,directive,maxlen)
      istart=intstr(directive,maxlen,inumchar)
      write(6,*) '#start = ',istart
    elseif(i==2)then
      call copystring(myarg,directive,maxlen)
      iend=intstr(directive,maxlen,inumchar)
      write(6,*) '#end   = ',iend
    elseif(i==3)then
      call copystring(myarg,directive,maxlen)
      ievery=intstr(directive,maxlen,inumchar)
      write(6,*) '#every = ',ievery
      
    elseif(i==4)then
      call copystring(myarg,directive,maxlen)
      nfilevtk=intstr(directive,maxlen,inumchar)
      write(6,*) '#narg  = ',nfilevtk
      if(nfilevtk>maxfilevtk)then
        write(6,*)'ERROR narg should be less than 5'
        goto 20
      endif
    elseif(i>=5)then
      call copystring(myarg,directive,maxlen)
      namefilein(i-4)=repeat(' ',maxlen)
      call copystring(myarg,namefilein(i-4),maxlen)
      call lowcase(namefilein(i-4),maxlen)
      call strip(namefilein(i-4),maxlen)
      write(6,*) '#file   = ',namefilein(i-4)
    endif
  enddo
  
  
  call init_random_seed(317)
  
  call test_little_endian(lelittle)
  
  do j=1,nfilevtk
    tempstring=namefilein(j)
    ltest=.false.
    do i=1,4
      if(tempstring(i:i)/=rootvtk(i:i))then
        ltest=.true.
        exit
      endif
    enddo
    if(ltest)goto 113
    string4(1:4)=tempstring(5:8)
  
    select case(string4)
    case('rho1')
      ndimvtk(j)=1
      vtkcase(j)=1
      filenamevar(j)='rho1    '
    case('rho2')
      ndimvtk(j)=1
      vtkcase(j)=2
      filenamevar(j)='rho2    '
    case('vel_')
      ndimvtk(j)=3
      vtkcase(j)=3
      filenamevar(j)='vel     '
    case('phas')
      ndimvtk(j)=1
      vtkcase(j)=4
      filenamevar(j)='phase   '
    case default
      goto 113
    end select
  enddo
  
  
  do ifile=istart,iend,ievery
    do jj=1,nfilevtk
      myfile=repeat(' ',maxlen)
      myfile=trim(filenamevtk)//'_'//trim(filenamevar(jj))// &
     '_'//trim(write_fmtnumb8(ifile)) // '.vti'
    
      write(6,'(2a)') 'serching file: ',trim(myfile)
      inquire(file=trim(myfile),exist=lexist)
      if(.not.lexist)then
        write(6,'(3a)') 'file: ',trim(myfile),' NOT FOUND'
        write(6,'(2a)') 'skip file ',trim(myfile)
        cycle
      endif
      write(6,'(2a)') 'reading file :',trim(myfile)
      open(unit=ifileio,file=trim(myfile),FORM='UNFORMATTED', &
       STATUS='old', ACTION='read', ACCESS='stream')
      myoffset=0
      indent=0
      call read_header_vtk(ifileio,headervtk,filenamevar(jj),extentvtk(jj), &
      ndimvtk(jj),0,iend,myoffset,new_myoffset,indent,inx,enx,iny,eny,inz,enz)
    
      
      
      enn=enx*eny*enz
      
      if(inx /=1)goto 20
      if(iny /=1)goto 20
      if(inz /=1)goto 20
      if(enx/=nx .or. eny/=ny .or. enz/=nz .or. enn/=nn)then
        do kk=1,nfilevtk
          select case(vtkcase(kk))
          case(1)
            if(allocated(rhoR))deallocate(rhoR)
            nx=enx
            ny=eny
            nz=enz
            allocate(rhoR(1:nx,1:ny,1:nz))
          case(2)
            if(allocated(rhoB))deallocate(rhoB)
            nx=enx
            ny=eny
            nz=enz
            allocate(rhoB(1:nx,1:ny,1:nz))
          case(3)
            if(allocated(u))deallocate(u)
            if(allocated(v))deallocate(v)
            if(allocated(w))deallocate(w)
            nx=enx
            ny=eny
            nz=enz
            allocate(u(1:nx,1:ny,1:nz))
            allocate(v(1:nx,1:ny,1:nz))
            allocate(w(1:nx,1:ny,1:nz))
          case(4)
            if(allocated(phase))deallocate(phase)
            nx=enx
            ny=eny
            nz=enz
            allocate(phase(1:nx,1:ny,1:nz))
          case default
           goto 113
          end select
        enddo
        if(any(vtkcase(1:nfilevtk)==3))then
          nn=nx*ny*nz
          if(allocated(myvar1d))deallocate(myvar1d)
          allocate(myvar1d(1:nn))
          if(allocated(myvar2d))deallocate(myvar2d)
          allocate(myvar2d(1:nn))
          if(allocated(myvar3d))deallocate(myvar3d)
          allocate(myvar3d(1:nn))
        else
          nn=nx*ny*nz
          if(allocated(myvar1d))deallocate(myvar1d)
          allocate(myvar1d(1:nn))
        endif
      endif
      
      select case(vtkcase(jj))
      case(1)
        read(ifileio)nbyte,(myvar1d(ii),ii=1,nn)
        ii=0
        do k=1,nz
          do j=1,ny
            do i=1,nx
              ii=ii+1
              rhoR(i,j,k)=dble(myvar1d(ii))
            enddo
          enddo
        enddo
      case(2)
        read(ifileio)nbyte,(myvar1d(ii),ii=1,nn)
        ii=0
        do k=1,nz
          do j=1,ny
            do i=1,nx
              ii=ii+1
              rhoB(i,j,k)=dble(myvar1d(ii))
            enddo
          enddo
        enddo
      case(3)
        read(ifileio)nbyte,(myvar1d(ii),myvar2d(ii), &
         myvar3d(ii),ii=1,nn)
        ii=0
        do k=1,nz
          do j=1,ny
            do i=1,nx
              ii=ii+1
              u(i,j,k)=dble(myvar1d(ii))
              v(i,j,k)=dble(myvar2d(ii))
              w(i,j,k)=dble(myvar3d(ii))
            enddo
          enddo
        enddo
      case(4)
        read(ifileio)nbyte,(myvar1d(ii),ii=1,nn)
        ii=0
        do k=1,nz
          do j=1,ny
            do i=1,nx
              ii=ii+1
              phase(i,j,k)=dble(myvar1d(ii))
            enddo
          enddo
        enddo
      case default
       goto 113
      end select
      
      close(ifileio)
      
    enddo
    
    call compute_analysis
    
  enddo
  
  stop
  
  
  write(6,'(a)')'Correctly closed!'
  
  stop
  
 20 continue
  
  write(6,'(a)')'Closed with error!'
  
  stop
  
 113 continue
  
  write(6,'(2a)')'Error in reading the input file: the file name is not regular'
  write(6,'(2a)')'the file name should be out_[var]_############.vti',&
   'where [var]=rho1,rho2,phase,vel and each # is an integer (twelve integers)'
  
  stop
  
 contains
 
 subroutine compute_analysis
 
!***********************************************************************
!     
!     PUT HERE YOUR ANALYSIS
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  return
  
 end subroutine compute_analysis
 
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
  
  integer, parameter :: idrank=0
          
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
 
 function dimenumb(inum)
 
!***********************************************************************
!     
!     LBsoft function for returning the number of digits
!     of an integer number
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************

  implicit none

  integer,intent(in) :: inum
  integer :: dimenumb
  integer :: i
  real(kind=8) :: tmp

  i=1
  tmp=real(inum,kind=8)
  do 
    if(tmp< 10.d0 )exit
    i=i+1
    tmp=tmp/ 10.d0
  enddo

  dimenumb=i

  return

 end function dimenumb

 function write_fmtnumb(inum)
 
!***********************************************************************
!     
!     LBsoft function for returning the string of six characters 
!     with integer digits and leading zeros to the left
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none

  integer,intent(in) :: inum
  character(len=6) :: write_fmtnumb
  integer :: numdigit,irest
  real(kind=8) :: tmp
  character(len=22) :: cnumberlabel

  numdigit=dimenumb(inum)
  irest=6-numdigit
  if(irest>0)then
    write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
    write(write_fmtnumb,fmt=cnumberlabel)repeat('0',irest),inum
  else
    write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
    write(write_fmtnumb,fmt=cnumberlabel)inum
  endif
  
  return

 end function write_fmtnumb
      
 function write_fmtnumb8(inum)
 
!***********************************************************************
!     
!     LBsoft function for returning the string of six characters 
!     with integer digits and leading zeros to the left
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none

  integer,intent(in) :: inum
  character(len=8) :: write_fmtnumb8
  integer :: numdigit,irest
  real(kind=8) :: tmp
  character(len=22) :: cnumberlabel

  numdigit=dimenumb(inum)
  irest=8-numdigit
  if(irest>0)then
    write(cnumberlabel,"(a,i8,a,i8,a)")"(a",irest,",i",numdigit,")"
    write(write_fmtnumb8,fmt=cnumberlabel)repeat('0',irest),inum
  else
    write(cnumberlabel,"(a,i8,a)")"(i",numdigit,")"
    write(write_fmtnumb8,fmt=cnumberlabel)inum
  endif
  
  return

 end function write_fmtnumb8
 
 subroutine read_header_vtk(iio,mystring500,namevar,extent,ncomp, &
  iinisub,iend,myoffset,new_myoffset,indent,inxtemp,nxtemp,inytemp, &
  nytemp,inztemp,nztemp)
  
  implicit none
  
  character(len=8),intent(in) :: namevar
  character(len=maxlen),intent(inout) :: extent
  integer, intent(in) :: iio,ncomp,iinisub,myoffset
  integer, intent(out) :: iend,new_myoffset
  integer, intent(inout) :: indent
  integer, intent(out) :: inxtemp,nxtemp,inytemp,nytemp,inztemp,nztemp
  !namevar='density1'
  
  character(len=500), intent(out) :: mystring500
  ! End-character for binary-record finalize.
  character(1), parameter:: end_rec = char(10) 
  character(1) :: string1
  character(len=*),parameter :: topology='ImageData' 
  integer :: ioffset,nele,bytechar,byteint,byter4,byter8,iini,i,j
  
  logical :: lstart,lact
  
  iini=iinisub
  bytechar=kind(end_rec)
  byteint=kind(iini)
  byter4  = 4
  byter8  = 8
  
  mystring500=repeat(' ',500)
  
  iend=iini
  
  iini=iend+1
  nele=22
  iend=iend+nele
  !mystring500(iini:iend)='<?xml version="1.0"?>'//end_rec
  read(iio)mystring500(iini:iend)
  
  new_myoffset=myoffset
  new_myoffset = new_myoffset + nele * bytechar
 
  
  iini=iend+1
  nele=67
  iend=iend+nele
!  if(lelittle)then  
!    mystring500(iini:iend) = '<VTKFile type="'//trim(topology)// &
!     '" version="0.1" byte_order="LittleEndian">'//end_rec
!  else
!    mystring500(iini:iend) = '<VTKFile type="'//trim(topology)// &
!     '" version="0.1" byte_order="BigEndian">   '//end_rec
!  endif
  read(iio)mystring500(iini:iend)
  
  new_myoffset = new_myoffset + 67 * bytechar
 
  
  indent = indent + 2
  iini=iend+1
  nele=70
  iend=iend+nele
!  mystring500(iini:iend) = repeat(' ',indent)//'<'//trim(topology)//' WholeExtent="'//&
!                 trim(extent)//'">'//end_rec
  read(iio)mystring500(iini:iend)

  new_myoffset = new_myoffset + 70 * bytechar
 
  
  indent = indent + 2
  iini=iend+1
  nele=63
  iend=iend+nele
!  mystring500(iini:iend) = repeat(' ',indent)//'<Piece Extent="'//trim(extent)//'">'//end_rec
  read(iio)mystring500(iini:iend)
  lstart=.false.
  extent=repeat(' ',maxlen)
  j=0
  do i=iini,iend
    if(mystring500(i:i)=='"')then
      if(lstart)then
        lact=.false.
      else
        lact=.true.
        lstart=.true.
      endif
    else
      if(lact)then
        j=j+1
        extent(j:j)=mystring500(i:i)
      endif
    endif
  enddo
  
  redstring=repeat(' ',maxlen)
  redstring(1:maxlen)=extent(1:maxlen)
  call strip(redstring,maxlen)
  call lowcase(redstring,maxlen)
  call copystring(redstring,directive,maxlen)
  inxtemp=intstr(directive,maxlen,idum)
  nxtemp=intstr(directive,maxlen,idum)
  inytemp=intstr(directive,maxlen,idum)
  nytemp=intstr(directive,maxlen,idum)
  inztemp=intstr(directive,maxlen,idum)
  nztemp=intstr(directive,maxlen,idum)
  
  new_myoffset = new_myoffset + 63 * bytechar
 
  
! initializing offset pointer
  ioffset = 0 
  
  indent = indent + 2
  iini=iend+1
  nele=18
  iend=iend+nele
  !mystring500(iini:iend)=repeat(' ',indent)//'<PointData>'//end_rec
  read(iio)mystring500(iini:iend)
  new_myoffset = new_myoffset + 18 * bytechar
  
  indent = indent + 2
  iini=iend+1
  nele=115
  iend=iend+nele
  
  if(ncomp/=1 .and. ncomp/=3)then
    write(6,'(a)')'ERROR in header_vtk'
    stop
  endif
  write(string1,'(i1)')ncomp
!  mystring500(iini:iend)=repeat(' ',indent)//'<DataArray type="Float32" Name="'// &
!   namevar//'" NumberOfComponents="'//string1// '" '//&
!   'format="appended" offset="'//space_fmtnumb12(ioffset)//'"/>'//end_rec
  read(iio)mystring500(iini:iend)
  new_myoffset = new_myoffset + 115 * bytechar
 
  
  indent = indent - 2
  iini=iend+1
  nele=19
  iend=iend+nele
!  mystring500(iini:iend)=repeat(' ',indent)//'</PointData>'//end_rec
  read(iio)mystring500(iini:iend)
  new_myoffset = new_myoffset + 19 * bytechar
  
  
  indent = indent - 2
  iini=iend+1
  nele=13
  iend=iend+nele
!  mystring500(iini:iend)=repeat(' ',indent)//'</Piece>'//end_rec
  read(iio)mystring500(iini:iend)
  
  new_myoffset = new_myoffset + 13 * bytechar
 
  
  indent = indent - 2
  iini=iend+1
  nele=15
  iend=iend+nele
!  mystring500(iini:iend)=repeat(' ',indent)//'</'//trim(topology)//'>'//end_rec
  read(iio)mystring500(iini:iend)
  new_myoffset = new_myoffset + 15 * bytechar
 

  iini=iend+1
  nele=32
  iend=iend+nele
!  mystring500(iini:iend)=repeat(' ',indent)//'<AppendedData encoding="raw">'//end_rec
  read(iio)mystring500(iini:iend)
  new_myoffset = new_myoffset + 32 * bytechar
  
  iini=iend+1
  nele=1
  iend=iend+nele
!  mystring500(iini:iend)='_'
  read(iio)mystring500(iini:iend)
  new_myoffset = new_myoffset + 1 * bytechar
  
  return
  
 end subroutine read_header_vtk
 
 subroutine test_little_endian(ltest)
 
!***********************************************************************
!     
!     LBsoft subroutine for checking if the computing architecture
!     is working in little-endian or big-endian
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************
 
  implicit none 
  integer, parameter :: ik1 = selected_int_kind(2) 
  integer, parameter :: ik4 = selected_int_kind(9) 
   
  logical, intent(out) :: ltest
   
  if(btest(transfer(int((/1,0,0,0/),ik1),1_ik4),0)) then 
    !it is little endian
    ltest=.true.
  else 
    !it is big endian
    ltest=.false.
  end if 
   
  return
   
 end subroutine test_little_endian 
 
 end program extract_vti
 
  
