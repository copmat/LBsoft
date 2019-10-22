#define ONLYPOS
 program run_test
 
  implicit none
  
  integer, parameter :: ntest=15
  
  integer, parameter :: orig_io=1017
  integer, parameter :: orig_ioxyz=3017
  integer, parameter :: orig_mapio=18
  integer, parameter :: map_io=19
  integer, parameter :: actual_io=2027
  integer, parameter :: actual_ioxyz=4017
  integer, parameter :: remove_file_io=37
  character(len=*), parameter :: orig_file='output000000.map.orig' 
  character(len=*), parameter :: orig_xyzfile='restart.xyz.orig'
  character(len=*), parameter :: actual_file='output000000.map'
  character(len=*), parameter :: actual_xyzfile='restart.xyz'
  character(len=*), parameter :: origmap_file='global.map.orig'
  character(len=*), parameter :: map_file='global.map'
  integer, parameter :: maxlen=120
  
  character(len=maxlen),dimension(ntest) :: labeltest
  character(len=maxlen) :: origpath,newpath,temppath,lbsoftpath, &
   arg,directive,nomefile
  character(len=1) :: delimiter
  integer :: imio1,imio2,inumchar
  integer, dimension(3) :: imio3,imio4
  double precision, dimension(:), allocatable :: arr1,arr2
  integer :: mxrank1,mxrank2,nx1,ny1,nz1,nx2,ny2,nz2
  integer :: minx1,maxx1,miny1,maxy1,minz1,maxz1,type1,iselect1
  integer :: minx2,maxx2,miny2,maxy2,minz2,maxz2,type2,iselect2
  
  logical :: ltest(ntest)
  logical :: lxyzfile(ntest)
  logical :: lskip(ntest)
  logical :: lsingle_fluid(ntest)
  
  integer :: ndims(3,ntest)
  integer :: nparticles(ntest)
  
  integer :: i,j,k,l,ii,jj,kk,idrank,ncpu,ndec,xdim,ydim,zdim,iii, &
   nmio,ifield

  character(len=4) :: atmfield(25)
  character(len=256) :: file_loc_proc
  logical :: lmpi=.false.
  character(len=8) :: mystring8,atmname
  logical :: lexist
  character(len=40) :: mystring40
  
  character(len=5), parameter :: mysave = '.save'
  character(len=5), parameter :: myreff = '.reff'
  
  character(len=maxlen) :: namerestfile(8)
  
  double precision, parameter :: mytol=1.d-14
  
  integer :: natms1,natms2,iatm,nargxyz1,nargxyz2,timestep1,timestep2
  integer :: ifile
  logical :: lrotate(ntest)
  
  character(len=maxlen) :: myorigfile
  character(len=maxlen) :: myactufile
  
  integer(kind=1) :: isf1,isf2
  
  logical :: lforcestop,lerase
  
  double precision :: rtemp1,rtemp2
  
  double precision :: testreport(3,8,ntest)
  
  if(iargc()/=0 .and. iargc()/=6)then
    write(6,*) 'error!'
    write(6,*) 'the command line should be'
    write(6,*) '[executable] [-mpi] [ncpu] [ndec] [xdim] [ydim] [zdim]'
    write(6,*) '-mpi  = if present the test will be done on the mpi version'
    write(6,*) 'ncpu  = number of cpus (only if -mpi is present)'
    write(6,*) 'ndec  = type of decomposition'
    write(6,*) 'xdim = dimension of the decomposition along x'
    write(6,*) 'ydim = dimension of the decomposition along y'
    write(6,*) 'zdim = dimension of the decomposition along z'
    write(6,*) 'STOP!'
    stop
  endif
  
  
  do i = 1, iargc()
    call getarg(i, arg)
    if(i==1)then
      call copystring(arg,directive,maxlen)
      nomefile=repeat(' ',maxlen)
      nomefile=trim(directive)
      call lowcase(nomefile,maxlen)
      if(nomefile(1:4)=='-mpi')lmpi=.true.
      if(lmpi)then
        write(6,*) 'mpi job = YES'
      else
        write(6,*) 'mpi job = NOK'
      endif
    elseif(i==2)then
      call copystring(arg,directive,maxlen)
      ncpu=intstr(directive,maxlen,inumchar)
      write(6,*) 'ncpu = ',ncpu
    elseif(i==3)then
      call copystring(arg,directive,maxlen)
      ndec=intstr(directive,maxlen,inumchar)
      write(6,*) 'ndec = ',ndec
    elseif(i==4)then
      call copystring(arg,directive,maxlen)
      xdim=intstr(directive,maxlen,inumchar)
      write(6,*) 'xdim = ',xdim
    elseif(i==5)then
      call copystring(arg,directive,maxlen)
      ydim=intstr(directive,maxlen,inumchar)
      write(6,*) 'ydim = ',ydim
    elseif(i==6)then
      call copystring(arg,directive,maxlen)
      zdim=intstr(directive,maxlen,inumchar)
      write(6,*) 'zdim = ',zdim
    endif
  enddo
  
  
  if(lmpi)then
    write(mystring8,'(i8)')ncpu
  endif
  
  labeltest(1:ntest)=repeat(' ',maxlen)
  origpath=repeat(' ',maxlen)
  
  namerestfile(1:8)=repeat(' ',maxlen)
  
  namerestfile(1)='dumpGlobal'
  namerestfile(2)='dumpisFluid'
  namerestfile(3)='dumpPart'
  namerestfile(4)='dumpRhoR'
  namerestfile(5)='dumpRhoB'
  namerestfile(6)='dumpU'
  namerestfile(7)='dumpV'
  namerestfile(8)='dumpW'
  
  atmfield(1) ='xxx '
  atmfield(2) ='yyy '
  atmfield(3) ='zzz '
  atmfield(4) ='xxo '
  atmfield(5) ='yyo '
  atmfield(6) ='zzo '
  atmfield(7) ='vxx '
  atmfield(8) ='vyy '
  atmfield(9) ='vzz '
  atmfield(10)='vxo '
  atmfield(11)='vyo '
  atmfield(12)='vzo '
  atmfield(13)='fxbo'
  atmfield(14)='fybo'
  atmfield(15)='fzbo'
  atmfield(16)='q0  '
  atmfield(17)='q1  '
  atmfield(18)='q2  '
  atmfield(19)='q3  '
  atmfield(20)='oxx '
  atmfield(21)='oyy '
  atmfield(22)='ozz '
  atmfield(23)='txbo'
  atmfield(24)='tybo'
  atmfield(25)='tzbo'
  
  
  labeltest(1)='2D_Poiseuille_gradP_xy'
  labeltest(2)='2D_Poiseuille_gradP_xz'
  labeltest(3)='2D_Poiseuille_gradP_yz'
  labeltest(4)='2D_Poiseuille_xy'
  labeltest(5)='2D_Poiseuille_xz'
  labeltest(6)='2D_Poiseuille_yz'
  labeltest(7)='2D_Shear_xy'
  labeltest(8)='2D_Shear_xz'
  labeltest(9)='2D_Shear_yz'
  labeltest(10)='3D_Spinodal'
  labeltest(11)='3D_Particle_pbc'
  labeltest(12)='3D_Rotating_Particle'
  labeltest(13)='3D_Shear_Particle'
  labeltest(14)='3D_Particle_SC'
  labeltest(15)='3D_Particles_box'
  
  testreport(1:3,1:8,1:ntest)=0.d0
  
  ndims(1:3,14)=32
  ndims(1:3,15)=64
  
  nparticles(1:10)=0
  nparticles(11:14)=1
  nparticles(15)=137
  
  lrotate(1:11)=.false.
  lrotate(12:15)=.true.
  
  lsingle_fluid(1:ntest)=.true.
  lsingle_fluid(10)=.false.
  lsingle_fluid(14:15)=.false.
  
  lxyzfile(1:10)=.false.
  lxyzfile(11:ntest)=.true.
  lskip(1:ntest)=.true.
  lskip(1:ntest)=.false.
  !lskip(15)=.false.
  lerase=.true.
  lforcestop=.false.
  
  call getcwd(origpath)
  origpath=trim(origpath)
  delimiter=origpath(1:1)
  
  lbsoftpath=repeat(' ',maxlen)
  
  
  if(lmpi)then
    lbsoftpath='mpirun -np '//adjustl(trim(mystring8))//' ..'//delimiter// &
     '..'//delimiter//'execute'//delimiter//'main_mpi.x'
  else
    lbsoftpath=' ..'//delimiter//'..'//delimiter// &
     'execute'//delimiter//'main.x'
  endif
  write(6,*)'execute command: ',trim(lbsoftpath)
  
  ltest(1:ntest)=.true.
  do i=1,ntest
    if(lskip(i))cycle
    call chdir(trim(origpath))
    newpath=repeat(' ',maxlen)
    newpath='..'//delimiter//labeltest(i)
    !write(6,*)'entro in: ',trim(newpath)
    call chdir(trim(newpath))
    if(lmpi)then
      open(unit=42,file='mysed.sh',status='replace')
      write(42,*)'#!/bin/bash'
      write(42,'(a,i2,a)')'sed -i "/decomposition type/ c\decomposition type ', &
       ndec,'" input.dat'
      write(42,'(a,3i4,a)')'sed -i "/decomposition dimension/ c\decomposition dimension ', &
       xdim,ydim,zdim,'" input.dat'
      write(42,*)'exit 0'
      close(42)
      call execute_command_line('chmod 777 mysed.sh',wait=.true.)
      call execute_command_line('./mysed.sh',wait=.true.)
    endif
    temppath=repeat(' ',maxlen)
    call getcwd(temppath)
    !write(6,*)'sono in: ',trim(temppath)
    write(6,'(2a)')'execute test: ',labeltest(i)
    call execute_command_line(trim(lbsoftpath),wait=.true.)
    write(6,'(2a)')'check test  : ',labeltest(i)
    
    !cycle
    
    ifile=1
    myorigfile=repeat(' ',maxlen)
    myactufile=repeat(' ',maxlen)
    myorigfile=trim(namerestfile(ifile))//myreff//'.dat'
    myactufile=trim(namerestfile(ifile))//mysave//'.dat'
    
    open(unit=orig_io,file=trim(myorigfile),status='old',action='read', &
     form='unformatted')
    read(orig_io)timestep1
    close(orig_io)
    
    inquire(file=trim(myactufile),exist=lexist)
    if(.not.lexist)then
      ltest(i)=.false.
      write(6,'(3a)')'file ',trim(myactufile),' not found!'
      !goto 100
    endif
    open(unit=actual_io,file=trim(myactufile),status='old',action='read', &
     form='unformatted')
    read(actual_io)timestep2
    close(actual_io)
    
    if(timestep1/=timestep2)then
      ltest(i)=.false.
      testreport(1,ifile,i)=1.d+6
      testreport(2,ifile,i)=dble(ifile)
      goto 100
    endif
    
    
    if(lxyzfile(i))then
      ifile=3
      myorigfile=repeat(' ',maxlen)
      myactufile=repeat(' ',maxlen)
      myorigfile=trim(namerestfile(ifile))//myreff//'.dat'
      myactufile=trim(namerestfile(ifile))//mysave//'.dat'
      
      nmio=nparticles(i)
      if(allocated(arr1))deallocate(arr1)
      if(allocated(arr2))deallocate(arr2)
      allocate(arr1(nmio),arr2(nmio))
      
      open(unit=orig_io,file=trim(myorigfile),status='old',action='read', &
       form='unformatted')
      
      inquire(file=trim(myactufile),exist=lexist)
      if(.not.lexist)then
        ltest(i)=.false.
        write(6,'(3a)')'file ',trim(myactufile),' not found!'
        !goto 100
      endif
      open(unit=actual_io,file=trim(myactufile),status='old',action='read', &
       form='unformatted')
       
#ifdef ONLYPOS
      do ifield=1,3
#else
      do ifield=1,15
#endif
        read(orig_io)(arr1(ii),ii=1,nmio)
        read(actual_io)(arr2(ii),ii=1,nmio)
        do ii=1,nmio
          if(abs(arr1(ii)-arr2(ii))>mytol)then
            ltest(i)=.false.
            testreport(1,ifile,i)=max(abs(arr1(ii)-arr2(ii)),testreport(1,ifile,i))
            if(testreport(1,ifile,i)==abs(arr1(ii)-arr2(ii)))then
              testreport(2,ifile,i)=dble(ifile)
              testreport(3,ifile,i)=dble(ifield)
            endif
            if(lforcestop)goto 120
            !goto 100
          endif
        enddo
      enddo
#ifndef ONLYPOS
      if(lrotate(i))then
        do ifield=16,25
          read(orig_io)(arr1(ii),ii=1,nmio)
          read(actual_io)(arr2(ii),ii=1,nmio)
          do ii=1,nmio
            if(abs(arr1(ii)-arr2(ii))>mytol)then
              !ltest(i)=.false.
              testreport(1,ifile,i)=max(abs(rtemp1-rtemp2),testreport(1,ifile,i))
              if(testreport(1,ifile,i)==abs(rtemp1-rtemp2))then
                testreport(2,ifile,i)=dble(ifile)
                testreport(3,ifile,i)=dble(ifield)
              endif
              if(lforcestop)goto 120
              !goto 100
            endif
          enddo
        enddo
      endif
#endif
      close(orig_io)
      close(actual_io)
      
    endif
    
    
    
    ifile=4
    myorigfile=repeat(' ',maxlen)
    myactufile=repeat(' ',maxlen)
    myorigfile=trim(namerestfile(ifile))//myreff//'.dat'
    myactufile=trim(namerestfile(ifile))//mysave//'.dat'
    
    open(unit=orig_io,file=trim(myorigfile),status='old',action='read', &
     form='unformatted',access='stream')
    
    inquire(file=trim(myactufile),exist=lexist)
    if(.not.lexist)then
      ltest(i)=.false.
      write(6,'(3a)')'file ',trim(myactufile),' not found!'
      !goto 100
    endif
    open(unit=actual_io,file=trim(myactufile),status='old',action='read', &
     form='unformatted',access='stream')
    
    do kk=1,ndims(3,i)
      do jj=1,ndims(2,i)
        do ii=1,ndims(1,i)
          read(orig_io)rtemp1
          read(actual_io)rtemp2
          if(abs(rtemp1-rtemp2)>mytol)then
            ltest(i)=.false.
            testreport(1,ifile,i)=max(abs(rtemp1-rtemp2),testreport(1,ifile,i))
            if(testreport(1,ifile,i)==abs(rtemp1-rtemp2))then
              testreport(2,ifile,i)=dble(ifile)
            endif
            if(lforcestop)goto 110
            !goto 100
          endif
        enddo
      enddo
    enddo
    
    if(.not. lsingle_fluid(i))then
      ifile=5
      myorigfile=repeat(' ',maxlen)
      myactufile=repeat(' ',maxlen)
      myorigfile=trim(namerestfile(ifile))//myreff//'.dat'
      myactufile=trim(namerestfile(ifile))//mysave//'.dat'
      
      open(unit=orig_io,file=trim(myorigfile),status='old',action='read', &
       form='unformatted',access='stream')
      
      inquire(file=trim(myactufile),exist=lexist)
      if(.not.lexist)then
        ltest(i)=.false.
        write(6,'(3a)')'file ',trim(myactufile),' not found!'
        !goto 100
      endif
      open(unit=actual_io,file=trim(myactufile),status='old',action='read', &
       form='unformatted',access='stream')
      
      do kk=1,ndims(3,i)
        do jj=1,ndims(2,i)
          do ii=1,ndims(1,i)
            read(orig_io)rtemp1
            read(actual_io)rtemp2
            if(abs(rtemp1-rtemp2)>mytol)then
              ltest(i)=.false.
              testreport(1,ifile,i)=max(abs(rtemp1-rtemp2),testreport(1,ifile,i))
              if(testreport(1,ifile,i)==abs(rtemp1-rtemp2))then
                testreport(2,ifile,i)=dble(ifile)
              endif
              if(lforcestop)goto 110
              !goto 100
            endif
          enddo
        enddo
      enddo
    endif
    
    ifile=6
    myorigfile=repeat(' ',maxlen)
    myactufile=repeat(' ',maxlen)
    myorigfile=trim(namerestfile(ifile))//myreff//'.dat'
    myactufile=trim(namerestfile(ifile))//mysave//'.dat'
    
    open(unit=orig_io,file=trim(myorigfile),status='old',action='read', &
     form='unformatted',access='stream')
    
    inquire(file=trim(myactufile),exist=lexist)
    if(.not.lexist)then
      ltest(i)=.false.
      write(6,'(3a)')'file ',trim(myactufile),' not found!'
      !goto 100
    endif
    open(unit=actual_io,file=trim(myactufile),status='old',action='read', &
     form='unformatted',access='stream')
    
    do kk=1,ndims(3,i)
      do jj=1,ndims(2,i)
        do ii=1,ndims(1,i)
          read(orig_io)rtemp1
          read(actual_io)rtemp2
          if(abs(rtemp1-rtemp2)>mytol)then
            ltest(i)=.false.
            testreport(1,ifile,i)=max(abs(rtemp1-rtemp2),testreport(1,ifile,i))
            if(testreport(1,ifile,i)==abs(rtemp1-rtemp2))then
              testreport(2,ifile,i)=dble(ifile)
            endif
            if(lforcestop)goto 110
            !goto 100
          endif
        enddo
      enddo
    enddo
    
    ifile=7
    myorigfile=repeat(' ',maxlen)
    myactufile=repeat(' ',maxlen)
    myorigfile=trim(namerestfile(ifile))//myreff//'.dat'
    myactufile=trim(namerestfile(ifile))//mysave//'.dat'
    
    open(unit=orig_io,file=trim(myorigfile),status='old',action='read', &
     form='unformatted',access='stream')
    
    inquire(file=trim(myactufile),exist=lexist)
    if(.not.lexist)then
      ltest(i)=.false.
      write(6,'(3a)')'file ',trim(myactufile),' not found!'
      !goto 100
    endif
    open(unit=actual_io,file=trim(myactufile),status='old',action='read', &
     form='unformatted',access='stream')
    
    do kk=1,ndims(3,i)
      do jj=1,ndims(2,i)
        do ii=1,ndims(1,i)
          read(orig_io)rtemp1
          read(actual_io)rtemp2
          if(abs(rtemp1-rtemp2)>mytol)then
            ltest(i)=.false.
            testreport(1,ifile,i)=max(abs(rtemp1-rtemp2),testreport(1,ifile,i))
            if(testreport(1,ifile,i)==abs(rtemp1-rtemp2))then
              testreport(2,ifile,i)=dble(ifile)
            endif
            if(lforcestop)goto 110
            !goto 100
          endif
        enddo
      enddo
    enddo
    
    ifile=8
    myorigfile=repeat(' ',maxlen)
    myactufile=repeat(' ',maxlen)
    myorigfile=trim(namerestfile(ifile))//myreff//'.dat'
    myactufile=trim(namerestfile(ifile))//mysave//'.dat'
    
    open(unit=orig_io,file=trim(myorigfile),status='old',action='read', &
     form='unformatted',access='stream')
    
    inquire(file=trim(myactufile),exist=lexist)
    if(.not.lexist)then
      ltest(i)=.false.
      write(6,'(3a)')'file ',trim(myactufile),' not found!'
      !goto 100
    endif
    open(unit=actual_io,file=trim(myactufile),status='old',action='read', &
     form='unformatted',access='stream')
    
    do kk=1,ndims(3,i)
      do jj=1,ndims(2,i)
        do ii=1,ndims(1,i)
          read(orig_io)rtemp1
          read(actual_io)rtemp2
          if(abs(rtemp1-rtemp2)>mytol)then
            ltest(i)=.false.
            testreport(1,ifile,i)=max(abs(rtemp1-rtemp2),testreport(1,ifile,i))
            if(testreport(1,ifile,i)==abs(rtemp1-rtemp2))then
              testreport(2,ifile,i)=dble(ifile)
            endif
            if(lforcestop)goto 110
            !goto 100
          endif
        enddo
      enddo
    enddo
    
100 continue
    
    close(orig_io)
    close(actual_io)
 
        
    open(unit=remove_file_io,file='rm.sh',status='replace')
    
    write(remove_file_io,'(a)')'#!/bin/bash'
    write(remove_file_io,'(a)')'if [ -f "time.dat" ]; then rm time.dat; fi'
    write(remove_file_io,'(a)')'if [ -f "mysed.sh" ]; then rm mysed.sh; fi'
    write(remove_file_io,'(a)')'if [ -f "restart.xyz" ]; then rm restart.xyz; fi'
    write(remove_file_io,'(a)')'if [ -f "restart.xyz.orig" ]; then rm restart.xyz.orig; fi'
    write(remove_file_io,'(a)')'if [ -f "global.map" ]; then rm global.map; fi'
    write(remove_file_io,'(a)')'if [ -f "global.map.orig" ]; then rm global.map.orig; fi'
    write(remove_file_io,'(a)')'if [ -f "dumpPopsR.save.dat" ]; then rm dumpPopsR.save.dat; fi'
    write(remove_file_io,'(a)')'if [ -f "dumpPopsB.save.dat" ]; then rm dumpPopsB.save.dat; fi'
    do ifile=1,8
      myactufile=trim(namerestfile(ifile))//mysave//'.dat'
      write(remove_file_io,'(5a)')'if [ -f "',trim(myactufile), &
       '" ]; then rm ',trim(myactufile),'; fi'
    enddo
    
    !write(remove_file_io,'(a)')'rm -f output*.map'
    !write(remove_file_io,'(a)')'rm -f output*.map.*'
    write(remove_file_io,'(a)')'if [ -f "rm.sh" ]; then rm rm.sh; fi'
    

    close(remove_file_io)
    call execute_command_line('chmod 777 rm.sh',wait=.true.)
    if(lerase)call execute_command_line('./rm.sh',wait=.true.)
    if(ltest(i))then
      write(6,'(3a)')'test ',labeltest(i),' : OK'
    else
     write(6,'(3a)')'test ',labeltest(i),' : NOK'
    endif
    cycle
    
110 continue
    write(6,*)'error test quantity: ',rtemp1,rtemp2
    write(6,*)'error test position: ',ii, jj, kk
    write(6,*)'error test file: ',trim(myactufile)
    stop
    
120 continue
    write(6,*)'error test quantity: ',arr1(ii),arr2(ii)
    write(6,*)'error test position: ',ii,ifield
    write(6,*)'error test file: ',trim(myactufile)
    stop
   
  enddo
  
  call chdir(trim(origpath))
  
  write(6,'(/,a,/)')repeat('*',80)
  write(6,'(a)')'list of test results: '
  do i=1,ntest
    if(ltest(i))then
      write(6,'(3a)')'test ',labeltest(i),' : OK'
    else
      write(6,'(3a)')'test ',labeltest(i),' : NOK'
      write(6,'(a)')'ERROR REPORT:'
      do ifile=3,8
        write(6,'(a,es16.8)')'max error = ',testreport(1,ifile,i)
        write(6,'(2a)')'in file : ',namerestfile(ifile)
        if(ifile==3)then
          write(6,'(2a)')'in field : ',atmfield(nint(testreport(3,ifile,i)))
        endif
      enddo
    endif
  enddo
  write(6,'(/,a,/)')repeat('*',80)
  stop
  
  
 contains
 
 function dimenumb(inum)
 
!***********************************************************************
!     
!     LBsoft function for returning the number of digits
!     of an integer number
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************

  implicit none

  integer,intent(in) :: inum
  integer :: dimenumb
  integer :: i
  real(kind=4) :: tmp

  i=1
  tmp=real(inum,kind=4)
  do 
    if(tmp< 10.e0 )exit
    i=i+1
    tmp=tmp/ 10.e0
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
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
 
  implicit none

  integer,intent(in) :: inum
  character(len=6) :: write_fmtnumb
  integer :: numdigit,irest
  real(kind=4) :: tmp
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
 
 function intstr(string,lenstring,laststring)
 
!***********************************************************************
!     
!     LBsoft subroutine for extracting integers from a character 
!     string
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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
!     LBsoft subroutine for extracting double precisions from a  
!     character string
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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

 subroutine copystring(oldstring,newstring,lenstring)
 
!***********************************************************************
!     
!     LBsoft subroutine to copy one character string into another
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
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
  
 function findstring(seek,string,here,lenstring)
 
!***********************************************************************
!     
!     LBsoft subroutine to find an explicit string in an input record
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification March 2015
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(in) :: seek
  character(len=*), intent(inout) :: string
  integer, intent(out) :: here
  integer, intent(in) :: lenstring
  
  logical :: findstring
  character*1, dimension(lenstring) :: word
  integer :: i,j,nseek,m
  logical :: findspace
  
  do j=1,lenstring
    word(j)=string(j:j)
  enddo 

  m=lenstring
  nseek=len(seek)
  findstring=.false.
  
  here=0
  findspace=.true.
  do while(here<m-nseek .and. (.not.findstring))
    
    findstring=.true.
    
    do i=1,nseek
      if(seek(i:i)/=word(here+i))findstring=.false.
    enddo
    findstring=(findstring.and.findspace)

    here=here+1
    
    if(word(here)==' ')then
      findspace=.true.
    else
      findspace=.false.
    endif

  enddo
  
  do j=1,lenstring
    string(j:j)=word(j)
  enddo

  return
  
 end function findstring
 
 subroutine strip(string,lenstring)
 
!***********************************************************************
!     
!     LBsoft subroutine to strip a character string
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************

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

 subroutine lowcase(string,lenstring)
 
!***********************************************************************
!     
!     LBsoft subroutine to lowercase a character string
!     originally written in JETSPIN by M. Lauricella et al.
!     
!     licensed under Open Software License v. 3.0 (OSL-3.0)
!     author: M. Lauricella
!     last modification July 2018
!     
!***********************************************************************
  
  implicit none
  
  character(len=*), intent(inout) :: string
  integer, intent(in) :: lenstring
  
  character*1, dimension(lenstring) :: word
  character*1 :: letter
  
  integer :: i,j
  
  do j=1,lenstring
    word(j)=string(j:j)
  enddo 
  
  do i=1,lenstring
    
    letter=word(i)
    
    if(letter=='A')then
      letter='a'
    else if(letter=='B')then
      letter='b'
    else if(letter=='C')then
      letter='c'
    else if(letter=='D')then
      letter='d'
    else if(letter=='E')then
      letter='e'
    else if(letter=='F')then
      letter='f'
    else if(letter=='G')then
      letter='g'
    else if(letter=='H')then
      letter='h'
    else if(letter=='I')then
      letter='i'
    else if(letter=='J')then
      letter='j'
    else if(letter=='K')then
      letter='k'
    else if(letter=='L')then
      letter='l'
    else if(letter=='M')then
      letter='m'
    else if(letter=='N')then
      letter='n'
    else if(letter=='O')then
      letter='o'
    else if(letter=='P')then
      letter='p'
    else if(letter=='Q')then
      letter='q'
    else if(letter=='R')then
      letter='r'
    else if(letter=='S')then
      letter='s'
    else if(letter=='T')then
      letter='t'
    else if(letter=='U')then
      letter='u'
    else if(letter=='V')then
      letter='v'
    else if(letter=='W')then
      letter='w'
    else if(letter=='X')then
      letter='x'
    else if(letter=='Y')then
      letter='y'
    else if(letter=='Z')then
      letter='z'
    endif
    
    word(i)=letter
    
  enddo
  
  do j=1,lenstring
    string(j:j)=word(j)
  enddo
  
  return
  
 end subroutine lowcase
  
 end program run_test
 
 
  
