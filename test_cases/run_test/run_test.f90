
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
  real(8), dimension(:,:,:,:), allocatable :: dmio1,dmio2
  real(8), dimension(:,:), allocatable :: dmio3,dmio4
  integer :: mxrank1,mxrank2,nx1,ny1,nz1,nx2,ny2,nz2
  integer :: minx1,maxx1,miny1,maxy1,minz1,maxz1,type1,iselect1
  integer :: minx2,maxx2,miny2,maxy2,minz2,maxz2,type2,iselect2
  
  logical :: ltest(ntest)
  logical :: lxyzfile(ntest)
  logical :: lskip(ntest)
  
  integer :: i,j,k,l,ii,jj,kk,idrank,ncpu,ndec,xdim,ydim,zdim,iii
  
  character(len=256) :: file_loc_proc
  logical :: lmpi=.false.
  character(len=8) :: mystring8,atmname
  logical :: lexist
  character(len=40) :: mystring40
  
  real, parameter :: mytol=5.e-15
  
  integer :: natms1,natms2,iatm,nargxyz1,nargxyz2
  logical :: lrotate1,lrotate2,lerase
  
  
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
  labeltest(14)='3D_Particles_box'
  labeltest(15)='3D_Particle_SC'
  
  lxyzfile(1:10)=.false.
  lxyzfile(11:ntest)=.true.
  !lskip(1:ntest)=.true.
  lskip(1:ntest)=.false.
  !lskip(14)=.false.
  lerase=.false.
  
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
    open(unit=orig_mapio,file=trim(origmap_file),status='old',action='read')
    read(orig_mapio,'(4i10)')mxrank1,nx1,ny1,nz1
    close(orig_mapio)
    
    inquire(file=trim(map_file),exist=lexist)
    if(.not.lexist)then
      ltest(i)=.false.
      cycle
    endif
    open(unit=map_io,file=trim(map_file),status='old',action='read')
    read(map_io,'(4i10)')mxrank2,nx2,ny2,nz2
    close(map_io)
    
    open(unit=orig_io,file=trim(orig_file),status='old',action='read')
    
    read(orig_io,'(i10)')imio1
    
    read(orig_io,'(8i10)')minx1,maxx1,miny1,maxy1,minz1,maxz1,type1,iselect1
    
    !if(any(imio1.ne.imio2)
    
    if(nx1/=nx2 .and. ny1/=ny2 .and. nz1/=nz2)then
      ltest(i)=.false.
    else
      if(iselect1==4)then
        allocate(dmio1(4,nx1,ny1,nz1))
        allocate(dmio2(4,nx2,ny2,nz2))
        
        dmio1(:,:,:,:)=0.e0
        dmio2(:,:,:,:)=0.e0
        
        do l=0,mxrank1-1
          file_loc_proc=repeat(' ',256)
          file_loc_proc = 'output'//trim(write_fmtnumb(l))//'.map.orig'
          open(unit=orig_io+l,file=trim(file_loc_proc),status='old',action='read')
          do kk=minz1,maxz1
            do jj=miny1,maxy1
              do ii=minx1,maxx1
                read(orig_io+l,'(3i8,4f24.16)',end=100,err=100)imio3(1:3),dmio1(1:4,ii,jj,kk)
              enddo
            enddo
          enddo
          close(orig_io+l)
        enddo
        
        do l=0,mxrank2-1
          file_loc_proc=repeat(' ',256)
          file_loc_proc = 'output'//trim(write_fmtnumb(l))//'.map'
          inquire(file=trim(file_loc_proc),exist=lexist)
          if(.not.lexist)goto 100
          inquire(unit=actual_io+l,opened=lexist)
          if(lexist)close(actual_io+l)
          open(unit=actual_io+l,file=trim(file_loc_proc),status='old',action='read')
          read(actual_io+l,'(i10)')imio2
          read(actual_io+l,'(8i10)')minx2,maxx2,miny2,maxy2,minz2,maxz2,type2,iselect2
          if(iselect1.ne.iselect2)then
            ltest(i)=.false.
          endif
          do kk=minz2,maxz2
            do jj=miny2,maxy2
              do ii=minx2,maxx2
                read(actual_io+l,'(3i8,4f24.16)',end=100,err=100)imio4(1:3),dmio2(1:4,ii,jj,kk)
              enddo
            enddo
          enddo
          close(actual_io+l)
        enddo
        
        do kk=1,nz1
          do jj=1,ny1
            do ii=1,nx1
              do l=1,4
                if(abs(dmio1(l,ii,jj,kk)-dmio2(l,ii,jj,kk))>mytol)then
                  ! write (*,*) "i,j,k, l", ii,jj,kk,l, dmio1(l,ii,jj,kk),dmio2(l,ii,jj,kk)
                  ltest(i)=.false.
                  cycle
                endif
              enddo
            enddo
          enddo
        enddo
        
        if(lxyzfile(i))then
          open(unit=orig_ioxyz,file=trim(orig_xyzfile),status='old',action='read')
          read(orig_ioxyz,*)natms1
          mystring40=repeat(' ',40)
          read(orig_ioxyz,'(a)')mystring40
          if(trim(mystring40)=='read list vxx vyy vzz phi theta psi')then
            lrotate1=.true.
            nargxyz1=9
          else
            lrotate1=.false.
            nargxyz1=6
          endif
          if(lrotate1)then
            allocate(dmio3(10,natms1))
            dmio3(:,:)=0.e0
            do iatm=1,natms1
               read(orig_ioxyz,'(a8,10g16.8)')atmname,(dmio3(iii,iatm),iii=1,nargxyz1)
            enddo
          else
            allocate(dmio3(6,natms1))
            dmio3(:,:)=0.e0
            do iatm=1,natms1
               read(orig_ioxyz,'(a8,6g16.8)')atmname,(dmio3(iii,iatm),iii=1,nargxyz1)
            enddo
          endif
          close(orig_ioxyz)
          
          open(unit=actual_ioxyz,file=trim(actual_xyzfile),status='old',action='read')
          read(actual_ioxyz,*)natms2
          mystring40=repeat(' ',40)
          read(actual_ioxyz,'(a)')mystring40
          if(trim(mystring40)=='read list vxx vyy vzz phi theta psi')then
            lrotate2=.true.
            nargxyz2=9
          else
            lrotate2=.false.
            nargxyz2=6
          endif
          if(lrotate2)then
            allocate(dmio4(10,natms2))
            dmio4(:,:)=0.e0
            do iatm=1,natms2
               read(actual_ioxyz,'(a8,10g16.8)')atmname,(dmio4(iii,iatm),iii=1,nargxyz2)
            enddo
          else
            allocate(dmio4(6,natms2))
            dmio4(:,:)=0.e0
            do iatm=1,natms2
               read(actual_ioxyz,'(a8,6g16.8)')atmname,(dmio4(iii,iatm),iii=1,nargxyz2)
            enddo
          endif
          close(actual_ioxyz)
          
          if((natms1/=natms2) .or. (nargxyz1/=nargxyz2) .or. &
           (lrotate1 .neqv. lrotate2))then
            ltest(i)=.false.
          else
            do iatm=1,natms1
              do iii=1,nargxyz1
                if(abs(dmio3(iii,iatm)-dmio4(iii,iatm))>mytol)then
                  ! write (*,*) "iii,iatm", iii,iatm, dmio3(iii,iatm),dmio4(iii,iatm)
                  ltest(i)=.false.
                  cycle
                endif
              enddo
            enddo
          endif
          
          deallocate(dmio3,dmio4)
        endif
        
        open(unit=remove_file_io,file='rm.sh',status='replace')
        if(lmpi)then
          if(lxyzfile(i))then
            write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
             'rm -f output00*.map global.map time.dat mysed.sh restart.xyz rm.sh'
          else
            write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
             'rm -f output00*.map global.map time.dat mysed.sh rm.sh'
          endif
        else
          if(lxyzfile(i))then
            write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
             'rm -f output00*.map global.map time.dat restart.xyz rm.sh'
          else
            write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
             'rm -f output00*.map global.map time.dat rm.sh'
          endif
        endif
        close(remove_file_io)
        call execute_command_line('chmod 777 rm.sh',wait=.true.)
        if(lerase)call execute_command_line('./rm.sh',wait=.true.)
        deallocate(dmio1)
        deallocate(dmio2)
        if(ltest(i))then
          write(6,'(3a)')'test ',labeltest(i),' : OK'
        else
          write(6,'(3a)')'test ',labeltest(i),' : NOK'
        endif
        cycle
 100    ltest(i)=.false.
        close(orig_io)
        close(actual_io)
        open(unit=remove_file_io,file='rm.sh',status='replace')
        if(lmpi)then
          write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
         'rm -f output00*.map global.map time.dat mysed.sh rm.sh'
        else 
          write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
           'rm -f output00*.map global.map time.dat rm.sh'
        endif
        close(remove_file_io)
        call execute_command_line('chmod 777 rm.sh',wait=.true.)
        if(lerase)call execute_command_line('./rm.sh',wait=.true.)
        deallocate(dmio1)
        deallocate(dmio2)
        if(ltest(i))then
          write(6,'(3a)')'test ',labeltest(i),' : OK'
        else
          write(6,'(3a)')'test ',labeltest(i),' : NOK'
        endif
        cycle
      elseif(iselect1==5)then
        
        allocate(dmio1(5,nx1,ny1,nz1))
        allocate(dmio2(5,nx2,ny2,nz2))
        
        dmio1(:,:,:,:)=0.e0
        dmio2(:,:,:,:)=0.e0
        
        do l=0,mxrank1-1
          file_loc_proc=repeat(' ',256)
          file_loc_proc = 'output'//trim(write_fmtnumb(l))//'.map.orig'
          open(unit=orig_io+l,file=trim(file_loc_proc),status='old',action='read')
          do kk=minz1,maxz1
            do jj=miny1,maxy1
              do ii=minx1,maxx1
                read(orig_io+l,'(3i8,5f24.16)',end=110,err=110)imio3(1:3),dmio1(1:5,ii,jj,kk)
              enddo
            enddo
          enddo
          close(orig_io+l)
        enddo
        
        do l=0,mxrank2-1
          file_loc_proc=repeat(' ',256)
          file_loc_proc = 'output'//trim(write_fmtnumb(l))//'.map'
          inquire(file=trim(file_loc_proc),exist=lexist)
          if(.not.lexist)goto 110
          inquire(unit=actual_io+l,opened=lexist)
          if(lexist)close(actual_io+l)
          open(unit=actual_io+l,file=trim(file_loc_proc),status='old',action='read')
          read(actual_io+l,'(i10)')imio2
          read(actual_io+l,'(8i10)')minx2,maxx2,miny2,maxy2,minz2,maxz2,type2,iselect2
          if(iselect1.ne.iselect2)then
            ltest(i)=.false.
          endif
          do kk=minz2,maxz2
            do jj=miny2,maxy2
              do ii=minx2,maxx2
                read(actual_io+l,'(3i8,5f24.16)',end=110,err=110)imio4(1:3),dmio2(1:5,ii,jj,kk)
              enddo
            enddo
          enddo
          close(actual_io+l)
        enddo
        
        do kk=1,nz1
          do jj=1,ny1
            do ii=1,nx1
              do l=1,5
                if(abs(dmio1(l,ii,jj,kk)-dmio2(l,ii,jj,kk))>mytol)then
                  ltest(i)=.false.
                  cycle
                endif
              enddo
            enddo
          enddo
        enddo
        
        if(lxyzfile(i))then
          open(unit=orig_ioxyz,file=trim(orig_xyzfile),status='old',action='read')
          read(orig_ioxyz,*)natms1
          mystring40=repeat(' ',40)
          read(orig_ioxyz,'(a)')mystring40
          if(trim(mystring40)=='read list vxx vyy vzz phi theta psi')then
            lrotate1=.true.
            nargxyz1=9
          else
            lrotate1=.false.
            nargxyz1=6
          endif
          if(lrotate1)then
            allocate(dmio3(10,natms1))
            dmio3(:,:)=0.e0
            do iatm=1,natms1
               read(orig_ioxyz,'(a8,10g16.8)')atmname,(dmio3(iii,iatm),iii=1,nargxyz1)
            enddo
          else
            allocate(dmio3(6,natms1))
            dmio3(:,:)=0.e0
            do iatm=1,natms1
               read(orig_ioxyz,'(a8,6g16.8)')atmname,(dmio3(iii,iatm),iii=1,nargxyz1)
            enddo
          endif
          close(orig_ioxyz)
          
          open(unit=actual_ioxyz,file=trim(actual_xyzfile),status='old',action='read')
          read(actual_ioxyz,*)natms2
          mystring40=repeat(' ',40)
          read(actual_ioxyz,'(a)')mystring40
          if(trim(mystring40)=='read list vxx vyy vzz phi theta psi')then
            lrotate2=.true.
            nargxyz2=9
          else
            lrotate2=.false.
            nargxyz2=6
          endif
          if(lrotate2)then
            allocate(dmio4(10,natms2))
            dmio4(:,:)=0.e0
            do iatm=1,natms2
               read(actual_ioxyz,'(a8,10g16.8)')atmname,(dmio4(iii,iatm),iii=1,nargxyz2)
            enddo
          else
            allocate(dmio4(6,natms2))
            dmio4(:,:)=0.e0
            do iatm=1,natms2
               read(actual_ioxyz,'(a8,6g16.8)')atmname,(dmio4(iii,iatm),iii=1,nargxyz2)
            enddo
          endif
          close(actual_ioxyz)
          
          if((natms1/=natms2) .or. (nargxyz1/=nargxyz2) .or. &
           (lrotate1 .neqv. lrotate2))then
            ltest(i)=.false.
          else
            do iatm=1,natms1
              do iii=1,nargxyz1
                if(abs(dmio3(iii,iatm)-dmio4(iii,iatm))>mytol)then
                  ltest(i)=.false.
                  cycle
                endif
              enddo
            enddo
          endif
          
          deallocate(dmio3,dmio4)
        endif
        
        open(unit=remove_file_io,file='rm.sh',status='replace')
        if(lmpi)then
         if(lxyzfile(i))then
            write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
             'rm -f output00*.map global.map time.dat mysed.sh restart.xyz rm.sh'
          else
            write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
             'rm -f output00*.map global.map time.dat mysed.sh rm.sh'
          endif
        else 
          if(lxyzfile(i))then
            write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
             'rm -f output00*.map global.map time.dat restart.xyz rm.sh'
          else
            write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
             'rm -f output00*.map global.map time.dat rm.sh'
          endif
        endif
        close(remove_file_io)
        call execute_command_line('chmod 777 rm.sh',wait=.true.)
        if(lerase)call execute_command_line('./rm.sh',wait=.true.)
        deallocate(dmio1)
        deallocate(dmio2)
        if(ltest(i))then
          write(6,'(3a)')'test ',labeltest(i),' : OK'
        else
          write(6,'(3a)')'test ',labeltest(i),' : NOK'
        endif
        cycle
 110    ltest(i)=.false.
        close(orig_io)
        close(actual_io)
        open(unit=remove_file_io,file='rm.sh',status='replace')
        if(lmpi)then
          if(lxyzfile(i))then
            write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
             'rm -f output00*.map global.map time.dat mysed.sh restart.xyz rm.sh'
          else
            write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
             'rm -f output00*.map global.map time.dat mysed.sh rm.sh'
          endif
        else 
          if(lxyzfile(i))then
            write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
             'rm -f output00*.map global.map time.dat restart.xyz rm.sh'
          else
            write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
             'rm -f output00*.map global.map time.dat rm.sh'
          endif
        endif
        close(remove_file_io)
        call execute_command_line('chmod 777 rm.sh',wait=.true.)
        if(lerase)call execute_command_line('./rm.sh',wait=.true.)
        deallocate(dmio1)
        deallocate(dmio2)
        if(ltest(i))then
          write(6,'(3a)')'test ',labeltest(i),' : OK'
        else
          write(6,'(3a)')'test ',labeltest(i),' : NOK'
        endif
        cycle
      else
        write(6,*)'error test',iselect1,iselect2
        write(6,*)'error test',nx1, ny1, nz1
        write(6,*)'error test',nx2, ny2, nz2
        stop
      endif
    endif
  enddo
  
  
  write(6,'(/,a,/)')repeat('*',80)
  write(6,'(a)')'list of test results: '
  do i=1,ntest
    if(ltest(i))then
      write(6,'(3a)')'test ',labeltest(i),' : OK'
    else
      write(6,'(3a)')'test ',labeltest(i),' : NOK'
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
 
 
  
