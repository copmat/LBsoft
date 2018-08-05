
 program run_test
 
  implicit none
  
  integer, parameter :: ntest=7
  
  integer, parameter :: orig_io=17
  integer, parameter :: actual_io=27
  integer, parameter :: remove_file_io=37
  character(len=*), parameter :: orig_file='output000000.map.orig' 
  character(len=*), parameter :: actual_file='output000000.map'
  integer, parameter :: maxlen=120
  
  character(len=maxlen),dimension(ntest) :: labeltest
  character(len=maxlen) :: origpath,newpath,temppath,lbsoftpath
  character(len=1) :: delimiter
  integer, dimension(6) :: imio1,imio2
  integer, dimension(3) :: imio3,imio4
  real, dimension(5) :: dmio1,dmio2
  
  logical :: ltest(ntest)
  
  integer :: i,ii,jj,kk
  
  labeltest(1:ntest)=repeat(' ',maxlen)
  origpath=repeat(' ',maxlen)
  
  labeltest(1)='2D_Poiseuille_xy'
  labeltest(2)='2D_Poiseuille_xz'
  labeltest(3)='2D_Poiseuille_yz'
  labeltest(4)='2D_Shear_xy'
  labeltest(5)='2D_Shear_xz'
  labeltest(6)='2D_Shear_yz'
  labeltest(7)='3D_Spinodal'
  
  call getcwd(origpath)
  origpath=trim(origpath)
  delimiter=origpath(1:1)
  
  lbsoftpath=repeat(' ',maxlen)
  
  lbsoftpath='..'//delimiter//'..'//delimiter//'execute'//delimiter//'main.x'
  !write(6,*)'sono in: ',trim(lbsoftpath)
  
  ltest(1:ntest)=.true.
  do i=1,ntest
    call chdir(trim(origpath))
    newpath=repeat(' ',maxlen)
    newpath='..'//delimiter//labeltest(i)
    !write(6,*)'entro in: ',trim(newpath)
    call chdir(trim(newpath))
    temppath=repeat(' ',maxlen)
    call getcwd(temppath)
    !write(6,*)'sono in: ',trim(temppath)
    write(6,'(2a)')'execute test: ',labeltest(i)
    call execute_command_line(trim(lbsoftpath),wait=.true.)
    write(6,'(2a)')'check test  : ',labeltest(i)
    open(unit=orig_io,file=trim(orig_file),status='old',action='read')
    open(unit=actual_io,file=trim(actual_file),status='old',action='read')
    read(orig_io,'(i10)')imio1(1)
    read(actual_io,'(i10)')imio2(1)
    read(orig_io,'(6i10)')imio1(1:6)
    read(actual_io,'(6i10)')imio2(1:6)
    !if(any(imio1.ne.imio2)
    if(any(imio1.ne.imio2))then
      ltest(i)=.false.
    else
      if(imio1(5)==4)then
        dmio1(:)=0.e0
        dmio2(:)=0.e0
        do kk=1,imio1(3)
          do jj=1,imio1(2)
            do ii=1,imio1(1)
              read(orig_io,'(3i8,4f20.10)',end=100,err=100)imio3(1:3),dmio1(1:4)
              read(actual_io,'(3i8,4f20.10)',end=100,err=100)imio4(1:3),dmio2(1:4)
              if(any(imio3.ne.imio4) .or. any(dmio1.ne.dmio2))then
                ltest(i)=.false.
                cycle
              endif
            enddo
          enddo
        enddo
        close(orig_io)
        close(actual_io)
        open(unit=remove_file_io,file='rm.sh',status='replace')
        write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
         'rm -f output000000.map time.dat rm.sh'
        close(remove_file_io)
        call execute_command_line('chmod 777 rm.sh',wait=.true.)
        call execute_command_line('./rm.sh',wait=.true.)
        cycle
 100    ltest(i)=.false.
        close(orig_io)
        close(actual_io)
        open(unit=remove_file_io,file='rm.sh',status='replace')
        write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
         'rm -f output000000.map time.dat rm.sh'
        close(remove_file_io)
        call execute_command_line('chmod 777 rm.sh',wait=.true.)
        call execute_command_line('./rm.sh',wait=.true.)
        cycle
      elseif(imio1(5)==5)then
        dmio1(:)=0.e0
        dmio2(:)=0.e0
        do kk=1,imio1(3)
          do jj=1,imio1(2)
            do ii=1,imio1(1)
              read(orig_io,'(3i8,5f20.10)',end=110,err=110)imio3(1:3),dmio1(1:5)
              read(actual_io,'(3i8,5f20.10)',end=110,err=110)imio4(1:3),dmio2(1:5)
              if(any(imio3.ne.imio4) .or. any(dmio1.ne.dmio2))then
                ltest(i)=.false.
                cycle
              endif
            enddo
          enddo
        enddo
        close(orig_io)
        close(actual_io)
        open(unit=remove_file_io,file='rm.sh',status='replace')
        write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
         'rm -f output000000.map time.dat rm.sh'
        close(remove_file_io)
        call execute_command_line('chmod 777 rm.sh',wait=.true.)
        call execute_command_line('./rm.sh',wait=.true.)
        cycle
 110    ltest(i)=.false.
        close(orig_io)
        close(actual_io)
        open(unit=remove_file_io,file='rm.sh',status='replace')
        write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
         'rm -f output000000.map time.dat rm.sh'
        close(remove_file_io)
        call execute_command_line('chmod 777 rm.sh',wait=.true.)
        call execute_command_line('./rm.sh',wait=.true.)
        cycle
      else
        write(6,*)'error test',imio1(1:6)
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
      write(6,'(3a)')'test ',labeltest(i),' : NO'
    endif
  enddo
  write(6,'(/,a,/)')repeat('*',80)
  stop
  
 end program run_test
 
 
  
