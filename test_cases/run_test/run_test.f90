
 program run_test
 
  implicit none
  
  integer, parameter :: ntest=7
  
  integer, parameter :: orig_io=1017
  integer, parameter :: orig_mapio=18
  integer, parameter :: map_io=19
  integer, parameter :: actual_io=2027
  integer, parameter :: remove_file_io=37
  character(len=*), parameter :: orig_file='output000000.map.orig' 
  character(len=*), parameter :: actual_file='output000000.map'
  character(len=*), parameter :: origmap_file='global.map.orig'
  character(len=*), parameter :: map_file='global.map'
  integer, parameter :: maxlen=120
  
  character(len=maxlen),dimension(ntest) :: labeltest
  character(len=maxlen) :: origpath,newpath,temppath,lbsoftpath
  character(len=1) :: delimiter
  integer :: imio1,imio2
  integer, dimension(3) :: imio3,imio4
  real, dimension(:,:,:,:),allocatable :: dmio1,dmio2
  integer :: mxrank1,mxrank2,nx1,ny1,nz1,nx2,ny2,nz2
  integer :: minx1,maxx1,miny1,maxy1,minz1,maxz1,type1,iselect1
  integer :: minx2,maxx2,miny2,maxy2,minz2,maxz2,type2,iselect2
  
  logical :: ltest(ntest)
  
  integer :: i,j,k,l,ii,jj,kk,idrank
  
  character(len=256) :: file_loc_proc
  
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
    
    open(unit=orig_mapio,file=trim(origmap_file),status='old',action='read')
    read(orig_mapio,'(4i10)')mxrank1,nx1,ny1,nz1
    close(orig_mapio)
    
    open(unit=map_io,file=trim(map_file),status='old',action='read')
    read(map_io,'(4i10)')mxrank2,nx2,ny2,nz2
    close(map_io)
    
    open(unit=orig_io,file=trim(orig_file),status='old',action='read')
    open(unit=actual_io,file=trim(actual_file),status='old',action='read')
    read(orig_io,'(i10)')imio1
    read(actual_io,'(i10)')imio2
    read(orig_io,'(8i10)')minx1,maxx1,miny1,maxy1,minz1,maxz1,type1,iselect1
    read(actual_io,'(8i10)')minx2,maxx2,miny2,maxy2,minz2,maxz2,type2,iselect2
    !if(any(imio1.ne.imio2)
    if(iselect1.ne.iselect2)then
      ltest(i)=.false.
    elseif(nx1/=nx2 .and. ny1/=ny2 .and. nz1/=nz2)then
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
                read(orig_io+l,'(3i8,4f20.10)',end=100,err=100)imio3(1:3),dmio1(1:4,ii,jj,kk)
              enddo
            enddo
          enddo
          close(orig_io+l)
        enddo
        
        do l=0,mxrank2-1
          file_loc_proc=repeat(' ',256)
          file_loc_proc = 'output'//trim(write_fmtnumb(l))//'.map'
          open(unit=actual_io+l,file=trim(file_loc_proc),status='old',action='read')
          do kk=minz2,maxz2
            do jj=miny2,maxy2
              do ii=minx2,maxx2
                read(actual_io+l,'(3i8,4f20.10)',end=100,err=100)imio4(1:3),dmio2(1:4,ii,jj,kk)
              enddo
            enddo
          enddo
          close(actual_io+l)
        enddo
        
        do kk=1,nz1
          do jj=1,ny1
            do ii=1,nx1
              do l=1,4
                if(dmio1(l,ii,jj,kk) .ne. dmio2(l,ii,jj,kk))then
                  ltest(i)=.false.
                  cycle
                endif
              enddo
            enddo
          enddo
        enddo
        
        open(unit=remove_file_io,file='rm.sh',status='replace')
        write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
         'rm -f output00*.map global.map time.dat rm.sh'
        close(remove_file_io)
        call execute_command_line('chmod 777 rm.sh',wait=.true.)
        call execute_command_line('./rm.sh',wait=.true.)
        deallocate(dmio1)
        deallocate(dmio2)
        cycle
 100    ltest(i)=.false.
        close(orig_io)
        close(actual_io)
        open(unit=remove_file_io,file='rm.sh',status='replace')
        write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
         'rm -f output00*.map global.map time.dat rm.sh'
        close(remove_file_io)
        call execute_command_line('chmod 777 rm.sh',wait=.true.)
        call execute_command_line('./rm.sh',wait=.true.)
        deallocate(dmio1)
        deallocate(dmio2)
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
                read(orig_io+l,'(3i8,5f20.10)',end=110,err=110)imio3(1:3),dmio1(1:5,ii,jj,kk)
              enddo
            enddo
          enddo
          close(orig_io+l)
        enddo
        
        do l=0,mxrank2-1
          file_loc_proc=repeat(' ',256)
          file_loc_proc = 'output'//trim(write_fmtnumb(l))//'.map'
          open(unit=actual_io+l,file=trim(file_loc_proc),status='old',action='read')
          do kk=minz2,maxz2
            do jj=miny2,maxy2
              do ii=minx2,maxx2
                read(actual_io+l,'(3i8,5f20.10)',end=110,err=110)imio4(1:3),dmio2(1:5,ii,jj,kk)
              enddo
            enddo
          enddo
          close(actual_io+l)
        enddo
        
        do kk=1,nz1
          do jj=1,ny1
            do ii=1,nx1
              do l=1,5
                if(dmio1(l,ii,jj,kk) .ne. dmio2(l,ii,jj,kk))then
                  ltest(i)=.false.
                  cycle
                endif
              enddo
            enddo
          enddo
        enddo
        
        open(unit=remove_file_io,file='rm.sh',status='replace')
        write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
         'rm -f output00*.map global.map time.dat rm.sh'
        close(remove_file_io)
        call execute_command_line('chmod 777 rm.sh',wait=.true.)
        call execute_command_line('./rm.sh',wait=.true.)
        deallocate(dmio1)
        deallocate(dmio2)
        cycle
 110    ltest(i)=.false.
        close(orig_io)
        close(actual_io)
        open(unit=remove_file_io,file='rm.sh',status='replace')
        write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
         'rm -f output00*.map global.map time.dat rm.sh'
        close(remove_file_io)
        call execute_command_line('chmod 777 rm.sh',wait=.true.)
        call execute_command_line('./rm.sh',wait=.true.)
        deallocate(dmio1)
        deallocate(dmio2)
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
      write(6,'(3a)')'test ',labeltest(i),' : NO'
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
  
 end program run_test
 
 
  
