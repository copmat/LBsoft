
 program compare
 
  implicit none
  
  
  
  integer, parameter :: orig_io=1017
  integer, parameter :: orig_mapio=18
  integer, parameter :: map_io=19
  integer, parameter :: actual_io=2027
  integer, parameter :: remove_file_io=37
  character(len=*), parameter :: orig_file='output000000.map.orig' 
  
  character(len=*), parameter :: origmap_file='global.map.orig'
  character(len=*), parameter :: map_file='global.map'
  integer, parameter :: maxlen=120
  
  
  character(len=maxlen) :: origpath,newpath,temppath,lbsoftpath
  character(len=1) :: delimiter
  integer :: imio1,imio2
  integer, dimension(3) :: imio3,imio4
  real, dimension(:,:,:,:),allocatable :: dmio1,dmio2
  integer :: mxrank1,mxrank2,nx1,ny1,nz1,nx2,ny2,nz2
  integer :: minx1,maxx1,miny1,maxy1,minz1,maxz1,type1,iselect1
  integer :: minx2,maxx2,miny2,maxy2,minz2,maxz2,type2,iselect2
  
  integer :: itest
  
  integer :: i,j,k,l,ii,jj,kk,idrank
  
  character(len=256) :: file_loc_proc
  
  real, parameter :: mytol=1.e-15
  
  
  origpath=repeat(' ',maxlen)
  
  
  call getcwd(origpath)
  origpath=trim(origpath)
  delimiter=origpath(1:1)
  
  
  itest=1

  write(6,'(a)')'read original output '

    
  open(unit=orig_mapio,file=trim(origmap_file),status='old',action='read')
  read(orig_mapio,'(4i10)')mxrank1,nx1,ny1,nz1
  close(orig_mapio)
    
  open(unit=map_io,file=trim(map_file),status='old',action='read')
  read(map_io,'(4i10)')mxrank2,nx2,ny2,nz2
  close(map_io)
    
  open(unit=orig_io,file=trim(orig_file),status='old',action='read')
  
  read(orig_io,'(i10)')imio1
  
  read(orig_io,'(8i10)')minx1,maxx1,miny1,maxy1,minz1,maxz1,type1,iselect1
  
  
  if(nx1/=nx2 .and. ny1/=ny2 .and. nz1/=nz2)then
    itest=0
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
              if(imio3(1)/=ii .or. imio3(2)/=jj .or. imio3(3)/=kk)then
                write(6,*)'problem in reading',trim(file_loc_proc),imio3(1:3),ii,jj,kk
                stop
              endif
            enddo
          enddo
        enddo
        close(orig_io+l)
      enddo
        
      do l=0,mxrank2-1
        file_loc_proc=repeat(' ',256)
        file_loc_proc = 'output'//trim(write_fmtnumb(l))//'.map'
        open(unit=actual_io+l,file=trim(file_loc_proc),status='old',action='read')
        read(actual_io+l,'(i10)')imio2
        read(actual_io+l,'(8i10)')minx2,maxx2,miny2,maxy2,minz2,maxz2,type2,iselect2
        if(iselect1.ne.iselect2)then
          itest=0
          goto 130
        endif
        do kk=minz2,maxz2
          do jj=miny2,maxy2
            do ii=minx2,maxx2
              read(actual_io+l,'(3i8,4f24.16)',end=100,err=100)imio4(1:3),dmio2(1:4,ii,jj,kk)
              if(imio4(1)/=ii .or. imio4(2)/=jj .or. imio4(3)/=kk)then
                write(6,*)'problem in reading',trim(file_loc_proc),imio4(1:3),ii,jj,kk
                stop
              endif
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
                itest=0
                write(6,'(2g20.10,4i6)')dmio1(l,ii,jj,kk),dmio2(l,ii,jj,kk),ii,jj,kk,l
                goto 130
              endif
            enddo
          enddo
        enddo
      enddo
        close(orig_io)
        close(actual_io)
!        open(unit=remove_file_io,file='rm.sh',status='replace')
!        write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
!         'rm -f output00*.map global.map time.dat rm.sh'
!        close(remove_file_io)
!        call execute_command_line('chmod 777 rm.sh',wait=.true.)
!        call execute_command_line('./rm.sh',wait=.true.)
        deallocate(dmio1)
        deallocate(dmio2)
        goto 130
 100    itest=0
        close(orig_io)
        close(actual_io)
!        open(unit=remove_file_io,file='rm.sh',status='replace')
!        write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
!         'rm -f output00*.map global.map time.dat rm.sh'
!        close(remove_file_io)
!        call execute_command_line('chmod 777 rm.sh',wait=.true.)
!        call execute_command_line('./rm.sh',wait=.true.)
        deallocate(dmio1)
        deallocate(dmio2)
        goto 130
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
                if(imio3(1)/=ii .or. imio3(2)/=jj .or. imio3(3)/=kk)then
                  write(6,*)'problem in reading',trim(file_loc_proc),imio3(1:3),ii,jj,kk
                  stop
                endif
              enddo
            enddo
          enddo
          close(orig_io+l)
        enddo
        
        do l=0,mxrank2-1
          file_loc_proc=repeat(' ',256)
          file_loc_proc = 'output'//trim(write_fmtnumb(l))//'.map'
          open(unit=actual_io+l,file=trim(file_loc_proc),status='old',action='read')
          read(actual_io+l,'(i10)')imio2
          read(actual_io+l,'(8i10)')minx2,maxx2,miny2,maxy2,minz2,maxz2,type2,iselect2
          if(iselect1.ne.iselect2)then
            itest=0
            goto 130
          endif
          do kk=minz2,maxz2
            do jj=miny2,maxy2
              do ii=minx2,maxx2
                read(actual_io+l,'(3i8,5f24.16)',end=110,err=110)imio4(1:3),dmio2(1:5,ii,jj,kk)
                if(imio4(1)/=ii .or. imio4(2)/=jj .or. imio4(3)/=kk)then
                  write(6,*)'problem in reading',trim(file_loc_proc),imio4(1:3),ii,jj,kk
                  stop
                endif
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
                  itest=0
                  write(6,'(2g20.10,4i6)')dmio1(l,ii,jj,kk),dmio2(l,ii,jj,kk),ii,jj,kk,l
                  goto 130
                endif
              enddo
            enddo
          enddo
        enddo
        
!        open(unit=remove_file_io,file='rm.sh',status='replace')
!        write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
!         'rm -f output00*.map global.map time.dat rm.sh'
!        close(remove_file_io)
!        call execute_command_line('chmod 777 rm.sh',wait=.true.)
!        call execute_command_line('./rm.sh',wait=.true.)
        close(orig_io)
        close(actual_io)
        deallocate(dmio1)
        deallocate(dmio2)
        goto 130
 110    itest=0
        close(orig_io)
        close(actual_io)
!        open(unit=remove_file_io,file='rm.sh',status='replace')
!        write(remove_file_io,'(a,/,a)')'#!/bin/bash', &
!         'rm -f output00*.map global.map time.dat rm.sh'
!        close(remove_file_io)
!        call execute_command_line('chmod 777 rm.sh',wait=.true.)
!        call execute_command_line('./rm.sh',wait=.true.)
        deallocate(dmio1)
        deallocate(dmio2)
        goto 130
      else
        write(6,*)'error test',iselect1,iselect2
        write(6,*)'error test',nx1, ny1, nz1
        write(6,*)'error test',nx2, ny2, nz2
        stop
      endif
    endif
 
  
130 continue
  
  write(6,'(/,a,/)')repeat('*',80)
  write(6,'(a)')'test results: '

    if(itest==1)then
      write(6,'(3a)')'test : OK'
    else
      write(6,'(3a)')'test : NO'
    endif
    
  open(unit=remove_file_io,file='test.out',status='replace')
  write(remove_file_io,*)itest
  close(remove_file_io)

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
  
 end program compare
 
 
  
