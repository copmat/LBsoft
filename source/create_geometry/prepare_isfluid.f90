
 program prepare
 
!***********************************************************************
!     
!     LBsoft tool for preparing the isfluid.dat file as a set
!     of 1 byte integers in natural order
!     
!     licensed under the 3-Clause BSD License (BSD-3-Clause)
!     author: M. Lauricella
!     last modification October 2019
!     
!***********************************************************************

  implicit none
  
  integer :: nx,ny,nz,nn,i,j,k,cx,cy,cz
  !isfluid should be 1 byte integer
  integer(kind=1), allocatable, dimension(:,:,:) :: isfluid
  
  integer, parameter :: iout=16
  double precision :: rdist,radius
  
  nx=64
  ny=64
  nz=64
  
  cx=32
  cy=32
  cz=32
  
  radius=12
  
  nn=nx*ny*nz
  
  allocate(isfluid(nx,ny,nz))
  
  open(unit=iout,file='isfluid.dat',action='write',status='replace', &
   access='stream',form='unformatted')
  
  !i should be the fastest (natural order of the 3d matrix)
  !remember that FORTRAN is column-major order
  do k=1,nz
    do j=1,ny
      do i=1,nx
        !insert a geometrical rule here
        rdist=dsqrt(dble((i-cx)**2+(j-cy)**2+(k-cz)**2))
        if(rdist<=radius)then
          isfluid(i,j,k)=3 !3 = solid
        else
          isfluid(i,j,k)=1 !1 = liquid
        endif
        !print on isfluid.dat
        write(iout)isfluid(i,j,k)
      enddo
    enddo
  enddo
  
  close(iout)
  
  stop
  
 end program prepare
        
        
        
