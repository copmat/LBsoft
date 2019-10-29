
 program prepare
 
  implicit none
  
  integer :: nx,ny,nz,nn,i,j,k,cx,cy,cz
  
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
   
  do k=1,nz
    do j=1,ny
      do i=1,nx
        rdist=dsqrt(dble((i-cx)**2+(j-cy)**2+(k-cz)**2))
        if(rdist<=radius)then
          isfluid(i,j,k)=3
        else
          isfluid(i,j,k)=1
        endif
        write(iout)isfluid(i,j,k)
      enddo
    enddo
  enddo
  
  close(iout)
  
  stop
  
 end program prepare
        
        
        
