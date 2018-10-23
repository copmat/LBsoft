
 program create_xyz
 
  implicit none
  
  integer :: nx=64
  integer :: ny=64
  integer :: nz=64
  
  integer :: ntot,i,j,k
  
  integer, parameter :: ioxyz=16
  character(len=*), parameter :: filexyz='input.xyz'
  
  ntot=nx*ny*nz/8
  
  open(unit=ioxyz,file=filexyz,status='replace',action='write')
  
  write(ioxyz,'(i8)')ntot
  
  write(ioxyz,'(a)')' read list vx vy vz '
  
  do k=1,nz,2
    do j=1,ny,2
      do i=1,nx,2
        write(ioxyz,'(a8,6f16.8)')'C       ',real(i),real(j),real(k), &
         real(i)*0.11d0,real(j)*0.12d0,real(k)*0.13d0
      enddo
    enddo
  enddo
  
  close(ioxyz)
  
  stop
  
 end program create_xyz
 
 
  
