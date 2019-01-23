#include <default_macro.h>
module mpi_comm
    use mpi
    use aop_mod
    use lbempi_mod,  only : domdec

    implicit none

    integer            :: ierr                                           ! error code
    integer            :: id_rank                                        ! rank of process
    integer            :: mx_rank                                        ! number of processes

    integer, parameter  :: npop=18
    !lattice vectors
    integer, dimension(0:npop), parameter :: &
    ex = (/ 0, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1, 0, 0, 0, 0/)
    !       0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18
    integer, dimension(0:npop), parameter :: &
    ey = (/ 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1/)
    integer, dimension(0:npop), parameter :: &
    ez = (/ 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1/)

private
    integer            :: globalDims(3), ldims(3), ldimxy,ldimxz,ldimyz
    integer             :: nx,ny,nz
    logical, dimension(3) :: period
    integer            :: mymin(3), mymax(3), minx,maxx, miny,maxy, minz,maxz
    integer            :: cartdims(3), cube_comm, mycoords(3), neigh(0:18) ! IDs in cart grid

    real(kind=PRC),allocatable, dimension(:,:)    :: buf_sendrecv1,buf_sendrecv2,buf_sendrecv3,  &
                            buf_sendrecv4,buf_sendrecv5,buf_sendrecv6,  &
                            buf_sendrecv7,buf_sendrecv8,buf_sendrecv9,  &
                            buf_sendrecv10,buf_sendrecv11,buf_sendrecv12,  &
                            buf_sendrecv13,buf_sendrecv14,buf_sendrecv15,  &
                            buf_sendrecv16,buf_sendrecv17,buf_sendrecv18
    real(kind=PRC),allocatable, dimension(:,:)    :: buf_sendrecv1_hvar,buf_sendrecv2_hvar,buf_sendrecv3_hvar,  &
                            buf_sendrecv4_hvar,buf_sendrecv5_hvar,buf_sendrecv6_hvar,  &
                            buf_sendrecv7_hvar,buf_sendrecv8_hvar,buf_sendrecv9_hvar,  &
                            buf_sendrecv10_hvar,buf_sendrecv11_hvar,buf_sendrecv12_hvar,  &
                            buf_sendrecv13_hvar,buf_sendrecv14_hvar,buf_sendrecv15_hvar,  &
                            buf_sendrecv16_hvar,buf_sendrecv17_hvar,buf_sendrecv18_hvar
    INTEGER(kind=1),allocatable, dimension(:,:)    :: buf_sendrecv1_isfluid,buf_sendrecv2_isfluid,buf_sendrecv3_isfluid,  &
                        buf_sendrecv4_isfluid,buf_sendrecv5_isfluid,buf_sendrecv6_isfluid,  &
                        buf_sendrecv7_isfluid,buf_sendrecv8_isfluid,buf_sendrecv9_isfluid,  &
                        buf_sendrecv10_isfluid,buf_sendrecv11_isfluid,buf_sendrecv12_isfluid,  &
                        buf_sendrecv13_isfluid,buf_sendrecv14_isfluid,buf_sendrecv15_isfluid,  &
                        buf_sendrecv16_isfluid,buf_sendrecv17_isfluid,buf_sendrecv18_isfluid

    integer            :: recv_req(18),send_req(18)
    integer            :: recv_req_hvar(18),send_req_hvar(18)
    integer            :: recv_req_isfluid(18),send_req_isfluid(18)
    integer, parameter :: movetag=300000
    integer, parameter :: hvartag=500000
    integer, parameter :: fldtag =700000
    integer, parameter :: halotag =900000

    public :: mpiInit
    public :: mpisendpops, mpirecvpops, mpibounceback
    public :: mpisend_hvar, mpirecv_hvar
    public :: mpisend_isfluid, mpirecv_isfluid
    public :: mpisendrecvhalopops

contains

    subroutine my_Cart_coords(me, mycoords, ierr)
    implicit none
    integer, intent(in) :: me
    integer, intent(out) :: mycoords(3), ierr
    integer i,j,k

    k = me / (cartdims(1)*cartdims(2))
    j = (me - k*cartdims(1)*cartdims(2)) / cartdims(1)
    i = mod(me, cartdims(1))

    mycoords = [ i,j,k ]

    end subroutine my_Cart_coords



    subroutine my_CART_RANK(mycoords, shiftpos, period, neigh, ierr)
    implicit none
    integer, intent(in) :: mycoords(3),shiftpos(3)
    logical, intent(in) :: period(3)
    integer, intent(out) :: neigh, ierr
    integer i,c, offs(3)

    offs = [ 1, cartdims(1), cartdims(1)*cartdims(2) ]

    neigh = 0
    do c=1,3
        i = shiftpos(c)
        if (i<0) then
            if (period(c)) then
                i = i + cartdims(c)
            else
		neigh = MPI_PROC_NULL
                return
            endif
        endif
        if (i>=cartdims(c)) then
            if (period(c)) then
                i = i - cartdims(c)
            else
		neigh = MPI_PROC_NULL
                return
            endif
        endif

        neigh = neigh + i*offs(c)
    enddo

    end subroutine my_CART_RANK


    subroutine mpiInit(lnx,lny,lnz, ibctype)

    implicit none
    integer, intent(in) :: lnx,lny,lnz, ibctype

    integer         :: i, mycoords(3), shift(3), tempDim, iswap
    integer         :: locindex
    logical         :: reorder



    globalDims(1) = lnx; globalDims(2) = lny; globalDims(3) = lnz;
    nx = lnx; ny = lny; nz = lnz;


    ! find out number of processes
    call MPI_Comm_size(MPI_COMM_WORLD, mx_rank, ierr)

    select case(domdec)
    case(1)
        cartdims(:) = [ 1, 1, 0]
    case(2)
        cartdims(:) = [ 1, 0, 1]
    case(3)
        cartdims(:) = [ 0, 1, 1]
    case(4)
        cartdims(:) = [ 1, 0, 0]
    case(5)
        cartdims(:) = [ 0, 1, 0]
    case(6)
        cartdims(:) = [ 0, 0, 1]
    case default
        cartdims(:) = [ 0, 0, 0]
    end select

    period(:) = .false.
    select case(ibctype)
    case(1)
       period(1) = .true.
    case(2)
       period(2) = .true.
    case(3)
       period(1:2) = .true.
    case(4)
       period(3) = .true.
    case(5)
       period(1) = .true.
       period(3) = .true.
    case(6)
       period(2:3) = .true.
    case(7)
       period(:) = .true.
    end select

    reorder = .false.


    call MPI_DIMS_CREATE(mx_rank, 3, cartdims, ierr)

    ! call MPI_CART_CREATE(MPI_COMM_WORLD, 3, cartdims, period, reorder, cube_comm, ierr)

    ! find out process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, id_rank, ierr)



    if(id_rank == 0)then
       write(*,*) "Global  Domain Dims: ", globalDims
       write(*,*) "Process Mesh   Dims: ", cartdims
    end if

    ! call MPI_Cart_coords(cube_comm, id_rank, 3, mycoords, ierr)
    call my_Cart_coords(id_rank, mycoords, ierr)
    cube_comm = MPI_COMM_WORLD

    do i=0, npop
        shift(1) = ex(i) + mycoords(1)
        shift(2) = ey(i) + mycoords(2)
        shift(3) = ez(i) + mycoords(3)
        !CALL MPI_CART_RANK(cube_comm, shift, neigh(i), ierr)
        CALL my_CART_RANK(mycoords, shift, period, neigh(i), ierr)
    enddo


    do i=1,3
        tempDim = globalDims(i) / cartdims(i)
        ldims(i) = tempDim
        if (tempDim * cartdims(i) /= globalDims(i) ) then
            write (6,*) "globalDims(", i, ")=", globalDims(i), "not divisible by", cartdims(i)
            call MPI_ABORT(cube_comm, 30, ierr)
        endif

        mymin(i) = tempDim * mycoords(i) + 1
        mymax(i) = tempDim * (mycoords(i) + 1)
    enddo

    do i=0, mx_rank-1
        if (i == id_rank) then
            write (6,'("id=", I3, " @(", 3I2,") arr=", 38I3 )') id_rank, mycoords, neigh
            write (6, *) "min max=", mymin(1),mymax(1), mymin(2),mymax(2), mymin(3),mymax(3)
        endif
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    enddo

    minx = mymin(1)
    miny = mymin(2)
    minz = mymin(3)
    maxx = mymax(1)
    maxy = mymax(2)
    maxz = mymax(3)


    ldimyz = ldims(2) * ldims(3)
    ldimxz = ldims(1) * ldims(3)
    ldimxy = ldims(1) * ldims(2)
    allocate( buf_sendrecv1(ldimyz * 19, 2 ) )     !! Parallel planes
    allocate( buf_sendrecv2(ldimyz * 19, 2 ) )
    allocate( buf_sendrecv3(ldimxz * 19, 2 ) )
    allocate( buf_sendrecv4(ldimxz * 19, 2 ) )
    allocate( buf_sendrecv5(ldimxy * 19, 2 ) )
    allocate( buf_sendrecv6(ldimxy * 19, 2 ) )
    allocate( buf_sendrecv7 (ldims(3), 2 ) )                    !! Diag lines forall z
    allocate( buf_sendrecv8 (ldims(3), 2 ) )
    allocate( buf_sendrecv9 (ldims(3), 2 ) )
    allocate( buf_sendrecv10(ldims(3), 2 ) )
    allocate( buf_sendrecv11(ldims(2), 2 ) )                  !! Diag lines forall y
    allocate( buf_sendrecv12(ldims(2), 2 ) )
    allocate( buf_sendrecv13(ldims(2), 2 ) )
    allocate( buf_sendrecv14(ldims(2), 2 ) )
    allocate( buf_sendrecv15(ldims(1), 2 ) )                    !! Diag lines forall x
    allocate( buf_sendrecv16(ldims(1), 2 ) )
    allocate( buf_sendrecv17(ldims(1), 2 ) )
    allocate( buf_sendrecv18(ldims(1), 2 ) )

    allocate( buf_sendrecv1_hvar(ldimyz, 2 ) )     !! Parallel planes
    allocate( buf_sendrecv2_hvar(ldimyz, 2 ) )
    allocate( buf_sendrecv3_hvar(ldimxz, 2 ) )
    allocate( buf_sendrecv4_hvar(ldimxz, 2 ) )
    allocate( buf_sendrecv5_hvar(ldimxy, 2 ) )
    allocate( buf_sendrecv6_hvar(ldimxy, 2 ) )
    allocate( buf_sendrecv7_hvar (ldims(3), 2 ) )                    !! Diag lines forall z
    allocate( buf_sendrecv8_hvar (ldims(3), 2 ) )
    allocate( buf_sendrecv9_hvar (ldims(3), 2 ) )
    allocate( buf_sendrecv10_hvar(ldims(3), 2 ) )
    allocate( buf_sendrecv11_hvar(ldims(2), 2 ) )                  !! Diag lines forall y
    allocate( buf_sendrecv12_hvar(ldims(2), 2 ) )
    allocate( buf_sendrecv13_hvar(ldims(2), 2 ) )
    allocate( buf_sendrecv14_hvar(ldims(2), 2 ) )
    allocate( buf_sendrecv15_hvar(ldims(1), 2 ) )                    !! Diag lines forall x
    allocate( buf_sendrecv16_hvar(ldims(1), 2 ) )
    allocate( buf_sendrecv17_hvar(ldims(1), 2 ) )
    allocate( buf_sendrecv18_hvar(ldims(1), 2 ) )

    allocate( buf_sendrecv1_isfluid(ldimyz, 2 ) )     !! Parallel planes
    allocate( buf_sendrecv2_isfluid(ldimyz, 2 ) )
    allocate( buf_sendrecv3_isfluid(ldimxz, 2 ) )
    allocate( buf_sendrecv4_isfluid(ldimxz, 2 ) )
    allocate( buf_sendrecv5_isfluid(ldimxy, 2 ) )
    allocate( buf_sendrecv6_isfluid(ldimxy, 2 ) )
    allocate( buf_sendrecv7_isfluid (ldims(3), 2 ) )                    !! Diag lines forall z
    allocate( buf_sendrecv8_isfluid (ldims(3), 2 ) )
    allocate( buf_sendrecv9_isfluid (ldims(3), 2 ) )
    allocate( buf_sendrecv10_isfluid(ldims(3), 2 ) )
    allocate( buf_sendrecv11_isfluid(ldims(2), 2 ) )                  !! Diag lines forall y
    allocate( buf_sendrecv12_isfluid(ldims(2), 2 ) )
    allocate( buf_sendrecv13_isfluid(ldims(2), 2 ) )
    allocate( buf_sendrecv14_isfluid(ldims(2), 2 ) )
    allocate( buf_sendrecv15_isfluid(ldims(1), 2 ) )                    !! Diag lines forall x
    allocate( buf_sendrecv16_isfluid(ldims(1), 2 ) )
    allocate( buf_sendrecv17_isfluid(ldims(1), 2 ) )
    allocate( buf_sendrecv18_isfluid(ldims(1), 2 ) )
    end subroutine mpiInit


    subroutine mpipostrecv
    implicit none

    ! Post MPI_IRECV
    call MPI_IRECV( buf_sendrecv1(1,2),  5*ldimyz, MPI_DOUBLE, neigh(2),  neigh(2)+movetag, cube_comm, recv_req(1),  ierr)
    call MPI_IRECV( buf_sendrecv2(1,2),  5*ldimyz, MPI_DOUBLE, neigh(1),  neigh(1)+movetag, cube_comm, recv_req(2),  ierr)
    call MPI_IRECV( buf_sendrecv3(1,2),  5*ldimxz, MPI_DOUBLE, neigh(4),  neigh(4)+movetag, cube_comm, recv_req(3),  ierr)
    call MPI_IRECV( buf_sendrecv4(1,2),  5*ldimxz, MPI_DOUBLE, neigh(3),  neigh(3)+movetag, cube_comm, recv_req(4),  ierr)
    call MPI_IRECV( buf_sendrecv5(1,2),  5*ldimxy, MPI_DOUBLE, neigh(6),  neigh(6)+movetag, cube_comm, recv_req(5),  ierr)
    call MPI_IRECV( buf_sendrecv6(1,2),  5*ldimxy, MPI_DOUBLE, neigh(5),  neigh(5)+movetag, cube_comm, recv_req(6),  ierr)

    call MPI_IRECV( buf_sendrecv7(1,2),  ldims(3), MPI_DOUBLE, neigh(8),  neigh(8)+movetag, cube_comm, recv_req(7),  ierr)
    call MPI_IRECV( buf_sendrecv8(1,2),  ldims(3), MPI_DOUBLE, neigh(7),  neigh(7)+movetag, cube_comm, recv_req(8),  ierr)
    call MPI_IRECV( buf_sendrecv9(1,2),  ldims(3), MPI_DOUBLE, neigh(10), neigh(10)+movetag, cube_comm, recv_req(9), ierr)
    call MPI_IRECV( buf_sendrecv10(1,2), ldims(3), MPI_DOUBLE, neigh(9),  neigh(9)+movetag, cube_comm, recv_req(10), ierr)
    call MPI_IRECV( buf_sendrecv11(1,2), ldims(2), MPI_DOUBLE, neigh(12), neigh(12)+movetag, cube_comm, recv_req(11), ierr)
    call MPI_IRECV( buf_sendrecv12(1,2), ldims(2), MPI_DOUBLE, neigh(11), neigh(11)+movetag, cube_comm, recv_req(12), ierr)
    call MPI_IRECV( buf_sendrecv13(1,2), ldims(2), MPI_DOUBLE, neigh(14), neigh(14)+movetag, cube_comm, recv_req(13), ierr)
    call MPI_IRECV( buf_sendrecv14(1,2), ldims(2), MPI_DOUBLE, neigh(13), neigh(13)+movetag, cube_comm, recv_req(14), ierr)
    call MPI_IRECV( buf_sendrecv15(1,2), ldims(1), MPI_DOUBLE, neigh(16), neigh(16)+movetag, cube_comm, recv_req(15), ierr)
    call MPI_IRECV( buf_sendrecv16(1,2), ldims(1), MPI_DOUBLE, neigh(15), neigh(15)+movetag, cube_comm, recv_req(16), ierr)
    call MPI_IRECV( buf_sendrecv17(1,2), ldims(1), MPI_DOUBLE, neigh(18), neigh(18)+movetag, cube_comm, recv_req(17), ierr)
    call MPI_IRECV( buf_sendrecv18(1,2), ldims(1), MPI_DOUBLE, neigh(17), neigh(17)+movetag, cube_comm, recv_req(18), ierr)

    end subroutine mpipostrecv


    subroutine mpisendpops(aoptp)
    implicit none
    type(REALPTR), dimension(0:npop), intent(in)   :: aoptp
    integer         :: i, j, k
    Logical, save   :: isFirst = .true.


    if (isFirst) then
        isFirst = .false.
        call mpipostrecv
    endif

    ! Post sends - X axis
    buf_sendrecv1( 0*ldimyz+1:1*ldimyz, 1) = RESHAPE(aoptp(1)%p(maxx, miny:maxy,minz:maxz), [ ldimyz ] )
    buf_sendrecv1( 1*ldimyz+1:2*ldimyz, 1) = RESHAPE(aoptp(7)%p(maxx, miny:maxy,minz:maxz), [ ldimyz ] )
    buf_sendrecv1( 2*ldimyz+1:3*ldimyz, 1) = RESHAPE(aoptp(10)%p(maxx, miny:maxy,minz:maxz), [ ldimyz ] )
    buf_sendrecv1( 3*ldimyz+1:4*ldimyz, 1) = RESHAPE(aoptp(11)%p(maxx, miny:maxy,minz:maxz), [ ldimyz ] )
    buf_sendrecv1( 4*ldimyz+1:5*ldimyz, 1) = RESHAPE(aoptp(14)%p(maxx, miny:maxy,minz:maxz), [ ldimyz ] )

    buf_sendrecv2( 0*ldimyz+1:1*ldimyz, 1) = RESHAPE(aoptp(2)%p(minx, miny:maxy,minz:maxz), [ ldimyz ] )
    buf_sendrecv2( 1*ldimyz+1:2*ldimyz, 1) = RESHAPE(aoptp(8)%p(minx, miny:maxy,minz:maxz), [ ldimyz ] )
    buf_sendrecv2( 2*ldimyz+1:3*ldimyz, 1) = RESHAPE(aoptp(9)%p(minx, miny:maxy,minz:maxz), [ ldimyz ] )
    buf_sendrecv2( 3*ldimyz+1:4*ldimyz, 1) = RESHAPE(aoptp(12)%p(minx, miny:maxy,minz:maxz), [ ldimyz ] )
    buf_sendrecv2( 4*ldimyz+1:5*ldimyz, 1) = RESHAPE(aoptp(13)%p(minx, miny:maxy,minz:maxz), [ ldimyz ] )


    call MPI_ISEND( buf_sendrecv1(1,1), 5*ldimyz, MPI_DOUBLE, neigh(1), id_rank+movetag, cube_comm, send_req(1), ierr)
    call MPI_ISEND( buf_sendrecv2(1,1), 5*ldimyz, MPI_DOUBLE, neigh(2), id_rank+movetag, cube_comm, send_req(2), ierr)

    ! Post sends - Y axis
    buf_sendrecv3( 0*ldimxz+1:1*ldimxz, 1) = RESHAPE(aoptp(3)%p(minx:maxx, maxy, minz:maxz), [ ldimxz ] )
    buf_sendrecv3( 1*ldimxz+1:2*ldimxz, 1) = RESHAPE(aoptp(7)%p(minx:maxx, maxy, minz:maxz), [ ldimxz ] )
    buf_sendrecv3( 2*ldimxz+1:3*ldimxz, 1) = RESHAPE(aoptp(9)%p(minx:maxx, maxy, minz:maxz), [ ldimxz ] )
    buf_sendrecv3( 3*ldimxz+1:4*ldimxz, 1) = RESHAPE(aoptp(15)%p(minx:maxx, maxy, minz:maxz), [ ldimxz ] )
    buf_sendrecv3( 4*ldimxz+1:5*ldimxz, 1) = RESHAPE(aoptp(18)%p(minx:maxx, maxy, minz:maxz), [ ldimxz ] )

    buf_sendrecv4( 0*ldimxz+1:1*ldimxz, 1) = RESHAPE(aoptp(4)%p(minx:maxx, miny, minz:maxz), [ ldimxz ] )
    buf_sendrecv4( 1*ldimxz+1:2*ldimxz, 1) = RESHAPE(aoptp(8)%p(minx:maxx, miny, minz:maxz), [ ldimxz ] )
    buf_sendrecv4( 2*ldimxz+1:3*ldimxz, 1) = RESHAPE(aoptp(10)%p(minx:maxx, miny, minz:maxz), [ ldimxz ] )
    buf_sendrecv4( 3*ldimxz+1:4*ldimxz, 1) = RESHAPE(aoptp(16)%p(minx:maxx, miny, minz:maxz), [ ldimxz ] )
    buf_sendrecv4( 4*ldimxz+1:5*ldimxz, 1) = RESHAPE(aoptp(17)%p(minx:maxx, miny, minz:maxz), [ ldimxz ] )
    call MPI_ISEND( buf_sendrecv3(1,1), 5*ldimxz, MPI_DOUBLE, neigh(3), id_rank+movetag, cube_comm, send_req(3), ierr)
    call MPI_ISEND( buf_sendrecv4(1,1), 5*ldimxz, MPI_DOUBLE, neigh(4), id_rank+movetag, cube_comm, send_req(4), ierr)

    ! Post sends - Z axis
    buf_sendrecv5( 0*ldimxy+1:1*ldimxy, 1) = RESHAPE(aoptp(5)%p(minx:maxx, miny:maxy, maxz), [ ldimxy ] )
    buf_sendrecv5( 1*ldimxy+1:2*ldimxy, 1) = RESHAPE(aoptp(11)%p(minx:maxx, miny:maxy, maxz), [ ldimxy ] )
    buf_sendrecv5( 2*ldimxy+1:3*ldimxy, 1) = RESHAPE(aoptp(13)%p(minx:maxx, miny:maxy, maxz), [ ldimxy ] )
    buf_sendrecv5( 3*ldimxy+1:4*ldimxy, 1) = RESHAPE(aoptp(15)%p(minx:maxx, miny:maxy, maxz), [ ldimxy ] )
    buf_sendrecv5( 4*ldimxy+1:5*ldimxy, 1) = RESHAPE(aoptp(17)%p(minx:maxx, miny:maxy, maxz), [ ldimxy ] )

    buf_sendrecv6( 0*ldimxy+1:1*ldimxy, 1) = RESHAPE(aoptp(6)%p(minx:maxx, miny:maxy, minz), [ ldimxy ] )
    buf_sendrecv6( 1*ldimxy+1:2*ldimxy, 1) = RESHAPE(aoptp(12)%p(minx:maxx, miny:maxy, minz), [ ldimxy ] )
    buf_sendrecv6( 2*ldimxy+1:3*ldimxy, 1) = RESHAPE(aoptp(14)%p(minx:maxx, miny:maxy, minz), [ ldimxy ] )
    buf_sendrecv6( 3*ldimxy+1:4*ldimxy, 1) = RESHAPE(aoptp(16)%p(minx:maxx, miny:maxy, minz), [ ldimxy ] )
    buf_sendrecv6( 4*ldimxy+1:5*ldimxy, 1) = RESHAPE(aoptp(18)%p(minx:maxx, miny:maxy, minz), [ ldimxy ] )
    call MPI_ISEND( buf_sendrecv5(1,1), 5*ldimxy, MPI_DOUBLE, neigh(5), id_rank+movetag, cube_comm, send_req(5), ierr)
    call MPI_ISEND( buf_sendrecv6(1,1), 5*ldimxy, MPI_DOUBLE, neigh(6), id_rank+movetag, cube_comm, send_req(6), ierr)

    ! Post sends - z edges
    buf_sendrecv7(1:ldims(3), 1) = aoptp(7)%p(maxx,maxy, minz:maxz)
    buf_sendrecv8(1:ldims(3), 1) = aoptp(8)%p(minx,miny, minz:maxz)
    call MPI_ISEND( buf_sendrecv7(1,1), ldims(3), MPI_DOUBLE, neigh(7), id_rank+movetag, cube_comm, send_req(7), ierr)
    call MPI_ISEND( buf_sendrecv8(1,1), ldims(3), MPI_DOUBLE, neigh(8), id_rank+movetag, cube_comm, send_req(8), ierr)

    buf_sendrecv9 (1:ldims(3), 1) = aoptp(9)%p(minx,maxy, minz:maxz)
    buf_sendrecv10(1:ldims(3), 1) = aoptp(10)%p(maxx,miny, minz:maxz)

    call MPI_ISEND( buf_sendrecv9(1,1),  ldims(3), MPI_DOUBLE, neigh(9),  id_rank+movetag, cube_comm, send_req(9), ierr)
    call MPI_ISEND( buf_sendrecv10(1,1), ldims(3), MPI_DOUBLE, neigh(10), id_rank+movetag, cube_comm, send_req(10), ierr)

    ! Post sends - y edges
    buf_sendrecv11(1:ldims(2), 1) = aoptp(11)%p(maxx, miny:maxy, maxz)
    buf_sendrecv12(1:ldims(2), 1) = aoptp(12)%p(minx, miny:maxy, minz)
    call MPI_ISEND( buf_sendrecv11(1,1), ldims(2), MPI_DOUBLE, neigh(11), id_rank+movetag, cube_comm, send_req(11), ierr)
    call MPI_ISEND( buf_sendrecv12(1,1), ldims(2), MPI_DOUBLE, neigh(12), id_rank+movetag, cube_comm, send_req(12), ierr)

    buf_sendrecv13(1:ldims(2), 1) = aoptp(13)%p(minx, miny:maxy, maxz)
    buf_sendrecv14(1:ldims(2), 1) = aoptp(14)%p(maxx, miny:maxy, minz)
    call MPI_ISEND( buf_sendrecv13(1,1), ldims(2), MPI_DOUBLE, neigh(13), id_rank+movetag, cube_comm, send_req(13), ierr)
    call MPI_ISEND( buf_sendrecv14(1,1), ldims(2), MPI_DOUBLE, neigh(14), id_rank+movetag, cube_comm, send_req(14), ierr)

    ! Post sends - x edges
    buf_sendrecv15(1:ldims(1), 1) = aoptp(15)%p(minx:maxx, maxy, maxz)
    buf_sendrecv16(1:ldims(1), 1) = aoptp(16)%p(minx:maxx, miny, minz)
    call MPI_ISEND( buf_sendrecv15(1,1), ldims(1), MPI_DOUBLE, neigh(15), id_rank+movetag, cube_comm, send_req(15), ierr)
    call MPI_ISEND( buf_sendrecv16(1,1), ldims(1), MPI_DOUBLE, neigh(16), id_rank+movetag, cube_comm, send_req(16), ierr)

    buf_sendrecv17(1:ldims(1), 1) = aoptp(17)%p(minx:maxx, miny, maxz)
    buf_sendrecv18(1:ldims(1), 1) = aoptp(18)%p(minx:maxx, maxy, minz)
    call MPI_ISEND( buf_sendrecv17(1,1), ldims(1), MPI_DOUBLE, neigh(17), id_rank+movetag, cube_comm, send_req(17), ierr)
    call MPI_ISEND( buf_sendrecv18(1,1), ldims(1), MPI_DOUBLE, neigh(18), id_rank+movetag, cube_comm, send_req(18), ierr)

    end subroutine mpisendpops


    subroutine mpirecvpops(aoptp)
    implicit none
    type(REALPTR), dimension(0:npop), intent(in)   :: aoptp
    integer         :: statusArray(MPI_STATUS_SIZE,0:17)

    CALL MPI_WAITALL(18,recv_req,statusArray,ierr)
    if (ierr /= 0) write (6,*) "mpirecvpops: procID=", id_rank, "ierr=", ierr

    ! Wait & collect recv
    !! CALL MPI_WAITALL(1,recv_req(1:1),statusArray,ierr)
    if (neigh(2) /= MPI_PROC_NULL) then
    aoptp( 1)%p(minx-1, miny:maxy,minz:maxz) = RESHAPE(buf_sendrecv1(0*ldimyz+1:1*ldimyz, 2), [ ldims(2),ldims(3) ])
    aoptp( 7)%p(minx-1, miny:maxy,minz:maxz) = RESHAPE(buf_sendrecv1(1*ldimyz+1:2*ldimyz, 2), [ ldims(2),ldims(3) ])
    aoptp(10)%p(minx-1, miny:maxy,minz:maxz) = RESHAPE(buf_sendrecv1(2*ldimyz+1:3*ldimyz, 2), [ ldims(2),ldims(3) ])
    aoptp(11)%p(minx-1, miny:maxy,minz:maxz) = RESHAPE(buf_sendrecv1(3*ldimyz+1:4*ldimyz, 2), [ ldims(2),ldims(3) ])
    aoptp(14)%p(minx-1, miny:maxy,minz:maxz) = RESHAPE(buf_sendrecv1(4*ldimyz+1:5*ldimyz, 2), [ ldims(2),ldims(3) ])
    endif

    !! CALL MPI_WAITALL(1,recv_req(2:2),statusArray,ierr)
    if (neigh(1) /= MPI_PROC_NULL) then
    aoptp( 2)%p(maxx+1,miny:maxy,minz:maxz) = RESHAPE(buf_sendrecv2(0*ldimyz+1:1*ldimyz, 2), [ ldims(2),ldims(3) ])
    aoptp( 8)%p(maxx+1,miny:maxy,minz:maxz) = RESHAPE(buf_sendrecv2(1*ldimyz+1:2*ldimyz, 2), [ ldims(2),ldims(3) ])
    aoptp( 9)%p(maxx+1,miny:maxy,minz:maxz) = RESHAPE(buf_sendrecv2(2*ldimyz+1:3*ldimyz, 2), [ ldims(2),ldims(3) ])
    aoptp(12)%p(maxx+1,miny:maxy,minz:maxz) = RESHAPE(buf_sendrecv2(3*ldimyz+1:4*ldimyz, 2), [ ldims(2),ldims(3) ])
    aoptp(13)%p(maxx+1,miny:maxy,minz:maxz) = RESHAPE(buf_sendrecv2(4*ldimyz+1:5*ldimyz, 2), [ ldims(2),ldims(3) ])
    endif

    !! CALL MPI_WAITALL(1,recv_req(3:3),statusArray,ierr)
    if (neigh(4) /= MPI_PROC_NULL) then
    aoptp( 3)%p(minx:maxx, miny-1, minz:maxz) = RESHAPE(buf_sendrecv3(0*ldimxz+1:1*ldimxz, 2), [ ldims(1),ldims(3) ])
    aoptp( 7)%p(minx:maxx, miny-1, minz:maxz) = RESHAPE(buf_sendrecv3(1*ldimxz+1:2*ldimxz, 2), [ ldims(1),ldims(3) ])
    aoptp( 9)%p(minx:maxx, miny-1, minz:maxz) = RESHAPE(buf_sendrecv3(2*ldimxz+1:3*ldimxz, 2), [ ldims(1),ldims(3) ])
    aoptp(15)%p(minx:maxx, miny-1, minz:maxz) = RESHAPE(buf_sendrecv3(3*ldimxz+1:4*ldimxz, 2), [ ldims(1),ldims(3) ])
    aoptp(18)%p(minx:maxx, miny-1, minz:maxz) = RESHAPE(buf_sendrecv3(4*ldimxz+1:5*ldimxz, 2), [ ldims(1),ldims(3) ])
    endif

    !! CALL MPI_WAITALL(1,recv_req(4:4),statusArray,ierr)
    if (neigh(3) /= MPI_PROC_NULL) then
    aoptp( 4)%p(minx:maxx, maxy+1, minz:maxz) = RESHAPE(buf_sendrecv4(0*ldimxz+1:1*ldimxz, 2), [ ldims(1),ldims(3) ])
    aoptp( 8)%p(minx:maxx, maxy+1, minz:maxz) = RESHAPE(buf_sendrecv4(1*ldimxz+1:2*ldimxz, 2), [ ldims(1),ldims(3) ])
    aoptp(10)%p(minx:maxx, maxy+1, minz:maxz) = RESHAPE(buf_sendrecv4(2*ldimxz+1:3*ldimxz, 2), [ ldims(1),ldims(3) ])
    aoptp(16)%p(minx:maxx, maxy+1, minz:maxz) = RESHAPE(buf_sendrecv4(3*ldimxz+1:4*ldimxz, 2), [ ldims(1),ldims(3) ])
    aoptp(17)%p(minx:maxx, maxy+1, minz:maxz) = RESHAPE(buf_sendrecv4(4*ldimxz+1:5*ldimxz, 2), [ ldims(1),ldims(3) ])
    endif

    !! CALL MPI_WAITALL(1,recv_req(5:5),statusArray,ierr)
    if (neigh(6) /= MPI_PROC_NULL) then
    aoptp( 5)%p(minx:maxx, miny:maxy, minz-1) = RESHAPE(buf_sendrecv5(0*ldimxy+1:1*ldimxy, 2), [ ldims(1),ldims(2) ])
    aoptp(11)%p(minx:maxx, miny:maxy, minz-1) = RESHAPE(buf_sendrecv5(1*ldimxy+1:2*ldimxy, 2), [ ldims(1),ldims(2) ])
    aoptp(13)%p(minx:maxx, miny:maxy, minz-1) = RESHAPE(buf_sendrecv5(2*ldimxy+1:3*ldimxy, 2), [ ldims(1),ldims(2) ])
    aoptp(15)%p(minx:maxx, miny:maxy, minz-1) = RESHAPE(buf_sendrecv5(3*ldimxy+1:4*ldimxy, 2), [ ldims(1),ldims(2) ])
    aoptp(17)%p(minx:maxx, miny:maxy, minz-1) = RESHAPE(buf_sendrecv5(4*ldimxy+1:5*ldimxy, 2), [ ldims(1),ldims(2) ])
    endif

    !! CALL MPI_WAITALL(1,recv_req(6:6),statusArray,ierr)
    if (neigh(5) /= MPI_PROC_NULL) then
    aoptp( 6)%p(minx:maxx, miny:maxy, maxz+1) = RESHAPE(buf_sendrecv6(0*ldimxy+1:1*ldimxy, 2), [ ldims(1),ldims(2) ])
    aoptp(12)%p(minx:maxx, miny:maxy, maxz+1) = RESHAPE(buf_sendrecv6(1*ldimxy+1:2*ldimxy, 2), [ ldims(1),ldims(2) ])
    aoptp(14)%p(minx:maxx, miny:maxy, maxz+1) = RESHAPE(buf_sendrecv6(2*ldimxy+1:3*ldimxy, 2), [ ldims(1),ldims(2) ])
    aoptp(16)%p(minx:maxx, miny:maxy, maxz+1) = RESHAPE(buf_sendrecv6(3*ldimxy+1:4*ldimxy, 2), [ ldims(1),ldims(2) ])
    aoptp(18)%p(minx:maxx, miny:maxy, maxz+1) = RESHAPE(buf_sendrecv6(4*ldimxy+1:5*ldimxy, 2), [ ldims(1),ldims(2) ])
    endif

    ! +1,+1 & -1,-1
    !! CALL MPI_WAITALL(1,recv_req(7:7),statusArray,ierr)
    if (neigh(8) /= MPI_PROC_NULL) then
    aoptp(7)%p(minx-1,miny-1, minz:maxz) = buf_sendrecv7(:,2)
    endif

    !! CALL MPI_WAITALL(1,recv_req(8:8),statusArray,ierr)
    if (neigh(7) /= MPI_PROC_NULL) then
    aoptp(8)%p(maxx+1,maxy+1, minz:maxz) = buf_sendrecv8(:,2)
    endif

    ! -1,+1 & +1,-1
    !! CALL MPI_WAITALL(1,recv_req(9:9),statusArray,ierr)
    if (neigh(10) /= MPI_PROC_NULL) then
    aoptp(9)%p(maxx+1,miny-1, minz:maxz) = buf_sendrecv9(:,2)
    endif

    !! CALL MPI_WAITALL(1,recv_req(10:10),statusArray,ierr)
    if (neigh(9) /= MPI_PROC_NULL) then
    aoptp(10)%p(minx-1,maxy+1, minz:maxz) = buf_sendrecv10(:,2)
    endif


    ! 1,0,1 & -1,0,-1
    !! CALL MPI_WAITALL(1,recv_req(11:11),statusArray,ierr)
    if (neigh(12) /= MPI_PROC_NULL) then
    aoptp(11)%p(minx-1,miny:maxy, minz-1) = buf_sendrecv11(:,2)
    endif

    !! CALL MPI_WAITALL(1,recv_req(12:12),statusArray,ierr)
    if (neigh(11) /= MPI_PROC_NULL) then
    aoptp(12)%p(maxx+1,miny:maxy, maxz+1) = buf_sendrecv12(:,2)
    endif

    ! -1,0,1 & 1,0,-1
    !! CALL MPI_WAITALL(1,recv_req(13:13),statusArray,ierr)
    if (neigh(14) /= MPI_PROC_NULL) then
    aoptp(13)%p(maxx+1,miny:maxy, minz-1) = buf_sendrecv13(:,2)
    endif

    !! CALL MPI_WAITALL(1,recv_req(14:14),statusArray,ierr)
    if (neigh(13) /= MPI_PROC_NULL) then
    aoptp(14)%p(minx-1,miny:maxy, maxz+1) = buf_sendrecv14(:,2)
    endif

    ! 0,1,1 & 0,-1,-1
    !! CALL MPI_WAITALL(1,recv_req(15:15),statusArray,ierr)
    if (neigh(16) /= MPI_PROC_NULL) then
    aoptp(15)%p(minx:maxx, miny-1, minz-1) = buf_sendrecv15(:,2)
    endif

    !! CALL MPI_WAITALL(1,recv_req(16:16),statusArray,ierr)
    if (neigh(15) /= MPI_PROC_NULL) then
    aoptp(16)%p(minx:maxx, maxy+1, maxz+1) = buf_sendrecv16(:,2)
    endif

    ! 0,-1,1 & 0,1,-1
    !! CALL MPI_WAITALL(1,recv_req(17:17),statusArray,ierr)
    if (neigh(18) /= MPI_PROC_NULL) then
    aoptp(17)%p(minx:maxx, maxy+1, minz-1) = buf_sendrecv17(:,2)
    endif

    !! CALL MPI_WAITALL(1,recv_req(18:18),statusArray,ierr)
    if (neigh(17) /= MPI_PROC_NULL) then
    aoptp(18)%p(minx:maxx, miny-1, maxz+1) = buf_sendrecv18(:,2)
    endif

    CALL MPI_WAITALL(18,send_req,statusArray,ierr)
    call mpipostrecv

    end subroutine mpirecvpops


    subroutine mpibounceback_x(aoptp)
    implicit none
    type(REALPTR), dimension(0:npop), intent(inout)   :: aoptp

    if(maxx == nx)then
        aoptp( 2)%p(nx+1, miny-1:maxy+1, minz-1:maxz+1) = aoptp( 1)%p(nx+1, miny-1:maxy+1, minz-1:maxz+1)
        aoptp( 8)%p(nx+1, miny-1:maxy+1, minz-1:maxz+1) = aoptp( 7)%p(nx+1, miny-1:maxy+1, minz-1:maxz+1)
        aoptp( 9)%p(nx+1, miny-1:maxy+1, minz-1:maxz+1) = aoptp(10)%p(nx+1, miny-1:maxy+1, minz-1:maxz+1)
        aoptp(12)%p(nx+1, miny-1:maxy+1, minz-1:maxz+1) = aoptp(11)%p(nx+1, miny-1:maxy+1, minz-1:maxz+1)
        aoptp(13)%p(nx+1, miny-1:maxy+1, minz-1:maxz+1) = aoptp(14)%p(nx+1, miny-1:maxy+1, minz-1:maxz+1)
    end if

    if(minx == 1)then
        aoptp( 1)%p(0, miny-1:maxy+1, minz-1:maxz+1) = aoptp( 2)%p(0, miny-1:maxy+1, minz-1:maxz+1)
        aoptp( 7)%p(0, miny-1:maxy+1, minz-1:maxz+1) = aoptp( 8)%p(0, miny-1:maxy+1, minz-1:maxz+1)
        aoptp(10)%p(0, miny-1:maxy+1, minz-1:maxz+1) = aoptp( 9)%p(0, miny-1:maxy+1, minz-1:maxz+1)
        aoptp(11)%p(0, miny-1:maxy+1, minz-1:maxz+1) = aoptp(12)%p(0, miny-1:maxy+1, minz-1:maxz+1)
        aoptp(14)%p(0, miny-1:maxy+1, minz-1:maxz+1) = aoptp(13)%p(0, miny-1:maxy+1, minz-1:maxz+1)
    end if
    end subroutine mpibounceback_x



    subroutine mpibounceback_y(aoptp)
    implicit none
    type(REALPTR), dimension(0:npop), intent(inout)   :: aoptp

    if(maxy == ny)then
        aoptp( 4)%p(minx-1:maxx+1, ny+1, minz-1:maxz+1) = aoptp( 3)%p(minx-1:maxx+1, ny+1, minz-1:maxz+1)
        aoptp( 8)%p(minx-1:maxx+1, ny+1, minz-1:maxz+1) = aoptp( 7)%p(minx-1:maxx+1, ny+1, minz-1:maxz+1)
        aoptp(10)%p(minx-1:maxx+1, ny+1, minz-1:maxz+1) = aoptp( 9)%p(minx-1:maxx+1, ny+1, minz-1:maxz+1)
        aoptp(16)%p(minx-1:maxx+1, ny+1, minz-1:maxz+1) = aoptp(15)%p(minx-1:maxx+1, ny+1, minz-1:maxz+1)
        aoptp(17)%p(minx-1:maxx+1, ny+1, minz-1:maxz+1) = aoptp(18)%p(minx-1:maxx+1, ny+1, minz-1:maxz+1)
    end if

    if(miny == 1)then
        aoptp( 3)%p(minx-1:maxx+1, 0, minz-1:maxz+1) = aoptp( 4)%p(minx-1:maxx+1, 0, minz-1:maxz+1)
        aoptp( 7)%p(minx-1:maxx+1, 0, minz-1:maxz+1) = aoptp( 8)%p(minx-1:maxx+1, 0, minz-1:maxz+1)
        aoptp( 9)%p(minx-1:maxx+1, 0, minz-1:maxz+1) = aoptp(10)%p(minx-1:maxx+1, 0, minz-1:maxz+1)
        aoptp(15)%p(minx-1:maxx+1, 0, minz-1:maxz+1) = aoptp(16)%p(minx-1:maxx+1, 0, minz-1:maxz+1)
        aoptp(18)%p(minx-1:maxx+1, 0, minz-1:maxz+1) = aoptp(17)%p(minx-1:maxx+1, 0, minz-1:maxz+1)
    end if
    end subroutine mpibounceback_y


    subroutine mpibounceback_z(aoptp)
    implicit none
    type(REALPTR), dimension(0:npop), intent(inout)   :: aoptp

    if(maxz == nz)then
        aoptp(6) %p(minx-1:maxx+1, miny-1:maxy+1, nz+1) = aoptp(5) %p(minx-1:maxx+1, miny-1:maxy+1, nz+1)
        aoptp(12)%p(minx-1:maxx+1, miny-1:maxy+1, nz+1) = aoptp(11)%p(minx-1:maxx+1, miny-1:maxy+1, nz+1)
        aoptp(14)%p(minx-1:maxx+1, miny-1:maxy+1, nz+1) = aoptp(13)%p(minx-1:maxx+1, miny-1:maxy+1, nz+1)
        aoptp(16)%p(minx-1:maxx+1, miny-1:maxy+1, nz+1) = aoptp(15)%p(minx-1:maxx+1, miny-1:maxy+1, nz+1)
        aoptp(18)%p(minx-1:maxx+1, miny-1:maxy+1, nz+1) = aoptp(17)%p(minx-1:maxx+1, miny-1:maxy+1, nz+1)
    end if

    if(minz == 1)then
        aoptp(5) %p(minx-1:maxx+1, miny-1:maxy+1, 0) = aoptp(6) %p(minx-1:maxx+1, miny-1:maxy+1, 0)
        aoptp(11)%p(minx-1:maxx+1, miny-1:maxy+1, 0) = aoptp(12)%p(minx-1:maxx+1, miny-1:maxy+1, 0)
        aoptp(13)%p(minx-1:maxx+1, miny-1:maxy+1, 0) = aoptp(14)%p(minx-1:maxx+1, miny-1:maxy+1, 0)
        aoptp(15)%p(minx-1:maxx+1, miny-1:maxy+1, 0) = aoptp(16)%p(minx-1:maxx+1, miny-1:maxy+1, 0)
        aoptp(17)%p(minx-1:maxx+1, miny-1:maxy+1, 0) = aoptp(18)%p(minx-1:maxx+1, miny-1:maxy+1, 0)
    end if
    end subroutine mpibounceback_z


    subroutine mpibounceback(aoptp)
    implicit none
    type(REALPTR), dimension(0:npop), intent(inout)   :: aoptp

    if (.not. period(1)) call mpibounceback_x(aoptp)
    if (.not. period(2)) call mpibounceback_y(aoptp)
    if (.not. period(3)) call mpibounceback_z(aoptp)

    end subroutine mpibounceback


    subroutine mpiwaitallrecv
    implicit none
    integer         :: statusArray(MPI_STATUS_SIZE,0:17)
    CALL MPI_WAITALL(18,recv_req,statusArray,ierr)
    end subroutine mpiwaitallrecv


    subroutine mpipostrecv_hvar
    implicit none

    ! Post MPI_IRECV
    call MPI_IRECV( buf_sendrecv1_hvar(1,2),  ldimyz, MPI_DOUBLE, neigh(2),  & 
            neigh(2)+hvartag, cube_comm, recv_req_hvar(1),  ierr)
    call MPI_IRECV( buf_sendrecv2_hvar(1,2),  ldimyz, MPI_DOUBLE, neigh(1),  & 
            neigh(1)+hvartag, cube_comm, recv_req_hvar(2),  ierr)
    call MPI_IRECV( buf_sendrecv3_hvar(1,2),  ldimxz, MPI_DOUBLE, neigh(4),  &
            neigh(4)+hvartag, cube_comm, recv_req_hvar(3),  ierr)
    call MPI_IRECV( buf_sendrecv4_hvar(1,2),  ldimxz, MPI_DOUBLE, neigh(3),  &
            neigh(3)+hvartag, cube_comm, recv_req_hvar(4),  ierr)
    call MPI_IRECV( buf_sendrecv5_hvar(1,2),  ldimxy, MPI_DOUBLE, neigh(6),  &
            neigh(6)+hvartag, cube_comm, recv_req_hvar(5),  ierr)
    call MPI_IRECV( buf_sendrecv6_hvar(1,2),  ldimxy, MPI_DOUBLE, neigh(5),  &
            neigh(5)+hvartag, cube_comm, recv_req_hvar(6),  ierr)

    call MPI_IRECV( buf_sendrecv7_hvar(1,2),  ldims(3), MPI_DOUBLE, neigh(8),  &
            neigh(8) +hvartag, cube_comm, recv_req_hvar(7),  ierr)
    call MPI_IRECV( buf_sendrecv8_hvar(1,2),  ldims(3), MPI_DOUBLE, neigh(7),  &
            neigh(7) +hvartag, cube_comm, recv_req_hvar(8),  ierr)
    call MPI_IRECV( buf_sendrecv9_hvar(1,2),  ldims(3), MPI_DOUBLE, neigh(10), &
            neigh(10)+hvartag, cube_comm, recv_req_hvar(9), ierr)
    call MPI_IRECV( buf_sendrecv10_hvar(1,2), ldims(3), MPI_DOUBLE, neigh(9),  &
            neigh(9) +hvartag, cube_comm, recv_req_hvar(10), ierr)
    call MPI_IRECV( buf_sendrecv11_hvar(1,2), ldims(2), MPI_DOUBLE, neigh(12), &
            neigh(12)+hvartag, cube_comm, recv_req_hvar(11), ierr)
    call MPI_IRECV( buf_sendrecv12_hvar(1,2), ldims(2), MPI_DOUBLE, neigh(11), &
            neigh(11)+hvartag, cube_comm, recv_req_hvar(12), ierr)
    call MPI_IRECV( buf_sendrecv13_hvar(1,2), ldims(2), MPI_DOUBLE, neigh(14), &
            neigh(14)+hvartag, cube_comm, recv_req_hvar(13), ierr)
    call MPI_IRECV( buf_sendrecv14_hvar(1,2), ldims(2), MPI_DOUBLE, neigh(13), &
            neigh(13)+hvartag, cube_comm, recv_req_hvar(14), ierr)
    call MPI_IRECV( buf_sendrecv15_hvar(1,2), ldims(1), MPI_DOUBLE, neigh(16), &
            neigh(16)+hvartag, cube_comm, recv_req_hvar(15), ierr)
    call MPI_IRECV( buf_sendrecv16_hvar(1,2), ldims(1), MPI_DOUBLE, neigh(15), &
            neigh(15)+hvartag, cube_comm, recv_req_hvar(16), ierr)
    call MPI_IRECV( buf_sendrecv17_hvar(1,2), ldims(1), MPI_DOUBLE, neigh(18), &
            neigh(18)+hvartag, cube_comm, recv_req_hvar(17), ierr)
    call MPI_IRECV( buf_sendrecv18_hvar(1,2), ldims(1), MPI_DOUBLE, neigh(17), &
            neigh(17)+hvartag, cube_comm, recv_req_hvar(18), ierr)

    end subroutine mpipostrecv_hvar


    subroutine mpisend_hvar(dtemp, isFirst)
    implicit none
    REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp
    Logical, intent(in)   :: isFirst


    if (isFirst) then
        call mpipostrecv_hvar
    endif

    ! Post sends - X axis
    buf_sendrecv1_hvar( 1:ldimyz, 1) = RESHAPE(dtemp(maxx, miny:maxy,minz:maxz), [ ldimyz ] )
    buf_sendrecv2_hvar( 1:ldimyz, 1) = RESHAPE(dtemp(minx, miny:maxy,minz:maxz), [ ldimyz ] )
    call MPI_ISEND( buf_sendrecv1_hvar(1,1), ldimyz, MPI_DOUBLE, neigh(1), id_rank+hvartag, cube_comm, send_req_hvar(1), ierr)
    call MPI_ISEND( buf_sendrecv2_hvar(1,1), ldimyz, MPI_DOUBLE, neigh(2), id_rank+hvartag, cube_comm, send_req_hvar(2), ierr)

    ! Post sends - Y axis
    buf_sendrecv3_hvar( 1:ldimxz, 1) = RESHAPE(dtemp(minx:maxx, maxy, minz:maxz), [ ldimxz ] )
    buf_sendrecv4_hvar( 1:ldimxz, 1) = RESHAPE(dtemp(minx:maxx, miny, minz:maxz), [ ldimxz ] )
    call MPI_ISEND( buf_sendrecv3_hvar(1,1), ldimxz, MPI_DOUBLE, neigh(3), id_rank+hvartag, cube_comm, send_req_hvar(3), ierr)
    call MPI_ISEND( buf_sendrecv4_hvar(1,1), ldimxz, MPI_DOUBLE, neigh(4), id_rank+hvartag, cube_comm, send_req_hvar(4), ierr)

    ! Post sends - Z axis
    buf_sendrecv5_hvar( 1:ldimxy, 1) = RESHAPE(dtemp(minx:maxx, miny:maxy, maxz), [ ldimxy ] )
    buf_sendrecv6_hvar( 1:ldimxy, 1) = RESHAPE(dtemp(minx:maxx, miny:maxy, minz), [ ldimxy ] )
    call MPI_ISEND( buf_sendrecv5_hvar(1,1), ldimxy, MPI_DOUBLE, neigh(5), id_rank+hvartag, cube_comm, send_req_hvar(5), ierr)
    call MPI_ISEND( buf_sendrecv6_hvar(1,1), ldimxy, MPI_DOUBLE, neigh(6), id_rank+hvartag, cube_comm, send_req_hvar(6), ierr)

    ! Post sends - z edges
    buf_sendrecv7_hvar(1:ldims(3), 1) = dtemp(maxx,maxy, minz:maxz)
    buf_sendrecv8_hvar(1:ldims(3), 1) = dtemp(minx,miny, minz:maxz)
    call MPI_ISEND( buf_sendrecv7_hvar(1,1), ldims(3), MPI_DOUBLE, neigh(7), id_rank+hvartag, cube_comm, send_req_hvar(7), ierr)
    call MPI_ISEND( buf_sendrecv8_hvar(1,1), ldims(3), MPI_DOUBLE, neigh(8), id_rank+hvartag, cube_comm, send_req_hvar(8), ierr)

    buf_sendrecv9_hvar (1:ldims(3), 1) = dtemp(minx,maxy, minz:maxz)
    buf_sendrecv10_hvar(1:ldims(3), 1) = dtemp(maxx,miny, minz:maxz)
    call MPI_ISEND( buf_sendrecv9_hvar(1,1),  ldims(3), MPI_DOUBLE, neigh(9),  id_rank+hvartag, cube_comm, send_req_hvar(9), ierr)
    call MPI_ISEND( buf_sendrecv10_hvar(1,1), ldims(3), MPI_DOUBLE, neigh(10), id_rank+hvartag, cube_comm, send_req_hvar(10), ierr)

    ! Post sends - y edges
    buf_sendrecv11_hvar(1:ldims(2), 1) = dtemp(maxx, miny:maxy, maxz)
    buf_sendrecv12_hvar(1:ldims(2), 1) = dtemp(minx, miny:maxy, minz)
    call MPI_ISEND( buf_sendrecv11_hvar(1,1), ldims(2), MPI_DOUBLE, neigh(11), id_rank+hvartag, cube_comm, send_req_hvar(11), ierr)
    call MPI_ISEND( buf_sendrecv12_hvar(1,1), ldims(2), MPI_DOUBLE, neigh(12), id_rank+hvartag, cube_comm, send_req_hvar(12), ierr)

    buf_sendrecv13_hvar(1:ldims(2), 1) = dtemp(minx, miny:maxy, maxz)
    buf_sendrecv14_hvar(1:ldims(2), 1) = dtemp(maxx, miny:maxy, minz)
    call MPI_ISEND( buf_sendrecv13_hvar(1,1), ldims(2), MPI_DOUBLE, neigh(13), id_rank+hvartag, cube_comm, send_req_hvar(13), ierr)
    call MPI_ISEND( buf_sendrecv14_hvar(1,1), ldims(2), MPI_DOUBLE, neigh(14), id_rank+hvartag, cube_comm, send_req_hvar(14), ierr)

    ! Post sends - x edges
    buf_sendrecv15_hvar(1:ldims(1), 1) = dtemp(minx:maxx, maxy, maxz)
    buf_sendrecv16_hvar(1:ldims(1), 1) = dtemp(minx:maxx, miny, minz)
    call MPI_ISEND( buf_sendrecv15_hvar(1,1), ldims(1), MPI_DOUBLE, neigh(15), id_rank+hvartag, cube_comm, send_req_hvar(15), ierr)
    call MPI_ISEND( buf_sendrecv16_hvar(1,1), ldims(1), MPI_DOUBLE, neigh(16), id_rank+hvartag, cube_comm, send_req_hvar(16), ierr)

    buf_sendrecv17_hvar(1:ldims(1), 1) = dtemp(minx:maxx, miny, maxz)
    buf_sendrecv18_hvar(1:ldims(1), 1) = dtemp(minx:maxx, maxy, minz)
    call MPI_ISEND( buf_sendrecv17_hvar(1,1), ldims(1), MPI_DOUBLE, neigh(17), id_rank+hvartag, cube_comm, send_req_hvar(17), ierr)
    call MPI_ISEND( buf_sendrecv18_hvar(1,1), ldims(1), MPI_DOUBLE, neigh(18), id_rank+hvartag, cube_comm, send_req_hvar(18), ierr)

    end subroutine mpisend_hvar


    subroutine mpirecv_hvar(dtemp)
    implicit none
    REAL(KIND=PRC), dimension(:,:,:), allocatable :: dtemp
    integer         :: statusArray(MPI_STATUS_SIZE,0:17)

    CALL MPI_WAITALL(18,recv_req_hvar,statusArray,ierr)
    if (ierr /= 0) write (6,*) "mpirecv_hvar: procID=", id_rank, "ierr=", ierr

    ! Wait & collect recv
    if (neigh(2) /= MPI_PROC_NULL) then
    dtemp(minx-1, miny:maxy, minz:maxz) = RESHAPE(buf_sendrecv1_hvar(1:ldimyz, 2), [ ldims(2),ldims(3) ])
    endif

    if (neigh(1) /= MPI_PROC_NULL) then
    dtemp(maxx+1, miny:maxy, minz:maxz) = RESHAPE(buf_sendrecv2_hvar(1:ldimyz, 2), [ ldims(2),ldims(3) ])
    endif

    if (neigh(4) /= MPI_PROC_NULL) then
    dtemp(minx:maxx, miny-1, minz:maxz) = RESHAPE(buf_sendrecv3_hvar(1:ldimxz, 2), [ ldims(1),ldims(3) ])
    endif

    if (neigh(3) /= MPI_PROC_NULL) then
    dtemp(minx:maxx, maxy+1, minz:maxz) = RESHAPE(buf_sendrecv4_hvar(1:ldimxz, 2), [ ldims(1),ldims(3) ])
    endif

    if (neigh(6) /= MPI_PROC_NULL) then
    dtemp(minx:maxx, miny:maxy, minz-1) = RESHAPE(buf_sendrecv5_hvar(1:ldimxy, 2), [ ldims(1),ldims(2) ])
    endif

    if (neigh(5) /= MPI_PROC_NULL) then
    dtemp(minx:maxx, miny:maxy, maxz+1) = RESHAPE(buf_sendrecv6_hvar(1:ldimxy, 2), [ ldims(1),ldims(2) ])
    endif

    ! +1,+1 & -1,-1
    if (neigh(8) /= MPI_PROC_NULL) then
    dtemp(minx-1,miny-1, minz:maxz) = buf_sendrecv7_hvar(:,2)
    endif

    if (neigh(7) /= MPI_PROC_NULL) then
    dtemp(maxx+1,maxy+1, minz:maxz) = buf_sendrecv8_hvar(:,2)
    endif

    ! -1,+1 & +1,-1
    if (neigh(10) /= MPI_PROC_NULL) then
    dtemp(maxx+1,miny-1, minz:maxz) = buf_sendrecv9_hvar(:,2)
    endif

    if (neigh(9) /= MPI_PROC_NULL) then
    dtemp(minx-1,maxy+1, minz:maxz) = buf_sendrecv10_hvar(:,2)
    endif


    ! 1,0,1 & -1,0,-1
    if (neigh(12) /= MPI_PROC_NULL) then
    dtemp(minx-1,miny:maxy, minz-1) = buf_sendrecv11_hvar(:,2)
    endif

    if (neigh(11) /= MPI_PROC_NULL) then
    dtemp(maxx+1,miny:maxy, maxz+1) = buf_sendrecv12_hvar(:,2)
    endif

    ! -1,0,1 & 1,0,-1
    if (neigh(14) /= MPI_PROC_NULL) then
    dtemp(maxx+1,miny:maxy, minz-1) = buf_sendrecv13_hvar(:,2)
    endif

    if (neigh(13) /= MPI_PROC_NULL) then
    dtemp(minx-1,miny:maxy, maxz+1) = buf_sendrecv14_hvar(:,2)
    endif

    ! 0,1,1 & 0,-1,-1
    if (neigh(16) /= MPI_PROC_NULL) then
    dtemp(minx:maxx, miny-1, minz-1) = buf_sendrecv15_hvar(:,2)
    endif

    if (neigh(15) /= MPI_PROC_NULL) then
    dtemp(minx:maxx, maxy+1, maxz+1) = buf_sendrecv16_hvar(:,2)
    endif

    ! 0,-1,1 & 0,1,-1
    if (neigh(18) /= MPI_PROC_NULL) then
    dtemp(minx:maxx, maxy+1, minz-1) = buf_sendrecv17_hvar(:,2)
    endif

    if (neigh(17) /= MPI_PROC_NULL) then
    dtemp(minx:maxx, miny-1, maxz+1) = buf_sendrecv18_hvar(:,2)
    endif

    CALL MPI_WAITALL(18,send_req_hvar,statusArray,ierr)
    call mpipostrecv_hvar

    end subroutine mpirecv_hvar


    subroutine mpipostrecv_isfluid
    implicit none

    ! Post MPI_IRECV
    call MPI_IRECV( buf_sendrecv1_isfluid(1,2),  ldimyz, MPI_INTEGER, neigh(2),  &
            neigh(2)+fldtag, cube_comm, recv_req_isfluid(1),  ierr)
    call MPI_IRECV( buf_sendrecv2_isfluid(1,2),  ldimyz, MPI_INTEGER, neigh(1),  &
            neigh(1)+fldtag, cube_comm, recv_req_isfluid(2),  ierr)
    call MPI_IRECV( buf_sendrecv3_isfluid(1,2),  ldimxz, MPI_INTEGER, neigh(4),  &
            neigh(4)+fldtag, cube_comm, recv_req_isfluid(3),  ierr)
    call MPI_IRECV( buf_sendrecv4_isfluid(1,2),  ldimxz, MPI_INTEGER, neigh(3),  &
            neigh(3)+fldtag, cube_comm, recv_req_isfluid(4),  ierr)
    call MPI_IRECV( buf_sendrecv5_isfluid(1,2),  ldimxy, MPI_INTEGER, neigh(6),  &
            neigh(6)+fldtag, cube_comm, recv_req_isfluid(5),  ierr)
    call MPI_IRECV( buf_sendrecv6_isfluid(1,2),  ldimxy, MPI_INTEGER, neigh(5),  &
            neigh(5)+fldtag, cube_comm, recv_req_isfluid(6),  ierr)

    call MPI_IRECV( buf_sendrecv7_isfluid(1,2),  ldims(3), MPI_INTEGER, neigh(8),  &
            neigh(8) +fldtag, cube_comm, recv_req_isfluid(7),  ierr)
    call MPI_IRECV( buf_sendrecv8_isfluid(1,2),  ldims(3), MPI_INTEGER, neigh(7),  &
            neigh(7) +fldtag, cube_comm, recv_req_isfluid(8),  ierr)
    call MPI_IRECV( buf_sendrecv9_isfluid(1,2),  ldims(3), MPI_INTEGER, neigh(10), &
            neigh(10)+fldtag, cube_comm, recv_req_isfluid(9), ierr)
    call MPI_IRECV( buf_sendrecv10_isfluid(1,2), ldims(3), MPI_INTEGER, neigh(9),  &
            neigh(9) +fldtag, cube_comm, recv_req_isfluid(10), ierr)
    call MPI_IRECV( buf_sendrecv11_isfluid(1,2), ldims(2), MPI_INTEGER, neigh(12), &
            neigh(12)+fldtag, cube_comm, recv_req_isfluid(11), ierr)
    call MPI_IRECV( buf_sendrecv12_isfluid(1,2), ldims(2), MPI_INTEGER, neigh(11), &
            neigh(11)+fldtag, cube_comm, recv_req_isfluid(12), ierr)
    call MPI_IRECV( buf_sendrecv13_isfluid(1,2), ldims(2), MPI_INTEGER, neigh(14), &
            neigh(14)+fldtag, cube_comm, recv_req_isfluid(13), ierr)
    call MPI_IRECV( buf_sendrecv14_isfluid(1,2), ldims(2), MPI_INTEGER, neigh(13), &
            neigh(13)+fldtag, cube_comm, recv_req_isfluid(14), ierr)
    call MPI_IRECV( buf_sendrecv15_isfluid(1,2), ldims(1), MPI_INTEGER, neigh(16), &
            neigh(16)+fldtag, cube_comm, recv_req_isfluid(15), ierr)
    call MPI_IRECV( buf_sendrecv16_isfluid(1,2), ldims(1), MPI_INTEGER, neigh(15), &
            neigh(15)+fldtag, cube_comm, recv_req_isfluid(16), ierr)
    call MPI_IRECV( buf_sendrecv17_isfluid(1,2), ldims(1), MPI_INTEGER, neigh(18), &
            neigh(18)+fldtag, cube_comm, recv_req_isfluid(17), ierr)
    call MPI_IRECV( buf_sendrecv18_isfluid(1,2), ldims(1), MPI_INTEGER, neigh(17), &
            neigh(17)+fldtag, cube_comm, recv_req_isfluid(18), ierr)

    end subroutine mpipostrecv_isfluid


    subroutine mpisend_isfluid(dtemp, isFirst)
    implicit none
    INTEGER(kind=1), dimension(:,:,:), allocatable :: dtemp
    Logical, intent(in)   :: isFirst


    if (isFirst) then
        call mpipostrecv_isfluid
    endif

    ! Post sends - X axis
    buf_sendrecv1_isfluid( 1:ldimyz, 1) = RESHAPE(dtemp(maxx, miny:maxy,minz:maxz), [ ldimyz ] )
    buf_sendrecv2_isfluid( 1:ldimyz, 1) = RESHAPE(dtemp(minx, miny:maxy,minz:maxz), [ ldimyz ] )
    call MPI_ISEND( buf_sendrecv1_isfluid(1,1), ldimyz, MPI_INTEGER, neigh(1), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(1), ierr)
    call MPI_ISEND( buf_sendrecv1_isfluid(1,1), ldimyz, MPI_INTEGER, neigh(1), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(1), ierr)
    call MPI_ISEND( buf_sendrecv2_isfluid(1,1), ldimyz, MPI_INTEGER, neigh(2), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(2), ierr)
    call MPI_ISEND( buf_sendrecv2_isfluid(1,1), ldimyz, MPI_INTEGER, neigh(2), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(2), ierr)

    ! Post sends - Y axis
    buf_sendrecv3_isfluid( 1:ldimxz, 1) = RESHAPE(dtemp(minx:maxx, maxy, minz:maxz), [ ldimxz ] )
    buf_sendrecv4_isfluid( 1:ldimxz, 1) = RESHAPE(dtemp(minx:maxx, miny, minz:maxz), [ ldimxz ] )
    call MPI_ISEND( buf_sendrecv3_isfluid(1,1), ldimxz, MPI_INTEGER, neigh(3), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(3), ierr)
    call MPI_ISEND( buf_sendrecv3_isfluid(1,1), ldimxz, MPI_INTEGER, neigh(3), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(3), ierr)
    call MPI_ISEND( buf_sendrecv4_isfluid(1,1), ldimxz, MPI_INTEGER, neigh(4), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(4), ierr)
    call MPI_ISEND( buf_sendrecv4_isfluid(1,1), ldimxz, MPI_INTEGER, neigh(4), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(4), ierr)

    ! Post sends - Z axis
    buf_sendrecv5_isfluid( 1:ldimxy, 1) = RESHAPE(dtemp(minx:maxx, miny:maxy, maxz), [ ldimxy ] )
    buf_sendrecv6_isfluid( 1:ldimxy, 1) = RESHAPE(dtemp(minx:maxx, miny:maxy, minz), [ ldimxy ] )
    call MPI_ISEND( buf_sendrecv5_isfluid(1,1), ldimxy, MPI_INTEGER, neigh(5), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(5), ierr)
    call MPI_ISEND( buf_sendrecv5_isfluid(1,1), ldimxy, MPI_INTEGER, neigh(5), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(5), ierr)
    call MPI_ISEND( buf_sendrecv6_isfluid(1,1), ldimxy, MPI_INTEGER, neigh(6), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(6), ierr)
    call MPI_ISEND( buf_sendrecv6_isfluid(1,1), ldimxy, MPI_INTEGER, neigh(6), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(6), ierr)

    ! Post sends - z edges
    buf_sendrecv7_isfluid(1:ldims(3), 1) = dtemp(maxx,maxy, minz:maxz)
    buf_sendrecv8_isfluid(1:ldims(3), 1) = dtemp(minx,miny, minz:maxz)
    call MPI_ISEND( buf_sendrecv7_isfluid(1,1), ldims(3), MPI_INTEGER, neigh(7), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(7), ierr)
    call MPI_ISEND( buf_sendrecv7_isfluid(1,1), ldims(3), MPI_INTEGER, neigh(7), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(7), ierr)
    call MPI_ISEND( buf_sendrecv8_isfluid(1,1), ldims(3), MPI_INTEGER, neigh(8), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(8), ierr)
    call MPI_ISEND( buf_sendrecv8_isfluid(1,1), ldims(3), MPI_INTEGER, neigh(8), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(8), ierr)

    buf_sendrecv9_isfluid (1:ldims(3), 1) = dtemp(minx,maxy, minz:maxz)
    buf_sendrecv10_isfluid(1:ldims(3), 1) = dtemp(maxx,miny, minz:maxz)
    call MPI_ISEND( buf_sendrecv9_isfluid(1,1),  ldims(3), MPI_INTEGER, neigh(9),  & 
            id_rank+fldtag, cube_comm, send_req_isfluid(9), ierr)
    call MPI_ISEND( buf_sendrecv9_isfluid(1,1),  ldims(3), MPI_INTEGER, neigh(9),  & 
            id_rank+fldtag, cube_comm, send_req_isfluid(9), ierr)
    call MPI_ISEND( buf_sendrecv10_isfluid(1,1), ldims(3), MPI_INTEGER, neigh(10), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(10), ierr)
    call MPI_ISEND( buf_sendrecv10_isfluid(1,1), ldims(3), MPI_INTEGER, neigh(10), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(10), ierr)

    ! Post sends - y edges
    buf_sendrecv11_isfluid(1:ldims(2), 1) = dtemp(maxx, miny:maxy, maxz)
    buf_sendrecv12_isfluid(1:ldims(2), 1) = dtemp(minx, miny:maxy, minz)
    call MPI_ISEND( buf_sendrecv11_isfluid(1,1), ldims(2), MPI_INTEGER, neigh(11), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(11), ierr)
    call MPI_ISEND( buf_sendrecv11_isfluid(1,1), ldims(2), MPI_INTEGER, neigh(11), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(11), ierr)
    call MPI_ISEND( buf_sendrecv12_isfluid(1,1), ldims(2), MPI_INTEGER, neigh(12), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(12), ierr)
    call MPI_ISEND( buf_sendrecv12_isfluid(1,1), ldims(2), MPI_INTEGER, neigh(12), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(12), ierr)

    buf_sendrecv13_isfluid(1:ldims(2), 1) = dtemp(minx, miny:maxy, maxz)
    buf_sendrecv14_isfluid(1:ldims(2), 1) = dtemp(maxx, miny:maxy, minz)
    call MPI_ISEND( buf_sendrecv13_isfluid(1,1), ldims(2), MPI_INTEGER, neigh(13), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(13), ierr)
    call MPI_ISEND( buf_sendrecv13_isfluid(1,1), ldims(2), MPI_INTEGER, neigh(13), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(13), ierr)
    call MPI_ISEND( buf_sendrecv14_isfluid(1,1), ldims(2), MPI_INTEGER, neigh(14), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(14), ierr)
    call MPI_ISEND( buf_sendrecv14_isfluid(1,1), ldims(2), MPI_INTEGER, neigh(14), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(14), ierr)

    ! Post sends - x edges
    buf_sendrecv15_isfluid(1:ldims(1), 1) = dtemp(minx:maxx, maxy, maxz)
    buf_sendrecv16_isfluid(1:ldims(1), 1) = dtemp(minx:maxx, miny, minz)
    call MPI_ISEND( buf_sendrecv15_isfluid(1,1), ldims(1), MPI_INTEGER, neigh(15), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(15), ierr)
    call MPI_ISEND( buf_sendrecv15_isfluid(1,1), ldims(1), MPI_INTEGER, neigh(15), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(15), ierr)
    call MPI_ISEND( buf_sendrecv16_isfluid(1,1), ldims(1), MPI_INTEGER, neigh(16), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(16), ierr)
    call MPI_ISEND( buf_sendrecv16_isfluid(1,1), ldims(1), MPI_INTEGER, neigh(16), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(16), ierr)

    buf_sendrecv17_isfluid(1:ldims(1), 1) = dtemp(minx:maxx, miny, maxz)
    buf_sendrecv18_isfluid(1:ldims(1), 1) = dtemp(minx:maxx, maxy, minz)
    call MPI_ISEND( buf_sendrecv17_isfluid(1,1), ldims(1), MPI_INTEGER, neigh(17), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(17), ierr)
    call MPI_ISEND( buf_sendrecv17_isfluid(1,1), ldims(1), MPI_INTEGER, neigh(17), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(17), ierr)
    call MPI_ISEND( buf_sendrecv18_isfluid(1,1), ldims(1), MPI_INTEGER, neigh(18), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(18), ierr)
    call MPI_ISEND( buf_sendrecv18_isfluid(1,1), ldims(1), MPI_INTEGER, neigh(18), & 
            id_rank+fldtag, cube_comm, send_req_isfluid(18), ierr)

    end subroutine mpisend_isfluid


    subroutine mpirecv_isfluid(dtemp)
    implicit none
    INTEGER(kind=1), dimension(:,:,:), allocatable :: dtemp
    integer         :: statusArray(MPI_STATUS_SIZE,0:17)

    CALL MPI_WAITALL(18,recv_req_isfluid,statusArray,ierr)
    if (ierr /= 0) write (6,*) "mpirecv_isfluid: procID=", id_rank, "ierr=", ierr

    ! Wait & collect recv
    if (neigh(2) /= MPI_PROC_NULL) then
    dtemp(minx-1, miny:maxy, minz:maxz) = RESHAPE(buf_sendrecv1_isfluid(1:ldimyz, 2), [ ldims(2),ldims(3) ])
    endif

    if (neigh(1) /= MPI_PROC_NULL) then
    dtemp(maxx+1, miny:maxy, minz:maxz) = RESHAPE(buf_sendrecv2_isfluid(1:ldimyz, 2), [ ldims(2),ldims(3) ])
    endif

    if (neigh(4) /= MPI_PROC_NULL) then
    dtemp(minx:maxx, miny-1, minz:maxz) = RESHAPE(buf_sendrecv3_isfluid(1:ldimxz, 2), [ ldims(1),ldims(3) ])
    endif

    if (neigh(3) /= MPI_PROC_NULL) then
    dtemp(minx:maxx, maxy+1, minz:maxz) = RESHAPE(buf_sendrecv4_isfluid(1:ldimxz, 2), [ ldims(1),ldims(3) ])
    endif

    if (neigh(6) /= MPI_PROC_NULL) then
    dtemp(minx:maxx, miny:maxy, minz-1) = RESHAPE(buf_sendrecv5_isfluid(1:ldimxy, 2), [ ldims(1),ldims(2) ])
    endif

    if (neigh(5) /= MPI_PROC_NULL) then
    dtemp(minx:maxx, miny:maxy, maxz+1) = RESHAPE(buf_sendrecv6_isfluid(1:ldimxy, 2), [ ldims(1),ldims(2) ])
    endif

    ! +1,+1 & -1,-1
    if (neigh(8) /= MPI_PROC_NULL) then
    dtemp(minx-1,miny-1, minz:maxz) = buf_sendrecv7_isfluid(:,2)
    endif

    if (neigh(7) /= MPI_PROC_NULL) then
    dtemp(maxx+1,maxy+1, minz:maxz) = buf_sendrecv8_isfluid(:,2)
    endif

    ! -1,+1 & +1,-1
    if (neigh(10) /= MPI_PROC_NULL) then
    dtemp(maxx+1,miny-1, minz:maxz) = buf_sendrecv9_isfluid(:,2)
    endif

    if (neigh(9) /= MPI_PROC_NULL) then
    dtemp(minx-1,maxy+1, minz:maxz) = buf_sendrecv10_isfluid(:,2)
    endif


    ! 1,0,1 & -1,0,-1
    if (neigh(12) /= MPI_PROC_NULL) then
    dtemp(minx-1,miny:maxy, minz-1) = buf_sendrecv11_isfluid(:,2)
    endif

    if (neigh(11) /= MPI_PROC_NULL) then
    dtemp(maxx+1,miny:maxy, maxz+1) = buf_sendrecv12_isfluid(:,2)
    endif

    ! -1,0,1 & 1,0,-1
    if (neigh(14) /= MPI_PROC_NULL) then
    dtemp(maxx+1,miny:maxy, minz-1) = buf_sendrecv13_isfluid(:,2)
    endif

    if (neigh(13) /= MPI_PROC_NULL) then
    dtemp(minx-1,miny:maxy, maxz+1) = buf_sendrecv14_isfluid(:,2)
    endif

    ! 0,1,1 & 0,-1,-1
    if (neigh(16) /= MPI_PROC_NULL) then
    dtemp(minx:maxx, miny-1, minz-1) = buf_sendrecv15_isfluid(:,2)
    endif

    if (neigh(15) /= MPI_PROC_NULL) then
    dtemp(minx:maxx, maxy+1, maxz+1) = buf_sendrecv16_isfluid(:,2)
    endif

    ! 0,-1,1 & 0,1,-1
    if (neigh(18) /= MPI_PROC_NULL) then
    dtemp(minx:maxx, maxy+1, minz-1) = buf_sendrecv17_isfluid(:,2)
    endif

    if (neigh(17) /= MPI_PROC_NULL) then
    dtemp(minx:maxx, miny-1, maxz+1) = buf_sendrecv18_isfluid(:,2)
    endif

    CALL MPI_WAITALL(18,send_req_isfluid,statusArray,ierr)
    call mpipostrecv_isfluid

    end subroutine mpirecv_isfluid


    subroutine mpisendrecvhalopops(aoptp)
    implicit none
    type(REALPTR), dimension(0:npop), intent(in)   :: aoptp
    integer         :: statusVar(MPI_STATUS_SIZE), i


    ! X axis
    do i=0,18
        buf_sendrecv1( i*ldimyz+1 : (i+1)*ldimyz, 1) = RESHAPE(aoptp(i)%p(maxx, miny:maxy,minz:maxz), [ ldimyz ] )
        buf_sendrecv2( i*ldimyz+1 : (i+1)*ldimyz, 1) = RESHAPE(aoptp(i)%p(minx, miny:maxy,minz:maxz), [ ldimyz ] )
    enddo

    call MPI_SENDRECV(  buf_sendrecv1(1,1), 19*ldimyz, MPI_DOUBLE, neigh(1), id_rank+halotag,   &
                        buf_sendrecv1(1,2), 19*ldimyz, MPI_DOUBLE, neigh(2), neigh(2)+halotag, &
        cube_comm, statusVar, ierr)
    call MPI_SENDRECV( buf_sendrecv2(1,1), 19*ldimyz, MPI_DOUBLE, neigh(2), id_rank+halotag,  &
                       buf_sendrecv2(1,2), 19*ldimyz, MPI_DOUBLE, neigh(1), neigh(1)+halotag, &
        cube_comm, statusVar, ierr)

    if (neigh(2) /= MPI_PROC_NULL) then
      do i=0,18
        aoptp(i)%p(minx-1, miny:maxy,minz:maxz) = RESHAPE(buf_sendrecv1(i*ldimyz+1:(i+1)*ldimyz, 2), [ ldims(2),ldims(3) ])
      enddo
    endif

    if (neigh(1) /= MPI_PROC_NULL) then
      do i=0,18
        aoptp(i)%p(maxx+1,miny:maxy,minz:maxz) = RESHAPE(buf_sendrecv2(i*ldimyz+1:(i+1)*ldimyz, 2), [ ldims(2),ldims(3) ])
      enddo
    endif

    ! Z axis
    do i=0,18
        buf_sendrecv5( i*ldimxy+1:(i+1)*ldimxy, 1) = RESHAPE(aoptp(i)%p(minx:maxx, miny:maxy, maxz), [ ldimxy ] )
        buf_sendrecv6( i*ldimxy+1:(i+1)*ldimxy, 1) = RESHAPE(aoptp(i)%p(minx:maxx, miny:maxy, minz), [ ldimxy ] )
    enddo

    call MPI_SENDRECV(  buf_sendrecv5(1,1), 19*ldimxy, MPI_DOUBLE, neigh(5), id_rank+halotag,   &
                        buf_sendrecv5(1,2), 19*ldimxy, MPI_DOUBLE, neigh(6), neigh(6)+halotag, &
        cube_comm, statusVar, ierr)
    call MPI_SENDRECV( buf_sendrecv6(1,1), 19*ldimxy, MPI_DOUBLE, neigh(6), id_rank+halotag,  &
                       buf_sendrecv6(1,2), 19*ldimxy, MPI_DOUBLE, neigh(5), neigh(5)+halotag, &
        cube_comm, statusVar, ierr)

    if (neigh(6) /= MPI_PROC_NULL) then
      do i=0,18
        aoptp(i)%p(minx:maxx, miny:maxy, minz-1) = RESHAPE(buf_sendrecv5(i*ldimxy+1:(i+1)*ldimxy, 2), [ ldims(1),ldims(2) ])
      enddo
    endif

    if (neigh(5) /= MPI_PROC_NULL) then
      do i=0,18
        aoptp(i)%p(minx:maxx, miny:maxy, maxz+1) = RESHAPE(buf_sendrecv6(i*ldimxy+1:(i+1)*ldimxy, 2), [ ldims(1),ldims(2) ])
      enddo
    endif

    end subroutine mpisendrecvhalopops

end module mpi_comm
