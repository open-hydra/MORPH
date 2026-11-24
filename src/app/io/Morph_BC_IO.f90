!>@brief Reader of the MORPH boundary-condition file written by ATLAS: one
!> record per boundary cell, coordinate line `block i j k face id` (6I8; 5I8
!> means k = 1), then the data lines of the code. Fills domain%bfaces and
!> domain%conns. Driver code: the library does no I/O.
!>
!> On a bad file it returns ok = .false. and a message with the file and the
!> line; a driver that reads it stops with exit status 1 (input).
module Morph_BC_IO
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_API

  implicit none
  private

  public :: Read_Morph_BC, Face_Local

  integer(I4), parameter :: UNSET = -1

contains

  !> Face-local cell indices (m, n) and extents (n1, n2) of a cell on face f.
  pure subroutine Face_Local ( f, dim, ijk, m, n, n1, n2 )
    integer(I4), intent(in)  :: f, dim(3), ijk(3)
    integer(I4), intent(out) :: m, n, n1, n2
    select case ( f )
    case ( 1, 2 ) ; m = ijk(2) ; n = ijk(3) ; n1 = dim(2) ; n2 = dim(3)
    case ( 3, 4 ) ; m = ijk(1) ; n = ijk(3) ; n1 = dim(1) ; n2 = dim(3)
    case default  ; m = ijk(1) ; n = ijk(2) ; n1 = dim(1) ; n2 = dim(2)
    end select
  end subroutine Face_Local


  subroutine Read_Morph_BC ( filename, domain, msg, ok )
    character(len=*),     intent(in)    :: filename
    type(morph_domain_t), intent(inout) :: domain
    character(len=*),     intent(out)   :: msg
    logical,              intent(out)   :: ok
    ! Local
    integer(I4) :: u, ios, nl, b, f, id, ijk(3), m, n, n1, n2, ic, nc, k2, mat, d(9)
    type :: face_index_t                                ! connection record of each cell face of one block face
      integer(I4), allocatable :: c(:,:)
    end type face_index_t
    type(face_index_t), allocatable :: cidx(:,:)        ! (nb, 6), sized face by face like domain%bfaces
    type(morph_conn_t), allocatable :: tmp(:), conns(:)
    character(len=512) :: line
    real(R8)    :: rate
    logical     :: two_d

    ok = .false. ; msg = ''
    if ( domain%nblocks < 1 ) then
      msg = 'Read_Morph_BC: the domain has no blocks' ; return
    endif
    call alloc_faces ()
    allocate ( conns(64) ) ; nc = 0

    open ( newunit = u, file = filename, status = 'old', action = 'read', iostat = ios )
    if ( ios /= 0 ) then
      msg = 'cannot open '//trim(filename) ; return
    endif
    nl = 0
    do
      read ( u, '(A)', iostat = ios ) line
      if ( ios /= 0 ) exit
      nl = nl + 1
      if ( len_trim(line) == 0 ) cycle
      read ( line, *, iostat = ios ) b, ijk(1), ijk(2), ijk(3), f, id
      two_d = .false.
      if ( ios /= 0 ) then
        read ( line, *, iostat = ios ) b, ijk(1), ijk(2), f, id
        ijk(3) = 1 ; two_d = .true.
      endif
      if ( ios /= 0 ) then
        call bad ( 'not a coordinate line' ) ; return
      endif
      if ( b < 1 .or. b > domain%nblocks .or. f < 1 .or. f > 6 ) then
        call bad ( 'block or face out of range' ) ; return
      endif
      if ( any( ijk < 1 ) .or. any( ijk > domain%blk(b)%dim ) .or. .not. on_face(b, f, ijk) ) then
        call bad ( 'cell indices out of the face' ) ; return
      endif
      call Face_Local ( f, domain%blk(b)%dim, ijk, m, n, n1, n2 )
      if ( domain%bfaces(b,f)%kind(m,n) /= UNSET ) then
        call bad ( 'two records for one cell face' ) ; return
      endif
      select case ( id )
      case ( 0 )
        domain%bfaces(b,f)%kind(m,n) = MORPH_BND_OTHER
      case ( 200, 300 )
        domain%bfaces(b,f)%kind(m,n) = MORPH_BND_PLANE
      case ( 102 )
        domain%bfaces(b,f)%kind(m,n) = MORPH_BND_CHIMERA
      case ( 703 )
        domain%bfaces(b,f)%kind(m,n) = MORPH_BND_FIXEDSURF
      case ( 101, 201 )
        call next_line ( ) ; if ( ios /= 0 ) return
        if ( two_d ) then
          read ( line, *, iostat = ios ) d(1:3), d(5:9) ; d(4) = 1
        else
          read ( line, *, iostat = ios ) d
        endif
        if ( ios /= 0 ) then
          call bad ( 'connection data: b2 i2 j2 k2 f2 c1 c2 c3 c4 expected' ) ; return
        endif
        if ( id == 201 ) then
          domain%bfaces(b,f)%kind(m,n) = MORPH_BND_PERIODIC
        else
          domain%bfaces(b,f)%kind(m,n) = MORPH_BND_CONN
          if ( nc == size(conns) ) then
            allocate ( tmp(2*nc) ) ; tmp(1:nc) = conns(1:nc) ; call move_alloc ( tmp, conns )
          endif
          nc = nc + 1
          conns(nc)%b = b ; conns(nc)%face = f ; conns(nc)%ijk = ijk
          conns(nc)%b2 = d(1) ; conns(nc)%ijk2 = d(2:4) ; conns(nc)%face2 = d(5)
          ! column by column, as MOSE reads them (Lib_Ghost.f90:536-583): (c1, c2) = image of our axis 1
          conns(nc)%C(1,1) = d(6) ; conns(nc)%C(2,1) = d(7) ; conns(nc)%C(1,2) = d(8) ; conns(nc)%C(2,2) = d(9)
          cidx(b,f)%c(m,n) = nc
        endif
      case ( 701 )
        call next_line ( ) ; if ( ios /= 0 ) return
        read ( line, *, iostat = ios ) mat
        if ( ios /= 0 .or. mat < 1 ) then
          call bad ( '701: a material >= 1 expected' ) ; return
        endif
        domain%bfaces(b,f)%kind(m,n) = MORPH_BND_RECEDE ; domain%bfaces(b,f)%mat(m,n) = mat
      case ( 702 )
        call next_line ( ) ; if ( ios /= 0 ) return
        read ( line, *, iostat = ios ) mat, rate
        if ( ios /= 0 .or. mat < 1 .or. rate < 0.0_R8 ) then
          call bad ( '702: a material >= 1 and a rate >= 0 expected' ) ; return
        endif
        domain%bfaces(b,f)%kind(m,n) = MORPH_BND_IMPOSED
        domain%bfaces(b,f)%mat(m,n) = mat ; domain%bfaces(b,f)%rate(m,n) = rate
      case default
        call bad ( 'unknown code' ) ; return
      end select
    enddo
    close ( u )

    ! Every boundary cell face has exactly one record
    do b = 1, domain%nblocks
      do f = 1, 6
        if ( any( domain%bfaces(b,f)%kind == UNSET ) ) then
          write ( msg, '(A,A,I0,A,I0)' ) trim(filename), ': a boundary cell without a record, block ', b, ' face ', f
          return
        endif
      enddo
    enddo

    ! The two records of a connection point at each other
    do ic = 1, nc
      associate ( c => conns(ic) )
        if ( c%b2 < 1 .or. c%b2 > domain%nblocks .or. c%face2 < 1 .or. c%face2 > 6 ) then
          write ( msg, '(A,A,I0)' ) trim(filename), ': connection to a block or face out of range, record ', ic
          return
        endif
        if ( any( c%ijk2 < 1 ) .or. any( c%ijk2 > domain%blk(c%b2)%dim ) ) then
          write ( msg, '(A,A,I0)' ) trim(filename), ': connection to a cell out of range, record ', ic
          return
        endif
        call Face_Local ( c%face2, domain%blk(c%b2)%dim, c%ijk2, m, n, n1, n2 )
        k2 = cidx(c%b2,c%face2)%c(m,n)
        if ( k2 == 0 ) then
          write ( msg, '(A,A,I0)' ) trim(filename), ': connection without the reverse record, record ', ic
          return
        endif
        if ( conns(k2)%b2 /= c%b .or. conns(k2)%face2 /= c%face .or. any( conns(k2)%ijk2 /= c%ijk ) ) then
          write ( msg, '(A,A,I0)' ) trim(filename), &
            ': the two records of a connection do not point at each other, record ', ic
          return
        endif
      end associate
    enddo

    if ( allocated(domain%conns) ) deallocate ( domain%conns )
    allocate ( domain%conns(nc) ) ; domain%conns = conns(1:nc)
    ok = .true.

  contains

    subroutine alloc_faces ()
      integer(I4) :: bb, ff, a1, a2
      if ( allocated(domain%bfaces) ) deallocate ( domain%bfaces )
      allocate ( domain%bfaces(domain%nblocks, 6), cidx(domain%nblocks, 6) )
      do bb = 1, domain%nblocks
        do ff = 1, 6
          call Face_Local ( ff, domain%blk(bb)%dim, [1, 1, 1], m, n, a1, a2 )
          allocate ( domain%bfaces(bb,ff)%kind(a1,a2), domain%bfaces(bb,ff)%mat(a1,a2), &
                     domain%bfaces(bb,ff)%rate(a1,a2) )
          domain%bfaces(bb,ff)%kind = UNSET
          domain%bfaces(bb,ff)%mat  = 0
          domain%bfaces(bb,ff)%rate = 0.0_R8
          allocate ( cidx(bb,ff)%c(a1,a2) ) ; cidx(bb,ff)%c = 0
        enddo
      enddo
    end subroutine alloc_faces

    logical function on_face ( bb, ff, c ) result ( res )
      integer(I4), intent(in) :: bb, ff, c(3)
      select case ( ff )
      case ( 1 ) ; res = ( c(1) == 1 )
      case ( 2 ) ; res = ( c(1) == domain%blk(bb)%dim(1) )
      case ( 3 ) ; res = ( c(2) == 1 )
      case ( 4 ) ; res = ( c(2) == domain%blk(bb)%dim(2) )
      case ( 5 ) ; res = ( c(3) == 1 )
      case default ; res = ( c(3) == domain%blk(bb)%dim(3) )
      end select
    end function on_face

    subroutine next_line ()
      read ( u, '(A)', iostat = ios ) line
      nl = nl + 1
      if ( ios /= 0 ) call bad ( 'data line missing' )
    end subroutine next_line

    subroutine bad ( what )
      character(len=*), intent(in) :: what
      write ( msg, '(A,A,I0,A,A)' ) trim(filename), ':', nl, ': ', what
      close ( u )
    end subroutine bad

  end subroutine Read_Morph_BC

end module Morph_BC_IO
