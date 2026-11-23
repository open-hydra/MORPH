!>@brief Test of the surface layer: Morph_Surface_Build and Morph_Fixed_Masks
!>       on small meshes, each read from a MORPH boundary file written here
!>       (the driver's reader, src/app/io/Morph_BC_IO.f90).
!>
!> The expected numbers are counted by hand from the geometry below.
!>  1. Two blocks of 2x2x2 cells, x in [0,1] and [1,2], joined at x = 1 by an
!>     interface rotated by 90 degrees (block 2: i -> +x, j -> -z, k -> +y), the
!>     wall on z = 0 across both, the plane y = 0 across both. Through the file:
!>     writer -> Read_Morph_BC -> Morph_Surface_Build.
!>  2. The same file with the orientation map written row by row (transposed):
!>     the build must refuse it.
!>  3. One block whose faces i0 and j0 are two planes meeting at 0.5 degrees
!>     along their common edge: two planes, not one.
!>  4. Fixed masks: the nodes of block 1 at x = 1 fixed; their copies in block 2.
program test_morph_surface
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_API
  use Morph_BC_IO, only: Read_Morph_BC
  implicit none

  integer(I4) :: nfail = 0
  real(R8), parameter :: H = 0.5_R8

  write(*,'(A)') '== Morph surface test =='
  write(*,*)

  call test_two_blocks ( .false. )
  call test_two_blocks ( .true. )
  call test_wedge_planes ()

  write(*,*)
  if ( nfail == 0 ) then
    write(*,'(A)') '   RESULT: PASS'
  else
    write(*,'(A)') '   RESULT: FAIL'
    error stop 1
  endif

contains

  subroutine check ( good, what )
    logical,          intent(in) :: good
    character(len=*), intent(in) :: what
    if ( good ) then
      write(*,'(A,A,A)') '   ', what, ': PASS'
    else
      write(*,'(A,A)') '   FAIL: ', what
      nfail = nfail + 1
    endif
  end subroutine check

  !> ' (text)' when there is a text, nothing when there is none.
  function note ( text ) result ( s )
    character(len=*), intent(in)  :: text
    character(len=:), allocatable :: s
    if ( len_trim(text) > 0 ) then
      s = ' ('//trim(text)//')'
    else
      s = ''
    endif
  end function note

  !> code of face f of block b in the two-block case
  integer(I4) function code2 ( b, f ) result ( c )
    integer(I4), intent(in) :: b, f
    integer(I4), parameter :: CODE_A(6) = [703, 101, 300, 703, 701, 703]
    integer(I4), parameter :: CODE_B(6) = [101, 703, 703, 701, 300, 703]
    if ( b == 1 ) then
      c = CODE_A(f)
    else
      c = CODE_B(f)
    endif
  end function code2

  subroutine test_two_blocks ( transposed )
    logical, intent(in) :: transposed
    type(morph_domain_t) :: dom
    type(morph_status_t) :: st
    type(morph_surface_stats_t) :: ss
    character(len=256) :: msg
    logical :: ok, ccw, owners
    integer(I4) :: dims(3,2), b, i, j, k, u, g, nkind(0:7), iv(7,1)
    real(R8) :: x(3,4), nrm(3)
    dims(:,1) = [2, 2, 2] ; dims(:,2) = [2, 2, 2]
    call dom%allocate_blocks ( dims, st )
    do k = 0, 2 ; do j = 0, 2 ; do i = 0, 2
      dom%blk(1)%node(i,j,k)%c = [ H*i, H*j, H*k ]
      dom%blk(2)%node(i,j,k)%c = [ 1.0_R8 + H*i, H*k, H*(2-j) ]
    enddo ; enddo ; enddo
    call write_two_blocks ( 'surf2.txt', transposed )
    call Read_Morph_BC ( 'surf2.txt', dom, msg, ok )
    call check ( ok, 'two blocks: file read'//note(msg) )
    if ( .not. ok ) return
    call Morph_Surface_Build ( dom, 1.0e-9_R8, st, ss )
    if ( transposed ) then
      call check ( st%code == MORPH_ERR_INPUT, 'transposed orientation map refused: '//trim(st%message) )
      call dom%free ()
      return
    endif
    call check ( st%ok(), 'two blocks: build ok'//note(st%message) )
    if ( .not. st%ok() ) return
    associate ( s => dom%surf )
      call check ( s%nface == 8 .and. s%nnode == 15 .and. s%first(s%nnode+1) - 1 == 18, &
                   'wall: 8 faces, 15 nodes, 18 copies (3 nodes at x = 1 with 2 copies)' )
      owners = .true.
      do g = 1, s%nnode
        if ( s%first(g+1) - s%first(g) == 2 .and. s%copy_b(s%owner(g)) /= 1 ) owners = .false.
      enddo
      call check ( owners, 'owner = first copy met (block 1)' )
      call check ( ss%edges_by_faces(1) == 12 .and. ss%edges_by_faces(2) == 10, 'edges: 12 boundary, 10 interior' )
      call check ( ss%nplane == 1 .and. abs(abs(dom%planes(1)%n(2)) - 1.0_R8) < 1.0e-12_R8 .and. &
                   abs(dom%planes(1)%p0(2)) < 1.0e-15_R8, 'one plane y = 0 across the two blocks' )
      nkind = 0
      do u = 1, s%nedge
        nkind(s%edge_kind(u)) = nkind(s%edge_kind(u)) + 1
      enddo
      call check ( nkind(MORPH_BND_PLANE) == 4 .and. nkind(MORPH_BND_FIXEDSURF) == 8, &
                   'boundary edges: 4 on the plane, 8 on fixed surfaces' )
      call check ( ss%wall_nodes_on_plane(1) == 5, 'wall nodes on the plane: 5' )
      ccw = .true.
      do u = 1, s%nface
        do i = 1, 4
          g = s%face_node(i,u)
          b = s%copy_b(s%owner(g))
          x(:,i) = dom%blk(b)%node(s%copy_ijk(1,s%owner(g)), s%copy_ijk(2,s%owner(g)), s%copy_ijk(3,s%owner(g)))%c
        enddo
        nrm = [ (x(2,3)-x(2,1))*(x(3,4)-x(3,2)) - (x(3,3)-x(3,1))*(x(2,4)-x(2,2)), &
                (x(3,3)-x(3,1))*(x(1,4)-x(1,2)) - (x(1,3)-x(1,1))*(x(3,4)-x(3,2)), &
                (x(1,3)-x(1,1))*(x(2,4)-x(2,2)) - (x(2,3)-x(2,1))*(x(1,4)-x(1,2)) ]
        if ( nrm(3) <= 0.0_R8 ) ccw = .false.   ! the fluid is above z = 0
      enddo
      call check ( ccw, 'wall faces counter-clockwise seen from the fluid' )
    end associate
    ! fixed masks: block 1 at x = 1 (i = 2), the whole face
    iv(:,1) = [ 1, 2, 2, -1, -1, -1, -1 ]
    call Morph_Fixed_Masks ( dom, iv, st )
    call check ( st%ok() .and. count(dom%fixed(1)%m) == 9 .and. count(dom%fixed(2)%m) == 9 .and. &
                 all( dom%fixed(2)%m(0,:,:) ), 'masks: 9 nodes fixed in block 1, their 9 copies in block 2' )
    call dom%free ()
  end subroutine test_two_blocks

  subroutine write_two_blocks ( name, transposed )
    character(len=*), intent(in) :: name
    logical,          intent(in) :: transposed
    integer(I4) :: u, b, f, c, m, n, ijk(3)
    open ( newunit = u, file = name, status = 'replace', action = 'write' )
    do b = 1, 2
      do f = 1, 6
        c = code2(b, f)
        do n = 1, 2
          do m = 1, 2
            select case ( f )
            case ( 1 ) ; ijk = [ 1, m, n ]
            case ( 2 ) ; ijk = [ 2, m, n ]
            case ( 3 ) ; ijk = [ m, 1, n ]
            case ( 4 ) ; ijk = [ m, 2, n ]
            case ( 5 ) ; ijk = [ m, n, 1 ]
            case default ; ijk = [ m, n, 2 ]
            end select
            write(u,'(6I8)') b, ijk, f, c
            if ( c == 101 .and. b == 1 ) then          ! A (2, j, k) -> B (1, 3-k, j); C_A columns (0,1) (-1,0)
              if ( transposed ) then
                write(u,'(9I8)') 2, 1, 3-ijk(3), ijk(2), 1, 0, -1, 1, 0
              else
                write(u,'(9I8)') 2, 1, 3-ijk(3), ijk(2), 1, 0, 1, -1, 0
              endif
            elseif ( c == 101 ) then                   ! B (1, jB, kB) -> A (2, kB, 3-jB); C_B = transpose of C_A
              if ( transposed ) then
                write(u,'(9I8)') 1, 2, ijk(3), 3-ijk(2), 2, 0, 1, -1, 0
              else
                write(u,'(9I8)') 1, 2, ijk(3), 3-ijk(2), 2, 0, -1, 1, 0
              endif
            elseif ( c == 701 ) then
              write(u,'(I8)') 1
            endif
          enddo
        enddo
      enddo
    enddo
    close ( u )
  end subroutine write_two_blocks

  !> One block: face j0 on y = 0, face i0 on the plane through the z axis at
  !> 0.5 degrees from it; they share the edge i = j = 0.
  subroutine test_wedge_planes ()
    type(morph_domain_t) :: dom
    type(morph_status_t) :: st
    type(morph_surface_stats_t) :: ss
    character(len=256) :: msg
    logical :: ok
    integer(I4) :: dims(3,1), i, j, k, u, f, c, m, n, ijk(3), nkind(0:7)
    real(R8), parameter :: TH = 0.5_R8 * 3.14159265358979323846_R8 / 180.0_R8
    integer(I4), parameter :: CODE_W(6) = [300, 703, 300, 703, 703, 701]
    dims(:,1) = [2, 2, 2]
    call dom%allocate_blocks ( dims, st )
    do k = 0, 2 ; do j = 0, 2 ; do i = 0, 2
      dom%blk(1)%node(i,j,k)%c = [ H*i + H*j*cos(TH), H*j*sin(TH), H*k ]
    enddo ; enddo ; enddo
    open ( newunit = u, file = 'surf1.txt', status = 'replace', action = 'write' )
    do f = 1, 6
      c = CODE_W(f)
      do n = 1, 2
        do m = 1, 2
          select case ( f )
          case ( 1 ) ; ijk = [ 1, m, n ]
          case ( 2 ) ; ijk = [ 2, m, n ]
          case ( 3 ) ; ijk = [ m, 1, n ]
          case ( 4 ) ; ijk = [ m, 2, n ]
          case ( 5 ) ; ijk = [ m, n, 1 ]
          case default ; ijk = [ m, n, 2 ]
          end select
          write(u,'(6I8)') 1, ijk, f, c
          if ( c == 701 ) write(u,'(I8)') 1
        enddo
      enddo
    enddo
    close ( u )
    call Read_Morph_BC ( 'surf1.txt', dom, msg, ok )
    call Morph_Surface_Build ( dom, 1.0e-9_R8, st, ss )
    call check ( ok .and. st%ok(), 'wedge 0.5 degrees: build ok'//note(st%message) )
    if ( .not. st%ok() ) return
    call check ( ss%nplane == 2, 'two planes meeting at 0.5 degrees stay two planes' )
    call check ( abs(ss%min_split_deg - 0.5_R8) < 1.0e-9_R8 .and. ss%max_joined_deg < 1.0e-9_R8, &
                 'diagnostics: 0.5 degrees between the planes, 0 inside each' )
    nkind = 0
    do u = 1, dom%surf%nedge
      nkind(dom%surf%edge_kind(u)) = nkind(dom%surf%edge_kind(u)) + 1
    enddo
    call check ( nkind(MORPH_BND_PLANE) == 4 .and. nkind(MORPH_BND_FIXEDSURF) == 4, &
                 'wall boundary: 2 + 2 edges on the two planes, 4 on fixed surfaces' )
    call dom%free ()
  end subroutine test_wedge_planes

end program test_morph_surface
