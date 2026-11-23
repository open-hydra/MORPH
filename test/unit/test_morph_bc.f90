!>@brief Test of the reader of the MORPH boundary-condition file (the driver's,
!>       src/app/io/Morph_BC_IO.f90).
!>
!> Two blocks of 2x1x1 cells joined by face 2 of block 1 and face 1 of block
!> 2; one valid file with every code, then five broken ones that must be
!> refused with their own message. The orientation map of a connection is read
!> column by column (c1 c2 = the image of axis 1), checked on a map that is not
!> symmetric.
program test_morph_bc
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_API
  use Morph_BC_IO, only: Read_Morph_BC
  implicit none

  integer(I4) :: nfail = 0
  character(len=512) :: msg
  logical :: ok
  type(morph_domain_t) :: dom
  type(morph_status_t) :: st

  write(*,'(A)') '== Morph boundary-file reader test =='
  write(*,*)

  call new_domain ()
  call write_file ( 'bc-ok.txt', 'none' )
  call Read_Morph_BC ( 'bc-ok.txt', dom, msg, ok )
  call check ( ok, 'valid file read'//note(msg) )
  if ( ok ) then
    call check ( size(dom%conns) == 2, 'two connection records' )
    call check ( all( dom%bfaces(1,3)%kind == MORPH_BND_RECEDE ) .and. all( dom%bfaces(1,3)%mat == 1 ), &
                 '701 -> recede, material 1' )
    call check ( all( dom%bfaces(2,4)%kind == MORPH_BND_IMPOSED ) .and. all( dom%bfaces(2,4)%rate == 1.5e-3_R8 ), &
                 '702 -> imposed, rate' )
    call check ( all( dom%bfaces(1,5)%kind == MORPH_BND_PLANE ), '300 -> plane' )
    call check ( dom%bfaces(1,2)%kind(1,1) == MORPH_BND_CONN .and. dom%conns(1)%C(1,1) == 0 .and. &
                 dom%conns(1)%C(2,1) == 1 .and. dom%conns(1)%C(1,2) == -1 .and. dom%conns(1)%C(2,2) == 0, &
                 '101 -> connection, map read column by column (c1 c2 = image of axis 1)' )
  endif

  call refused ( 'missing',   'without a record' )
  call refused ( 'duplicate', 'two records' )
  call refused ( 'wrongconn', 'point at each other' )
  call refused ( 'unknown',   'unknown code' )
  call refused ( 'offface',   'out of the face' )

  call dom%free ()

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

  subroutine new_domain ()
    integer(I4) :: dims(3,2)
    dims(:,1) = [2, 1, 1] ; dims(:,2) = [2, 1, 1]
    call dom%allocate_blocks ( dims, st )
  end subroutine new_domain

  subroutine refused ( kind, expect )
    character(len=*), intent(in) :: kind, expect
    call write_file ( 'bc-'//kind//'.txt', kind )
    call Read_Morph_BC ( 'bc-'//kind//'.txt', dom, msg, ok )
    call check ( .not. ok .and. index(msg, expect) > 0, 'refused ('//kind//'): '//trim(msg) )
  end subroutine refused

  !> The records of both blocks; `mut` breaks the file in one way.
  subroutine write_file ( name, mut )
    character(len=*), intent(in) :: name, mut
    integer(I4) :: u, b, f, i, j, k, code
    logical :: done_dup
    open ( newunit = u, file = name, status = 'replace', action = 'write' )
    done_dup = .false.
    do b = 1, 2
      do f = 1, 6
        do k = 1, 1; do j = 1, 1; do i = 1, 2
          if ( .not. on_face(f, i, j, k) ) cycle
          code = code_of(b, f)
          if ( mut == 'missing' .and. b == 2 .and. f == 6 .and. i == 2 ) cycle
          if ( mut == 'unknown' .and. b == 1 .and. f == 5 .and. i == 1 ) code = 999
          if ( mut == 'offface' .and. b == 1 .and. f == 2 ) then
            write(u,'(6I8)') b, 1, j, k, f, code       ! i = 1 is not on face 2 of a 2-cell block
          else
            write(u,'(6I8)') b, i, j, k, f, code
          endif
          select case ( code )
          case ( 101 )
            if ( b == 1 ) then
              if ( mut == 'wrongconn' ) then
                write(u,'(9I8)') 2, 2, 1, 1, 1, 1, 0, 0, 1   ! cell 2 of block 2 is not on its face 1
              else
                write(u,'(9I8)') 2, 1, 1, 1, 1, 0, 1, -1, 0   ! a map that is not symmetric: column order checked
              endif
            else
              write(u,'(9I8)') 1, 2, 1, 1, 2, 0, -1, 1, 0
            endif
          case ( 701 ) ; write(u,'(I8)') 1
          case ( 702 )
            if ( b == 1 ) then
              write(u,'(I8,ES24.16)') 2, 0.0_R8
            else
              write(u,'(I8,ES24.16)') 1, 1.5e-3_R8
            endif
          end select
          if ( mut == 'duplicate' .and. .not. done_dup .and. code == 300 ) then
            write(u,'(6I8)') b, i, j, k, f, code ; done_dup = .true.
          endif
        enddo; enddo; enddo
      enddo
    enddo
    close ( u )
  end subroutine write_file

  integer(I4) function code_of ( b, f ) result ( c )
    integer(I4), intent(in) :: b, f
    select case ( f )
    case ( 1 ) ; c = merge( 703, 101, b == 1 )
    case ( 2 ) ; c = merge( 101, 703, b == 1 )
    case ( 3 ) ; c = 701
    case ( 4 ) ; c = 702
    case default ; c = 300
    end select
  end function code_of

  logical function on_face ( f, i, j, k ) result ( res )
    integer(I4), intent(in) :: f, i, j, k
    select case ( f )
    case ( 1 ) ; res = ( i == 1 )
    case ( 2 ) ; res = ( i == 2 )
    case ( 3, 4 ) ; res = ( j == 1 )
    case default ; res = ( k == 1 )
    end select
  end function on_face

end program test_morph_bc
