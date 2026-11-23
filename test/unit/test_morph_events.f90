!>@brief Test of the event list (Morph_Events).
!>
!> It pins: an event with a code below MORPH_EV_STOP_LIMIT stops and one at
!> the limit or above is reported only (both sides of the boundary); `count`
!> defaults to 1; the list grows past its first allocation; `clear` empties it.
program test_morph_events
  use iso_fortran_env, only: I4 => int32
  use Morph_API
  implicit none
  type(morph_event_list_t) :: ev
  integer(I4) :: n
  logical     :: ok

  ok = .true.
  write(*,'(A)') '== Morph event-list test =='
  write(*,*)

  call ev%add ( MORPH_EV_THIN, 'thin' )
  call expect ( ev%n == 1 .and. .not. ev%any_stop(), 'a report-only event does not stop' )
  call expect ( ev%ev(1)%count == 1, 'count defaults to 1' )

  call ev%clear ()
  call ev%add ( MORPH_EV_STOP_LIMIT, 'at the limit' )
  call expect ( .not. ev%ev(1)%stops .and. .not. ev%any_stop(), 'a code equal to the limit does not stop' )
  call ev%add ( MORPH_EV_STOP_LIMIT - 1, 'below the limit' )
  call expect ( ev%ev(2)%stops .and. ev%any_stop(), 'a code just below the limit stops' )

  call ev%clear ()
  call ev%add ( MORPH_EV_THIN, 'thin' )
  do n = 1, 20
    call ev%add ( MORPH_EV_MATERIAL, 'material', count = n )
  enddo
  call expect ( ev%n == 21 .and. ev%ev(21)%count == 20 .and. ev%ev(21)%code == MORPH_EV_MATERIAL, &
                'the list grows past its first allocation' )
  call expect ( ev%any_stop(), 'a stopping event among report-only ones is seen' )

  call ev%clear ()
  call expect ( ev%n == 0 .and. .not. ev%any_stop(), 'clear empties the list' )

  write(*,*)
  if ( ok ) then
    write(*,'(A)') '   RESULT: PASS'
  else
    write(*,'(A)') '   RESULT: FAIL'
    error stop 1
  endif

contains

  subroutine expect ( good, what )
    logical,          intent(in) :: good
    character(len=*), intent(in) :: what
    if ( good ) then
      write(*,'(A,A,A)') '   ', what, ': PASS'
    else
      write(*,'(A,A)') '   FAIL: ', what
      ok = .false.
    endif
  end subroutine expect

end program test_morph_events
