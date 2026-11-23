!>@brief Events: what happened during a call that the caller must know about.
!>
!> An event is not an error: the status stays MORPH_OK. Events with a code
!> below MORPH_EV_STOP_LIMIT stop the call (the regression returns the fraction
!> of the step it applied); the others are reported only. Each code is
!> described where it is defined. In base/ because the level set, the volume
!> law and the nodal engine all emit them.
!>
!> MORPH must never `use` a FUSS module -- see Morph_Types_m.
module Morph_Events
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m, only: MORPH_MSGLEN

  implicit none
  private

  integer(I4), parameter, public :: MORPH_EV_FACE_INERT     = 1   !> faces became inert (case contact, end of burn)
  integer(I4), parameter, public :: MORPH_EV_SPAN_BLOCKEDGE = 2   !> a collapsed span is a block edge
  integer(I4), parameter, public :: MORPH_EV_PATCH_CONSUMED = 3   !> contact node without a root
  integer(I4), parameter, public :: MORPH_EV_FRONTS_MERGE   = 4   !> two walls become one
  integer(I4), parameter, public :: MORPH_EV_LINE_CONSUMED  = 5   !> solid line below the minimum thickness
  integer(I4), parameter, public :: MORPH_EV_BREAKTHROUGH   = 6   !> wall reaches a fixed interface in front
  integer(I4), parameter, public :: MORPH_EV_STEP_MIDWALL   = 7   !> step born in the middle of the wall
  integer(I4), parameter, public :: MORPH_EV_QUALITY        = 8   !> scaled Jacobian below the threshold
  integer(I4), parameter, public :: MORPH_EV_MATERIAL       = 9   !> buried material uncovered (follow-up: reserved)
  integer(I4), parameter, public :: MORPH_EV_THIN           = 101 !> thickness below the level-set resolution (report only)
  integer(I4), parameter, public :: MORPH_EV_ANCHOR_BORN    = 102 !> a convex corner born from a fillet (report only)

  !> Codes below this stop the call.
  integer(I4), parameter, public :: MORPH_EV_STOP_LIMIT = 100

  type, public :: morph_event_t
    integer(I4) :: code = 0
    logical     :: stops = .false.
    integer(I4) :: b = 0, i = 0, j = 0, k = 0, face = 0, count = 0
    real(R8)    :: stroke = 0.0_R8              !> accumulated stroke when it happened
    character(len=MORPH_MSGLEN) :: message = ''
  end type morph_event_t

  type, public :: morph_event_list_t
    integer(I4) :: n = 0
    type(morph_event_t), allocatable :: ev(:)
  contains
    procedure :: add      => events_add
    procedure :: clear    => events_clear
    procedure :: any_stop => events_any_stop
  end type morph_event_list_t

contains

  !> Append one event. Faces of one sub-step that share the code are given as
  !> one event, with their number in `count`.
  pure subroutine events_add ( self, code, message, b, i, j, k, face, count, stroke )
    class(morph_event_list_t), intent(inout) :: self
    integer(I4),               intent(in)    :: code
    character(len=*),          intent(in)    :: message
    integer(I4), optional,     intent(in)    :: b, i, j, k, face, count
    real(R8),    optional,     intent(in)    :: stroke
    ! Local
    type(morph_event_t), allocatable :: tmp(:)

    if ( .not. allocated(self%ev) ) allocate( self%ev(8) )
    if ( self%n == size(self%ev) ) then
      allocate( tmp(2*size(self%ev)) )
      tmp(1:self%n) = self%ev(1:self%n)
      call move_alloc ( tmp, self%ev )
    endif
    self%n = self%n + 1
    associate ( e => self%ev(self%n) )
      e%code    = code
      e%stops   = ( code < MORPH_EV_STOP_LIMIT )
      e%message = message
      e%b = 0; e%i = 0; e%j = 0; e%k = 0; e%face = 0; e%count = 1; e%stroke = 0.0_R8
      if ( present(b) )      e%b      = b
      if ( present(i) )      e%i      = i
      if ( present(j) )      e%j      = j
      if ( present(k) )      e%k      = k
      if ( present(face) )   e%face   = face
      if ( present(count) )  e%count  = count
      if ( present(stroke) ) e%stroke = stroke
    end associate
  end subroutine events_add


  pure subroutine events_clear ( self )
    class(morph_event_list_t), intent(inout) :: self
    self%n = 0
  end subroutine events_clear


  pure logical function events_any_stop ( self ) result ( res )
    class(morph_event_list_t), intent(in) :: self
    integer(I4) :: e
    res = .false.
    do e = 1, self%n
      if ( self%ev(e)%stops ) res = .true.
    enddo
  end function events_any_stop

end module Morph_Events
