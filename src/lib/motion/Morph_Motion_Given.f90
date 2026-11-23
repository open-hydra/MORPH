!>@brief The `given` law: the new node positions of a block are decided
!> elsewhere (the regression: the level set on the wall, the multi-block volume
!> law inside) and this law applies them, so metrics, swept volumes, the GCL
!> and the quality checks go through the same Morph_Update path as for every
!> other law.
!>
!> One step: the decider hands the positions over with `fill`; the caller then
!> calls Morph_Update (directly, or through domain%advance), whose `apply` uses
!> them once. The `pending` flag pins that order:
!>   * `fill` refuses (MORPH_ERR_STATE) while the previous positions were never
!>     applied: the mesh would lag behind whoever decided them;
!>   * `apply` refuses (MORPH_ERR_STATE) when nothing is pending: a second
!>     update with the same positions would move nothing and still count the
!>     step's time.
!> Morph_Update rolls node_old and vol_old before it calls `apply`, so a
!> refused `apply` has already rolled the block's time levels, as any failing
!> law has.
!>
!> Copy rule: a node not marked in `set` takes node_old bit for bit; adding a
!> zero displacement would turn -0.0 into +0.0. With nothing set, the block is
!> bit-identical to the previous step and its swept volumes are exactly zero.
!>
!> `pos`, `set` and `pending` are public so that a solver keeping its own node
!> arrays can read the target positions; write them only through `fill`.
!>
!> MORPH must never `use` a FUSS module -- see Morph_Types_m.
module Morph_Motion_Given
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m
  use Morph_Motion_m

  implicit none
  private

  public :: morph_motion_given_t

  type, extends(morph_motion_t) :: morph_motion_given_t
    type(morph_vec3_t), allocatable :: pos(:,:,:)   !> (0:dim(1), 0:dim(2), 0:dim(3)): target positions
    logical,            allocatable :: set(:,:,:)   !> .true. where pos holds a new position
    logical :: pending = .false.                    !> filled and not applied yet
  contains
    procedure :: apply => given_apply
    procedure :: fill  => given_fill
  end type morph_motion_given_t

contains

  !> Hand over the target positions of one step. The arrays are moved, not
  !> copied: on success `pos` and `set` come back deallocated. When the call
  !> is refused, the law and the arrays are left as they were.
  subroutine given_fill ( self, pos, set, status )
    class(morph_motion_given_t),     intent(inout) :: self
    type(morph_vec3_t), allocatable, intent(inout) :: pos(:,:,:)
    logical,            allocatable, intent(inout) :: set(:,:,:)
    type(morph_status_t),            intent(out)   :: status

    call status%clear()

    if ( self%pending ) then
      call status%fail ( MORPH_ERR_STATE, 'given: the positions of the previous fill were never applied' )
      return
    endif
    if ( .not. allocated(pos) .or. .not. allocated(set) ) then
      call status%fail ( MORPH_ERR_SHAPE, 'given: fill without target positions' )
      return
    endif
    if ( any( lbound(pos) /= lbound(set) ) .or. any( ubound(pos) /= ubound(set) ) ) then
      call status%fail ( MORPH_ERR_SHAPE, 'given: target positions and mask of different shape' )
      return
    endif

    call move_alloc ( pos, self%pos )
    call move_alloc ( set, self%set )
    self%pending = .true.

  end subroutine given_fill


  subroutine given_apply ( self, node, node_old, dim, t, dt, status )
    class(morph_motion_given_t), intent(inout) :: self
    type(morph_vec3_t),          intent(inout) :: node(0:,0:,0:)
    type(morph_vec3_t),          intent(in)    :: node_old(0:,0:,0:)
    integer(I4),                 intent(in)    :: dim(3)
    real(R8),                    intent(in)    :: t, dt
    type(morph_status_t),        intent(out)   :: status
    ! Local
    integer(I4) :: i, j, k

    call status%clear()

    if ( .not. self%pending ) then
      call status%fail ( MORPH_ERR_STATE, 'given: no positions waiting (fill before each update)' )
      return
    endif
    if ( .not. allocated(self%pos) .or. .not. allocated(self%set) ) then
      call status%fail ( MORPH_ERR_SHAPE, 'given: no target positions' )
      return
    endif
    if ( any( ubound(self%pos) /= dim ) .or. any( lbound(self%pos) /= 0 ) .or. &
         any( ubound(self%set) /= dim ) .or. any( lbound(self%set) /= 0 ) ) then
      call status%fail ( MORPH_ERR_SHAPE, 'given: the target arrays do not match the block' )
      return
    endif

    !$omp parallel do collapse(3) private(i,j,k)
    do k = 0, dim(3)
    do j = 0, dim(2)
    do i = 0, dim(1)
      if ( self%set(i,j,k) ) then
        node(i,j,k) = self%pos(i,j,k)
      else
        node(i,j,k) = node_old(i,j,k)
      endif
    enddo; enddo; enddo
    !$omp end parallel do

    self%pending = .false.

  end subroutine given_apply

end module Morph_Motion_Given
