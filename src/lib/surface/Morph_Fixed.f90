!>@brief Fixed nodes: index ranges declared by the caller, plus the frozen
!> blocks, closed over the connections: a node fixed in one block is fixed in
!> every block that shares it. The closure uses the node classes of
!> Morph_Node_Classes, so one pass gives the fixed point.
!>
!> MORPH must never `use` a FUSS module -- see Morph_Types_m.
module Morph_Fixed
  use iso_fortran_env, only: I4 => int32
  use Morph_Types_m
  use Morph_Surface_Types
  use Morph_Domain_m,  only: morph_domain_t
  use Morph_Surface,   only: Morph_Node_Classes, Morph_Node_Offsets

  implicit none
  private

  public :: Morph_Fixed_Masks

contains

  !> intervals(7, n): (b, i0, i1, j0, j1, k0, k1), node indices from 0, -1 = the
  !> whole axis. Fills domain%fixed(b)%m; MORPH_ERR_INPUT on an interval out of
  !> range.
  subroutine Morph_Fixed_Masks ( domain, intervals, status )
    type(morph_domain_t), intent(inout) :: domain
    integer(I4),          intent(in)    :: intervals(:,:)
    type(morph_status_t), intent(out)   :: status
    integer(I4), allocatable :: parent(:), off(:)
    logical,     allocatable :: froot(:)
    integer(I4) :: nb, b, t, a, lo(3), hi(3), i, j, k, id, r, ntot

    call status%clear()
    nb = domain%nblocks
    if ( size(intervals,1) /= 7 .and. size(intervals,2) > 0 ) then
      call status%fail ( MORPH_ERR_SHAPE, 'Morph_Fixed_Masks: intervals must be (7, n)' )
      return
    endif
    if ( allocated(domain%fixed) ) then
      do b = 1, size(domain%fixed)
        if ( allocated(domain%fixed(b)%m) ) deallocate ( domain%fixed(b)%m )
      enddo
      deallocate ( domain%fixed )
    endif
    allocate ( domain%fixed(nb) )
    do b = 1, nb
      allocate ( domain%fixed(b)%m(0:domain%blk(b)%dim(1), 0:domain%blk(b)%dim(2), 0:domain%blk(b)%dim(3)) )
      domain%fixed(b)%m = .false.
      if ( allocated(domain%frozen) ) then
        if ( domain%frozen(b) ) domain%fixed(b)%m = .true.
      endif
    enddo
    do t = 1, size(intervals,2)
      b = intervals(1,t)
      if ( b < 1 .or. b > nb ) then
        call status%fail ( MORPH_ERR_INPUT, 'Morph_Fixed_Masks: interval on a block out of range', b = b )
        return
      endif
      do a = 1, 3
        lo(a) = intervals(2*a,t) ; hi(a) = intervals(2*a+1,t)
        if ( lo(a) == -1 ) lo(a) = 0
        if ( hi(a) == -1 ) hi(a) = domain%blk(b)%dim(a)
        if ( lo(a) < 0 .or. hi(a) > domain%blk(b)%dim(a) .or. lo(a) > hi(a) ) then
          call status%fail ( MORPH_ERR_INPUT, 'Morph_Fixed_Masks: interval out of the block', b = b )
          return
        endif
      enddo
      domain%fixed(b)%m(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) = .true.
    enddo
    ! closure over the node classes: one pass, since a class is complete
    call Morph_Node_Classes ( domain, parent, status )
    if ( .not. status%ok() ) return
    call Morph_Node_Offsets ( domain, off, ntot )
    allocate ( froot(ntot) ) ; froot = .false.
    do b = 1, nb
      do k = 0, domain%blk(b)%dim(3) ; do j = 0, domain%blk(b)%dim(2) ; do i = 0, domain%blk(b)%dim(1)
        if ( .not. domain%fixed(b)%m(i,j,k) ) cycle
        id = off(b) + 1 + i + (domain%blk(b)%dim(1)+1) * ( j + (domain%blk(b)%dim(2)+1) * k )
        call root_of ( id, r )
        froot(r) = .true.
      enddo ; enddo ; enddo
    enddo
    do b = 1, nb
      do k = 0, domain%blk(b)%dim(3) ; do j = 0, domain%blk(b)%dim(2) ; do i = 0, domain%blk(b)%dim(1)
        id = off(b) + 1 + i + (domain%blk(b)%dim(1)+1) * ( j + (domain%blk(b)%dim(2)+1) * k )
        call root_of ( id, r )
        if ( froot(r) ) domain%fixed(b)%m(i,j,k) = .true.
      enddo ; enddo ; enddo
    enddo

  contains

    subroutine root_of ( x, rt )
      integer(I4), intent(in)  :: x
      integer(I4), intent(out) :: rt
      rt = x
      do while ( parent(rt) /= rt )
        parent(rt) = parent(parent(rt))
        rt = parent(rt)
      enddo
    end subroutine root_of

  end subroutine Morph_Fixed_Masks

end module Morph_Fixed
