!>@brief The multi-block mesh state MORPH owns when it runs on its own (plan 09
!>       section 3: the state type deferred to the extraction).
!>
!> A solver that embeds MORPH keeps its own node arrays and calls Morph_API one
!> block at a time, as FUSS's adapter does. This type is for the standalone
!> driver, and for any caller that wants MORPH to hold the mesh: per block it
!> owns the nodes at the current and at the previous time level, the geometry
!> and the motion law. The previous positions (FUSS's blk%node_old, an open
!> ownership question in plan 09 section 11) are therefore MORPH's here.
!>
!> INTERIOR GEOMETRY ONLY. Morph_Init and Morph_Update fill the cells 1..dim of
!> each block; the ghost layer of geom (index 0 and dim+1) is never written.
!> Ghost-cell metrics need the block connectivity, which belongs to the solver
!> (user decision 2026-10-24: they stay in FUSS).
!>
!> Like the rest of MORPH it never stops: every procedure returns a
!> morph_status_t and, for a failure inside a block, the index of that block,
!> which morph_status_t itself does not carry.
!>
!> MORPH must never `use` a FUSS module -- see Morph_Types_m.
module Morph_Domain_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_API

  implicit none
  private

  public :: morph_block_t, morph_domain_t

  !> One structured block of hexahedra.
  type :: morph_block_t
    integer(I4)                        :: dim(3) = 0       !> cells per direction
    type(morph_vec3_t),    allocatable :: node(:,:,:)      !> (0:dim(1), 0:dim(2), 0:dim(3)), positions at domain%t
    type(morph_vec3_t),    allocatable :: node_old(:,:,:)  !> positions at the start of the last step
    type(morph_geom_t)                 :: geom             !> geometry of the cells 1..dim
    !> The block's own law: the prescribed law caches the reference mesh of the
    !> block it is first applied to, so one instance cannot serve two blocks.
    class(morph_motion_t), allocatable :: law
  end type morph_block_t

  type :: morph_domain_t
    integer(I4)                      :: nblocks = 0
    type(morph_block_t), allocatable :: blk(:)
    real(R8)                         :: t    = 0.0_R8     !> time of the current node positions
    integer(I4)                      :: step = 0          !> steps taken since init
  contains
    procedure :: allocate_blocks => domain_allocate_blocks
    procedure :: init            => domain_init
    procedure :: advance         => domain_advance
    procedure :: gcl_residual    => domain_gcl_residual
    procedure :: outflow_ratio   => domain_outflow_ratio
    procedure :: free            => domain_free
  end type morph_domain_t

contains

  !> Allocate one block per column of dims (cells per direction), every node at
  !> the origin and the static law on every block. The caller then fills the
  !> nodes and replaces the laws it wants to move.
  subroutine domain_allocate_blocks ( self, dims, status )
    class(morph_domain_t), intent(inout) :: self
    integer(I4),           intent(in)    :: dims(:,:)      !> (3, nblocks)
    type(morph_status_t),  intent(out)   :: status
    ! Local
    integer(I4) :: b

    call status%clear()
    call self%free()

    if ( size(dims,1) /= 3 .or. size(dims,2) < 1 ) then
      call status%fail ( MORPH_ERR_SHAPE, 'domain: dims must be (3, nblocks) with at least one block' )
      return
    endif
    if ( any( dims < 1 ) ) then
      call status%fail ( MORPH_ERR_SHAPE, 'domain: every block needs at least one cell in each direction' )
      return
    endif

    self%nblocks = size(dims,2)
    allocate( self%blk(self%nblocks) )
    do b = 1, self%nblocks
      self%blk(b)%dim = dims(:,b)
      allocate( self%blk(b)%node    ( 0:dims(1,b), 0:dims(2,b), 0:dims(3,b) ) )
      allocate( self%blk(b)%node_old( 0:dims(1,b), 0:dims(2,b), 0:dims(3,b) ) )
      allocate( morph_motion_static_t :: self%blk(b)%law )
      self%blk(b)%law%name = 'static'
    enddo

  end subroutine domain_allocate_blocks


  !> Geometry of the mesh as it stands, Morph_Init on every block, at time t0.
  !> The GCL residual taken straight after is exactly zero.
  subroutine domain_init ( self, t0, status, ib )
    class(morph_domain_t), intent(inout) :: self
    real(R8),              intent(in)    :: t0
    type(morph_status_t),  intent(out)   :: status
    integer(I4),           intent(out)   :: ib            !> the failing block, 0 if none
    ! Local
    integer(I4) :: b

    call status%clear()
    ib = 0

    if ( self%nblocks < 1 .or. .not. allocated(self%blk) ) then
      call status%fail ( MORPH_ERR_SHAPE, 'domain: init called before allocate_blocks' )
      return
    endif

    self%t    = t0
    self%step = 0
    do b = 1, self%nblocks
      self%blk(b)%node_old = self%blk(b)%node
      call Morph_Init ( self%blk(b)%node, self%blk(b)%dim, self%blk(b)%geom, status )
      if ( .not. status%ok() ) then
        ib = b
        return
      endif
    enddo

  end subroutine domain_init


  !> Advance every block by dt under its own law, Morph_Update one block at a
  !> time (each call is OpenMP-parallel inside, so this must not be called from
  !> a parallel region). On a failure in block ib the blocks before it have
  !> moved and the others have not: the domain is left part-way through the
  !> step, and t and step are not advanced. The caller should stop.
  subroutine domain_advance ( self, dt, status, ib )
    class(morph_domain_t), intent(inout) :: self
    real(R8),              intent(in)    :: dt
    type(morph_status_t),  intent(out)   :: status
    integer(I4),           intent(out)   :: ib            !> the failing block, 0 if none
    ! Local
    integer(I4) :: b

    call status%clear()
    ib = 0

    do b = 1, self%nblocks
      call Morph_Update ( self%blk(b)%node, self%blk(b)%node_old, self%blk(b)%dim, &
                          self%t, dt, self%blk(b)%law, self%blk(b)%geom, status )
      if ( .not. status%ok() ) then
        ib = b
        return
      endif
    enddo

    self%t    = self%t + dt
    self%step = self%step + 1

  end subroutine domain_advance


  !> The largest discrete-GCL residual of the last step over all blocks,
  !> absolute [m^3] and relative to the cell volume (Morph_GCL_Residual, block
  !> by block), and the block where each maximum is found (0 when it is zero).
  subroutine domain_gcl_residual ( self, r_abs, r_rel, ib_abs, ib_rel )
    class(morph_domain_t), intent(in)  :: self
    real(R8),              intent(out) :: r_abs, r_rel
    integer(I4),           intent(out) :: ib_abs, ib_rel
    ! Local
    integer(I4) :: b
    real(R8)    :: a, r

    r_abs = 0.0_R8 ; r_rel = 0.0_R8
    ib_abs = 0     ; ib_rel = 0
    do b = 1, self%nblocks
      call Morph_GCL_Residual ( self%blk(b)%geom, self%blk(b)%dim, a, r )
      if ( a > r_abs ) then
        r_abs = a ; ib_abs = b
      endif
      if ( r > r_rel ) then
        r_rel = r ; ib_rel = b
      endif
    enddo

  end subroutine domain_gcl_residual


  !> Largest fraction of its own volume that any cell swept out in the last
  !> step (Morph_Swept_Outflow_Ratio, block by block), and where.
  subroutine domain_outflow_ratio ( self, ratio, ib, ci, cj, ck )
    class(morph_domain_t), intent(in)  :: self
    real(R8),              intent(out) :: ratio
    integer(I4),           intent(out) :: ib, ci, cj, ck
    ! Local
    integer(I4) :: b, i, j, k
    real(R8)    :: r

    ratio = 0.0_R8
    ib = 0 ; ci = 0 ; cj = 0 ; ck = 0
    do b = 1, self%nblocks
      call Morph_Swept_Outflow_Ratio ( self%blk(b)%geom, self%blk(b)%dim, r, i, j, k )
      if ( r > ratio ) then
        ratio = r
        ib = b ; ci = i ; cj = j ; ck = k
      endif
    enddo

  end subroutine domain_outflow_ratio


  !> Release every block (nodes, geometry and law with it).
  subroutine domain_free ( self )
    class(morph_domain_t), intent(inout) :: self

    if ( allocated(self%blk) ) deallocate( self%blk )
    self%nblocks = 0
    self%t       = 0.0_R8
    self%step    = 0

  end subroutine domain_free

end module Morph_Domain_m
