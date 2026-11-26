!>@brief The level-set field on a 2-D grid of blocks split among ranks. The
!> grid has nblk(1) x nblk(2) blocks of nbn x nbn nodes, spacing h; block
!> b = ib + (jb-1) nblk(1) holds the nodes i = (ib-1) nbn + 1 .. ib nbn,
!> j = (jb-1) nbn + 1 .. jb nbn of the whole grid and lives on rank owner(b).
!> A rank keeps its own active blocks in block order, its slots:
!> phi(nbn, nbn, slot).
!>
!> THE BAND. Only the active blocks hold phi and move. Which blocks are
!> active is the same on every rank: all of them, or the set given to
!> Morph_LS_Grid_Init, then the set of each Morph_LS_Grid_Update, where a
!> block is active when at least one rank asks for it (the rank whose block
!> the front approaches asks for the blocks beyond it, whoever owns them).
!> The ghosts toward an inactive neighbour come from the block's own nodes
!> (THE GHOSTS), and no message carries them. Every active node moves at its
!> full speed. There is no window that
!> slows the nodes far from the front: a node whose speed is zero never
!> changes, so it never comes back into the window, and the levels ahead of a
!> moving front pile up against it. Without one, the field near the front is
!> the field of the whole grid up to two errors that come from far from it:
!> the extrapolated ghosts at the edge of the band and the values a block
!> starts from when it becomes active. Both fall as the margin between the
!> front and the inactive blocks grows. The caller chooses that margin, asks
!> for the blocks and fills the new ones: Morph_LS_Grid_Remap moves each of
!> its fields to the new slots and leaves a NaN in the slots of the blocks
!> just activated, so that a slot left unfilled shows in the results.
!>
!> THE GHOSTS. A step works on w(1-NG:nbn+NG, 1-NG:nbn+NG, slot), the layout
!> of the kernels of Morph_LS_HJ for each block. Morph_LS_Halo_Fill fills the
!> NG layers along each side with the nodes of the active block beyond it (a
!> block has at least NG nodes per side, so one neighbour has them all). At
!> the edge of the grid it extrapolates linearly (Morph_LS_Ghosts_Extrap), as
!> on the whole grid. Toward an inactive block, given the speed V (the step
!> gives it), Morph_LS_Ghosts_Upwind chooses line by line: linear where the
!> characteristics enter the block, quadratic where they leave it; without V,
!> linear. Where they leave, the quadratic follows the curvature of the level
!> sets and brings the band's front closer to the whole grid's; where they
!> enter, what the ghosts carry moves into the band, and with a quadratic
!> there a band drifted away from the whole grid (the README's narrow band
!> gives the numbers). The corners are neither filled nor read.
!> With every block active, no side faces an inactive block and every value
!> a stencil reads is the one it reads on the whole grid, so
!> Morph_LS_HJ_Step gives the bits of Morph_LS_HJ_Step_Full for any blocks,
!> any split among ranks and any number of threads.
!>
!> THE MESSAGES. A strip is the NG x nbn ghost values along one side of one
!> block, in the order of that block's array (i fastest). Each rank lists the
!> strips it sends and receives from the grid and its active set alone, in
!> one order: by rank, then by the block that receives, then by its side
!> (W, E, S, N). Only strips between two active blocks travel.
!>
!> THE OTHER RANKS. Morph_LS_HJ_Step, Morph_LS_Halo_Fill and
!> Morph_LS_Grid_Update are collectives: every rank calls them, in the same
!> order, outside any OpenMP region. The step and the update decide with all
!> the ranks whether they refuse (the arguments of each), and the step
!> whether it moves anything (V = 0 on every block), so either every rank
!> leaves early or none does. Halo_Fill on its own refuses a w of another
!> shape on that rank only, before the exchange, and the other ranks then
!> wait in it: the step checks the shapes first, together. Init, Remap,
!> Scatter and Gather talk to no other rank.
!>
!> Never stops, no I/O: refused arguments return MORPH_ERR_SHAPE,
!> MORPH_ERR_INPUT or MORPH_ERR_EXCHANGE in the status.
module Morph_LS_Grid
  use iso_fortran_env, only: I4 => int32, I8 => int64, R8 => real64
  use Morph_Types_m
  use Morph_Exchange,  only: morph_exchange_t
  use Morph_LS_HJ,     only: MORPH_LS_NG, Morph_LS_Ghosts_Extrap, Morph_LS_Ghosts_Upwind, Morph_LS_Rate, &
                             Morph_LS_RK3_Combine

  implicit none
  private

  integer(I4), parameter :: NG = MORPH_LS_NG
  !> A quiet NaN from its IEEE 754 bits, 0x7FF8000000000000 (the library uses no
  !> module beyond iso_fortran_env, omp_lib and its own).
  real(R8),    parameter :: QNAN = transfer( 9221120237041090560_I8, 0.0_R8 )

  type, public :: morph_ls_grid_t
    integer(I4) :: nbn = 0                   !> nodes along each side of a block
    integer(I4) :: nblk(2) = 0               !> blocks along each axis
    real(R8)    :: h = 0.0_R8                !> node spacing
    integer(I4) :: me = 0                    !> this rank
    integer(I4) :: nranks = 0
    integer(I4), allocatable :: owner(:)     !> the rank of each block
    logical,     allocatable :: active(:)    !> each block active or not, the same on every rank
    integer(I4), allocatable :: blk(:)       !> slot -> block: this rank's active blocks in block order
    integer(I4), allocatable :: slot(:)      !> block -> its slot here, 0 if inactive or on another rank
    integer(I4), allocatable :: nbr(:,:)     !> (4, slot): the active block beyond side W, E, S, N; 0 at an edge
    integer(I4), allocatable :: rpos(:,:)    !> (4, slot): the strip's place among those received, 0 if not received
    integer(I4), allocatable :: send(:,:)    !> (2, strip): the slot that gives it, the side of the block that receives
    integer(I4), allocatable :: scount(:)    !> (0:nranks-1): values sent to each rank
    integer(I4), allocatable :: rcount(:)    !> (0:nranks-1): values received from each rank
  end type morph_ls_grid_t

  public :: Morph_LS_Grid_Init, Morph_LS_Grid_Update, Morph_LS_Grid_Remap
  public :: Morph_LS_Grid_Scatter, Morph_LS_Grid_Gather
  public :: Morph_LS_Halo_Fill, Morph_LS_HJ_Step

contains

  !> The blocks beyond the sides W, E, S, N of block b; 0 at the edge.
  pure function neighbours ( nblk, b ) result ( nb )
    integer(I4), intent(in) :: nblk(2), b
    integer(I4) :: nb(4), ib, jb
    ib = mod(b-1, nblk(1)) + 1
    jb = (b-1)/nblk(1) + 1
    nb = 0
    if ( ib > 1 )       nb(1) = b - 1
    if ( ib < nblk(1) ) nb(2) = b + 1
    if ( jb > 1 )       nb(3) = b - nblk(1)
    if ( jb < nblk(2) ) nb(4) = b + nblk(1)
  end function neighbours


  !> The ghosts of side d of a block: (i1, i2, j1, j2).
  pure function ghost_box ( nbn, d ) result ( g )
    integer(I4), intent(in) :: nbn, d
    integer(I4) :: g(4)
    select case ( d )
    case ( 1 )
      g = [ 1-NG,  0,      1,     nbn    ]
    case ( 2 )
      g = [ nbn+1, nbn+NG, 1,     nbn    ]
    case ( 3 )
      g = [ 1,     nbn,    1-NG,  0      ]
    case default
      g = [ 1,     nbn,    nbn+1, nbn+NG ]
    end select
  end function ghost_box


  !> The nodes of the block beyond side d that fill the ghosts of side d: the
  !> ghost box moved by one block along the axis.
  pure function source_box ( nbn, d ) result ( g )
    integer(I4), intent(in) :: nbn, d
    integer(I4) :: g(4)
    g = ghost_box ( nbn, d )
    select case ( d )
    case ( 1 )
      g(1:2) = g(1:2) + nbn
    case ( 2 )
      g(1:2) = g(1:2) - nbn
    case ( 3 )
      g(3:4) = g(3:4) + nbn
    case default
      g(3:4) = g(3:4) - nbn
    end select
  end function source_box


  !> The number of this rank's active blocks; -1 for a grid not initialised,
  !> so that no shape matches it.
  pure integer(I4) function nslots ( grid )
    type(morph_ls_grid_t), intent(in) :: grid
    nslots = -1
    if ( allocated(grid%blk) ) nslots = size(grid%blk)
  end function nslots


  !> This rank's slots and the strips it sends and receives, from the owners
  !> and the active set of the grid.
  subroutine build_schedule ( grid )
    type(morph_ls_grid_t), intent(inout) :: grid
    integer(I4), allocatable :: tmp(:,:)
    integer(I4) :: nb, nr, me, b, s, d, p, k, nmine, nbs(4)

    nb = size(grid%owner)
    nr = grid%nranks
    me = grid%me
    grid%blk = pack( [ (b, b = 1, nb) ], grid%owner == me .and. grid%active )
    nmine = size(grid%blk)
    if ( allocated(grid%slot) )   deallocate( grid%slot )
    if ( allocated(grid%nbr) )    deallocate( grid%nbr )
    if ( allocated(grid%rpos) )   deallocate( grid%rpos )
    if ( allocated(grid%scount) ) deallocate( grid%scount )
    if ( allocated(grid%rcount) ) deallocate( grid%rcount )
    allocate( grid%slot(nb), grid%nbr(4,nmine), grid%rpos(4,nmine) )
    allocate( grid%scount(0:nr-1), grid%rcount(0:nr-1) )
    grid%slot   = 0
    grid%rpos   = 0
    grid%scount = 0
    grid%rcount = 0
    do s = 1, nmine
      grid%slot(grid%blk(s)) = s
      nbs = neighbours ( grid%nblk, grid%blk(s) )
      do d = 1, 4
        if ( nbs(d) /= 0 ) then
          if ( .not. grid%active(nbs(d)) ) nbs(d) = 0      ! the edge of the band
        endif
      enddo
      grid%nbr(:,s) = nbs
    enddo

    ! received: from each other rank in rank order, my blocks in order, their sides
    k = 0
    do p = 0, nr-1
      if ( p == me ) cycle
      do s = 1, nmine
        do d = 1, 4
          if ( grid%nbr(d,s) == 0 ) cycle
          if ( grid%owner(grid%nbr(d,s)) /= p ) cycle
          k = k + 1
          grid%rpos(d,s) = k
          grid%rcount(p) = grid%rcount(p) + NG*grid%nbn
        enddo
      enddo
    enddo

    ! sent: to each other rank in rank order, its active blocks in order, their
    ! sides beyond which the block is mine and active (one per side of my blocks)
    allocate( tmp(2, 4*nmine) )
    k = 0
    do p = 0, nr-1
      if ( p == me ) cycle
      do b = 1, nb
        if ( grid%owner(b) /= p .or. .not. grid%active(b) ) cycle
        nbs = neighbours ( grid%nblk, b )
        do d = 1, 4
          if ( nbs(d) == 0 ) cycle
          if ( grid%owner(nbs(d)) /= me .or. .not. grid%active(nbs(d)) ) cycle
          k = k + 1
          tmp(:,k) = [ grid%slot(nbs(d)), d ]
          grid%scount(p) = grid%scount(p) + NG*grid%nbn
        enddo
      enddo
    enddo
    grid%send = tmp(:,1:k)
  end subroutine build_schedule


  !> The grid, its split among ranks, its active set (active, one value per
  !> block; all the blocks when absent) and the list of the strips this rank
  !> sends and receives. Talks to no other rank: exch gives this rank and
  !> their number, and active must be the same on every rank. Refused, grid
  !> not initialised: nbn < NG or fewer than one block along an axis, owner or
  !> active without one value per block (MORPH_ERR_SHAPE); h <= 0 or a block
  !> on a rank that does not exist (MORPH_ERR_INPUT); an exchange that is not
  !> rank 0..n-1 of n >= 1 (MORPH_ERR_EXCHANGE).
  subroutine Morph_LS_Grid_Init ( grid, nbn, nblk, h, owner, exch, status, active )
    type(morph_ls_grid_t),   intent(out) :: grid
    integer(I4),             intent(in)  :: nbn, nblk(2)
    real(R8),                intent(in)  :: h
    integer(I4),             intent(in)  :: owner(:)
    class(morph_exchange_t), intent(in)  :: exch
    type(morph_status_t),    intent(out) :: status
    logical, optional,       intent(in)  :: active(:)
    integer(I4) :: nb, me, nr
    character(len=MORPH_MSGLEN) :: msg

    call status%clear()
    if ( nbn < NG .or. nblk(1) < 1 .or. nblk(2) < 1 ) then
      write(msg,'(A,I0,A,I0,A,I0,A,I0,A)') 'level-set grid: blocks of ', nbn, ' x ', nbn, ' nodes, ', &
            nblk(1), ' x ', nblk(2), ' blocks (at least 3 nodes and one block per axis)'
      call status%fail ( MORPH_ERR_SHAPE, msg )
      return
    endif
    nb = nblk(1)*nblk(2)
    if ( size(owner) /= nb ) then
      write(msg,'(A,I0,A,I0,A)') 'level-set grid: ', size(owner), ' owners for ', nb, ' blocks'
      call status%fail ( MORPH_ERR_SHAPE, msg )
      return
    endif
    if ( present(active) ) then
      if ( size(active) /= nb ) then
        write(msg,'(A,I0,A,I0,A)') 'level-set grid: ', size(active), ' active flags for ', nb, ' blocks'
        call status%fail ( MORPH_ERR_SHAPE, msg )
        return
      endif
    endif
    if ( .not. (h > 0.0_R8) ) then
      write(msg,'(A,ES10.3,A)') 'level-set grid: h = ', h, ' (h > 0)'
      call status%fail ( MORPH_ERR_INPUT, msg )
      return
    endif
    me = exch%rank()
    nr = exch%nranks()
    if ( me < 0 .or. me >= nr ) then
      write(msg,'(A,I0,A,I0,A)') 'level-set grid: the exchange is rank ', me, ' of ', nr, ' ranks'
      call status%fail ( MORPH_ERR_EXCHANGE, msg )
      return
    endif
    if ( any( owner < 0 .or. owner >= nr ) ) then
      write(msg,'(A,I0,A,I0,A,I0)') 'level-set grid: a block on rank ', minval(owner), ' .. ', &
            maxval(owner), ', ranks 0 .. ', nr-1
      call status%fail ( MORPH_ERR_INPUT, msg )
      return
    endif

    grid%nbn    = nbn
    grid%nblk   = nblk
    grid%h      = h
    grid%me     = me
    grid%nranks = nr
    grid%owner  = owner
    allocate( grid%active(nb) )
    grid%active = .true.
    if ( present(active) ) grid%active = active
    call build_schedule ( grid )
  end subroutine Morph_LS_Grid_Init


  !> A new active set, decided with all the ranks: block b is active when at
  !> least one rank asks for it (want(b) true there), inactive when none does.
  !> from(s) is the slot that the block of the new slot s had before, 0 for a
  !> block that was inactive: Morph_LS_Grid_Remap moves any field of slots
  !> with it. A collective (see the header). First every rank checks its
  !> arguments and all share the verdicts: refused on every rank, the grid
  !> untouched and from(s) = s, when one rank refuses (its own message on that
  !> rank, the first refusing rank's code elsewhere): want without one value
  !> per block or a grid not initialised (MORPH_ERR_SHAPE), an exchange that
  !> is not the grid's (MORPH_ERR_EXCHANGE). Grids of different sizes on
  !> different ranks are refused on all of them (MORPH_ERR_EXCHANGE).
  subroutine Morph_LS_Grid_Update ( grid, exch, want, from, status )
    type(morph_ls_grid_t),    intent(inout) :: grid
    class(morph_exchange_t),  intent(inout) :: exch
    logical,                  intent(in)    :: want(:)
    integer(I4), allocatable, intent(out)   :: from(:)
    type(morph_status_t),     intent(out)   :: status
    real(R8),    allocatable :: got(:), asked(:)
    integer(I4), allocatable :: counts(:), before(:)
    logical,     allocatable :: next(:)
    integer(I4) :: nb, nr, code, p, b, s, k, nask
    character(len=MORPH_MSGLEN) :: msg

    call status%clear()
    from = [ (s, s = 1, nslots(grid)) ]
    nb   = -1
    if ( nslots(grid) >= 0 ) nb = size(grid%owner)
    code = MORPH_OK
    msg  = ''
    if ( nslots(grid) < 0 ) then
      code = MORPH_ERR_SHAPE
      msg  = 'level-set band update: the grid is not initialised'
    else if ( size(want) /= nb ) then
      code = MORPH_ERR_SHAPE
      write(msg,'(A,I0,A,I0,A)') 'level-set band update: ', size(want), ' requests for ', nb, ' blocks'
    else if ( exch%rank() /= grid%me .or. exch%nranks() /= grid%nranks ) then
      code = MORPH_ERR_EXCHANGE
      write(msg,'(A,I0,A,I0,A,I0,A,I0,A)') 'level-set band update: the exchange is rank ', exch%rank(), &
            ' of ', exch%nranks(), ', the grid rank ', grid%me, ' of ', grid%nranks
    endif
    nask = count(want)

    ! the verdicts, the sizes of the grids and the number of requests of every rank
    nr = exch%nranks()
    allocate( counts(0:nr-1), got(3*nr) )
    counts = 3
    call exch%allgather ( [ real(code,R8), real(nb,R8), real(nask,R8) ], got, counts, status )
    if ( .not. status%ok() ) return
    if ( code /= MORPH_OK ) then
      call status%fail ( code, msg )
      return
    endif
    do p = 0, nr-1
      if ( got(3*p+1) /= 0.0_R8 ) then
        write(msg,'(A,I0,A,I0,A)') 'level-set band update: refused on rank ', p, ' (status ', nint(got(3*p+1)), ')'
        call status%fail ( nint(got(3*p+1), I4), msg )
        return
      endif
    enddo
    if ( any( got(2::3) /= real(nb,R8) ) ) then
      write(msg,'(A,I0,A)') 'level-set band update: the ranks hold grids of different sizes (', nb, ' blocks here)'
      call status%fail ( MORPH_ERR_EXCHANGE, msg )
      return
    endif

    ! the requests of every rank, then the new set: what any rank asks for
    do p = 0, nr-1
      counts(p) = nint(got(3*p+3), I4)
    enddo
    allocate( asked(sum(counts)) )
    call exch%allgather ( real( pack( [ (b, b = 1, nb) ], want ), R8 ), asked, counts, status )
    if ( .not. status%ok() ) return
    allocate( next(nb) )
    next = .false.
    do k = 1, size(asked)
      next(nint(asked(k), I4)) = .true.
    enddo

    before = grid%slot
    grid%active = next
    call build_schedule ( grid )
    deallocate( from )
    allocate( from(size(grid%blk)) )
    do s = 1, size(grid%blk)
      from(s) = before(grid%blk(s))
    enddo
  end subroutine Morph_LS_Grid_Update


  !> A field of this rank's slots, a(:, :, slot) (phi, V or any other), moved
  !> to the slots of an update: the new slot s takes the old slot from(s) bit
  !> for bit, or a quiet NaN where from(s) = 0 (a block just activated, which
  !> the caller fills). Talks to no other rank. Refused, a untouched: a not
  !> allocated, or from naming a slot that a does not have (MORPH_ERR_SHAPE).
  subroutine Morph_LS_Grid_Remap ( from, a, status )
    integer(I4),           intent(in)    :: from(:)
    real(R8), allocatable, intent(inout) :: a(:,:,:)
    type(morph_status_t),  intent(out)   :: status
    real(R8), allocatable :: b(:,:,:)
    integer(I4) :: s
    character(len=MORPH_MSGLEN) :: msg

    call status%clear()
    if ( .not. allocated(a) ) then
      call status%fail ( MORPH_ERR_SHAPE, 'level-set remap: the field is not allocated' )
      return
    endif
    if ( minval(from) < 0 .or. maxval(from) > size(a,3) ) then     ! never on no slot: +huge, -huge
      write(msg,'(A,I0,A,I0,A,I0,A)') 'level-set remap: slots ', minval(from), ' .. ', maxval(from), &
            ' of a field of ', size(a,3), ' slots'
      call status%fail ( MORPH_ERR_SHAPE, msg )
      return
    endif
    allocate( b(size(a,1), size(a,2), size(from)) )
    do s = 1, size(from)
      if ( from(s) == 0 ) then
        b(:,:,s) = QNAN
      else
        b(:,:,s) = a(:,:,from(s))
      endif
    enddo
    call move_alloc ( b, a )
  end subroutine Morph_LS_Grid_Remap


  !> Where block b starts in the whole grid: its nodes are i0 + 1 .. i0 + nbn,
  !> j0 + 1 .. j0 + nbn.
  pure subroutine block_origin ( grid, b, i0, j0 )
    type(morph_ls_grid_t), intent(in)  :: grid
    integer(I4),           intent(in)  :: b
    integer(I4),           intent(out) :: i0, j0
    i0 = mod(b-1, grid%nblk(1)) * grid%nbn
    j0 = ((b-1)/grid%nblk(1)) * grid%nbn
  end subroutine block_origin


  !> True when full is the whole grid and phi this rank's blocks; else a
  !> MORPH_ERR_SHAPE in the status.
  logical function shapes_ok ( grid, full, phi, what, status ) result ( ok )
    type(morph_ls_grid_t), intent(in)    :: grid
    real(R8),              intent(in)    :: full(:,:), phi(:,:,:)
    character(len=*),      intent(in)    :: what
    type(morph_status_t),  intent(inout) :: status
    character(len=MORPH_MSGLEN) :: msg
    ok = size(full,1) == grid%nblk(1)*grid%nbn .and. size(full,2) == grid%nblk(2)*grid%nbn .and. &
         size(phi,1) == grid%nbn .and. size(phi,2) == grid%nbn .and. size(phi,3) == nslots(grid)
    if ( .not. ok ) then
      write(msg,'(A,A,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'level-set grid, ', what, ': whole grid ', &
            size(full,1), ' x ', size(full,2), ', blocks ', size(phi,1), ' x ', size(phi,2), ' x ', &
            size(phi,3), ' (the grid has ', grid%nbn, ' nodes per block side and ', nslots(grid), &
            ' blocks here)'
      call status%fail ( MORPH_ERR_SHAPE, msg )
    endif
  end function shapes_ok


  !> This rank's active blocks from the whole grid: phi(:,:,s) = the nodes of
  !> block blk(s) of full(nblk(1) nbn, nblk(2) nbn).
  subroutine Morph_LS_Grid_Scatter ( grid, full, phi, status )
    type(morph_ls_grid_t), intent(in)    :: grid
    real(R8),              intent(in)    :: full(:,:)
    real(R8),              intent(inout) :: phi(:,:,:)
    type(morph_status_t),  intent(out)   :: status
    integer(I4) :: s, i0, j0, n

    call status%clear()
    if ( .not. shapes_ok ( grid, full, phi, 'scatter', status ) ) return
    n = grid%nbn
    do s = 1, size(grid%blk)
      call block_origin ( grid, grid%blk(s), i0, j0 )
      phi(:,:,s) = full(i0+1:i0+n, j0+1:j0+n)
    enddo
  end subroutine Morph_LS_Grid_Scatter


  !> This rank's active blocks into the whole grid. Only the nodes of its
  !> active blocks are written, so that every rank can gather into one array
  !> shared in memory.
  subroutine Morph_LS_Grid_Gather ( grid, phi, full, status )
    type(morph_ls_grid_t), intent(in)    :: grid
    real(R8),              intent(in)    :: phi(:,:,:)
    real(R8),              intent(inout) :: full(:,:)
    type(morph_status_t),  intent(out)   :: status
    integer(I4) :: s, i0, j0, n

    call status%clear()
    if ( .not. shapes_ok ( grid, full, phi, 'gather', status ) ) return
    n = grid%nbn
    do s = 1, size(grid%blk)
      call block_origin ( grid, grid%blk(s), i0, j0 )
      full(i0+1:i0+n, j0+1:j0+n) = phi(:,:,s)
    enddo
  end subroutine Morph_LS_Grid_Gather


  !> The ghosts of this rank's blocks, w(1-NG:nbn+NG, 1-NG:nbn+NG, slot), from
  !> their interior nodes and those of the other ranks: the strips for the
  !> other ranks are packed, swapped and unpacked; those between two of my
  !> blocks are copied; at the edge of the grid, extrapolated linearly; toward
  !> an inactive block, by Morph_LS_Ghosts_Upwind from vls(nbn, nbn, slot),
  !> the speed at the nodes, when it is given, else linearly. A collective
  !> (see the header); w or vls of another shape is refused on this rank only.
  subroutine Morph_LS_Halo_Fill ( grid, exch, w, status, vls )
    type(morph_ls_grid_t),   intent(in)    :: grid
    class(morph_exchange_t), intent(inout) :: exch
    real(R8),                intent(inout) :: w(1-NG:, 1-NG:, :)
    type(morph_status_t),    intent(out)   :: status
    real(R8), optional,      intent(in)    :: vls(:,:,:)
    real(R8), allocatable :: sbuf(:), rbuf(:)
    integer(I4) :: m, k, s, d, nb, g(4), r(4), nbs(4)
    character(len=MORPH_MSGLEN) :: msg

    call status%clear()
    if ( size(w,1) /= grid%nbn + 2*NG .or. size(w,2) /= grid%nbn + 2*NG .or. size(w,3) /= nslots(grid) ) then
      write(msg,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'level-set ghosts: w is ', size(w,1), ' x ', size(w,2), &
            ' x ', size(w,3), ' (the grid: ', grid%nbn, ' nodes per block side, ', nslots(grid), &
            ' blocks here, 3 ghost layers)'
      call status%fail ( MORPH_ERR_SHAPE, msg )
      return
    endif
    if ( present(vls) ) then
      if ( size(vls,1) /= grid%nbn .or. size(vls,2) /= grid%nbn .or. size(vls,3) /= nslots(grid) ) then
        write(msg,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'level-set ghosts: V is ', size(vls,1), ' x ', size(vls,2), &
              ' x ', size(vls,3), ' (the grid: ', grid%nbn, ' nodes per block side, ', nslots(grid), &
              ' blocks here)'
        call status%fail ( MORPH_ERR_SHAPE, msg )
        return
      endif
    endif
    m = NG*grid%nbn                                       ! values per strip
    allocate( sbuf(m*size(grid%send,2)), rbuf(sum(grid%rcount)) )
    !$omp parallel do private(k,g)
    do k = 1, size(grid%send,2)
      g = source_box ( grid%nbn, grid%send(2,k) )
      sbuf((k-1)*m+1:k*m) = reshape( w(g(1):g(2), g(3):g(4), grid%send(1,k)), [m] )
    enddo
    !$omp end parallel do
    call exch%swap ( grid%scount, sbuf, grid%rcount, rbuf, status )
    if ( .not. status%ok() ) return

    !$omp parallel do private(s,d,nb,k,g,r,nbs)
    do s = 1, size(grid%blk)
      nbs = neighbours ( grid%nblk, grid%blk(s) )
      do d = 1, 4
        nb = grid%nbr(d,s)
        g  = ghost_box ( grid%nbn, d )
        if ( nb == 0 ) then                               ! the edge of the grid or of the band
          if ( nbs(d) /= 0 .and. present(vls) ) then     ! an inactive block, V given
            call Morph_LS_Ghosts_Upwind ( w(:,:,s), vls(:,:,s), (d+1)/2, 2 - mod(d,2) )
          else
            call Morph_LS_Ghosts_Extrap ( w(:,:,s), (d+1)/2, 2 - mod(d,2) )
          endif
        else if ( grid%owner(nb) == grid%me ) then        ! one of my blocks
          r = source_box ( grid%nbn, d )
          w(g(1):g(2), g(3):g(4), s) = w(r(1):r(2), r(3):r(4), grid%slot(nb))
        else                                              ! received
          k = grid%rpos(d,s)
          w(g(1):g(2), g(3):g(4), s) = reshape( rbuf((k-1)*m+1:k*m), [g(2)-g(1)+1, g(4)-g(3)+1] )
        endif
      enddo
    enddo
    !$omp end parallel do
  end subroutine Morph_LS_Halo_Fill


  !> One step of length dt on this rank's blocks phi(nbn, nbn, slot), V at
  !> their nodes vls(nbn, nbn, slot): Morph_LS_HJ_Step_Full on blocks, the
  !> ghosts filled at every stage by Morph_LS_Halo_Fill with V (toward an
  !> inactive block, Morph_LS_Ghosts_Upwind). A collective.
  !> First every rank checks its arguments and all share the verdicts:
  !> refused on every rank, phi untouched, when one rank refuses (its own
  !> message on that rank, the first refusing rank's code elsewhere): shapes
  !> that are not the grid's (MORPH_ERR_SHAPE), an exchange that is not the
  !> grid's (MORPH_ERR_EXCHANGE), dt < 0 (MORPH_ERR_INPUT). V = 0 on every
  !> block of every rank leaves phi untouched, bit for bit, as on the whole grid.
  subroutine Morph_LS_HJ_Step ( grid, exch, phi, vls, dt, status )
    type(morph_ls_grid_t),   intent(in)    :: grid
    class(morph_exchange_t), intent(inout) :: exch
    real(R8),                intent(inout) :: phi(:,:,:)
    real(R8),                intent(in)    :: vls(:,:,:)
    real(R8),                intent(in)    :: dt
    type(morph_status_t),    intent(out)   :: status
    real(R8), allocatable :: w(:,:,:), u0(:,:,:), us(:,:,:), rate(:,:,:), got(:)
    integer(I4), allocatable :: counts(:)
    integer(I4) :: n, ns, nr, code, p, s, stage
    character(len=MORPH_MSGLEN) :: msg

    call status%clear()
    n  = grid%nbn
    ns = nslots(grid)
    code = MORPH_OK
    msg  = ''
    if ( size(phi,1) /= n .or. size(phi,2) /= n .or. size(phi,3) /= ns .or. &
         size(vls,1) /= n .or. size(vls,2) /= n .or. size(vls,3) /= ns ) then
      code = MORPH_ERR_SHAPE
      write(msg,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'level-set step: phi is ', size(phi,1), ' x ', &
            size(phi,2), ' x ', size(phi,3), ', V ', size(vls,1), ' x ', size(vls,2), ' x ', size(vls,3), &
            ' (the grid: ', n, ' nodes per block side, ', ns, ' blocks here)'
    else if ( exch%rank() /= grid%me .or. exch%nranks() /= grid%nranks ) then
      code = MORPH_ERR_EXCHANGE
      write(msg,'(A,I0,A,I0,A,I0,A,I0,A)') 'level-set step: the exchange is rank ', exch%rank(), ' of ', &
            exch%nranks(), ', the grid rank ', grid%me, ' of ', grid%nranks
    else if ( .not. (dt >= 0.0_R8) ) then
      code = MORPH_ERR_INPUT
      write(msg,'(A,ES10.3,A)') 'level-set step: dt = ', dt, ' (dt >= 0)'
    endif

    ! the verdicts and whether V moves anything, shared by all the ranks
    nr = exch%nranks()
    allocate( counts(0:nr-1), got(2*nr) )
    counts = 2
    call exch%allgather ( [ real(code,R8), merge(1.0_R8, 0.0_R8, any(vls /= 0.0_R8)) ], got, counts, status )
    if ( .not. status%ok() ) return
    if ( code /= MORPH_OK ) then
      call status%fail ( code, msg )
      return
    endif
    do p = 0, nr-1
      if ( got(2*p+1) /= 0.0_R8 ) then
        write(msg,'(A,I0,A,I0,A)') 'level-set step: refused on rank ', p, ' (status ', nint(got(2*p+1)), ')'
        call status%fail ( nint(got(2*p+1), I4), msg )
        return
      endif
    enddo
    if ( all( got(2::2) == 0.0_R8 ) ) return

    allocate( w(1-NG:n+NG, 1-NG:n+NG, ns), rate(n,n,ns) )
    w  = 0.0_R8
    u0 = phi
    us = phi
    do stage = 1, 3
      w(1:n,1:n,:) = us
      call Morph_LS_Halo_Fill ( grid, exch, w, status, vls )
      if ( .not. status%ok() ) return
      !$omp parallel do private(s)
      do s = 1, ns
        call Morph_LS_Rate ( w(:,:,s), vls(:,:,s), grid%h, rate(:,:,s) )
        call Morph_LS_RK3_Combine ( stage, u0(:,:,s), rate(:,:,s), dt, us(:,:,s) )
      enddo
      !$omp end parallel do
    enddo
    phi = us
  end subroutine Morph_LS_HJ_Step

end module Morph_LS_Grid
