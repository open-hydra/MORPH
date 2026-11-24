!>@brief The level-set field on a 2-D grid of blocks split among ranks. The
!> grid has nblk(1) x nblk(2) blocks of nbn x nbn nodes, spacing h; block
!> b = ib + (jb-1) nblk(1) holds the nodes i = (ib-1) nbn + 1 .. ib nbn,
!> j = (jb-1) nbn + 1 .. jb nbn of the whole grid and lives on rank owner(b).
!> A rank keeps its own blocks in block order, its slots: phi(nbn, nbn, slot).
!>
!> THE GHOSTS. A step works on w(1-NG:nbn+NG, 1-NG:nbn+NG, slot), the layout
!> of the kernels of Morph_LS_HJ for each block. Morph_LS_Halo_Fill fills the
!> NG layers along each side with the nodes of the block beyond it (a block
!> has at least NG nodes per side, so one neighbour has them all), and at the
!> edge of the grid with the linear extrapolation of Morph_LS_Ghosts_Extrap,
!> as on the whole grid. The corners are neither filled nor read. Every value
!> a stencil reads is then the one it reads on the whole grid, so
!> Morph_LS_HJ_Step gives the bits of Morph_LS_HJ_Step_Full for any blocks,
!> any split among ranks and any number of threads.
!>
!> THE MESSAGES. A strip is the NG x nbn ghost values along one side of one
!> block, in the order of that block's array (i fastest). Each rank lists the
!> strips it sends and receives from the grid alone, in one order: by rank,
!> then by the block that receives, then by its side (W, E, S, N).
!>
!> THE OTHER RANKS. Morph_LS_HJ_Step and Morph_LS_Halo_Fill are collectives:
!> every rank calls them, in the same order, outside any OpenMP region. The
!> step decides with all the ranks whether it refuses (the arguments of each)
!> and whether it moves anything (V = 0 on every block), so either every rank
!> leaves early or none does. Halo_Fill on its own refuses a w of another
!> shape on that rank only, before the exchange, and the other ranks then
!> wait in it: the step checks the shapes first, together. Init, Scatter and
!> Gather talk to no other rank.
!>
!> Never stops, no I/O: refused arguments return MORPH_ERR_SHAPE,
!> MORPH_ERR_INPUT or MORPH_ERR_EXCHANGE in the status.
module Morph_LS_Grid
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m
  use Morph_Exchange,  only: morph_exchange_t
  use Morph_LS_HJ,     only: MORPH_LS_NG, Morph_LS_Ghosts_Extrap, Morph_LS_Rate, Morph_LS_RK3_Combine

  implicit none
  private

  integer(I4), parameter :: NG = MORPH_LS_NG

  type, public :: morph_ls_grid_t
    integer(I4) :: nbn = 0                   !> nodes along each side of a block
    integer(I4) :: nblk(2) = 0               !> blocks along each axis
    real(R8)    :: h = 0.0_R8                !> node spacing
    integer(I4) :: me = 0                    !> this rank
    integer(I4) :: nranks = 0
    integer(I4), allocatable :: owner(:)     !> the rank of each block
    integer(I4), allocatable :: blk(:)       !> slot -> block: this rank's blocks in block order
    integer(I4), allocatable :: slot(:)      !> block -> its slot here, 0 on another rank
    integer(I4), allocatable :: nbr(:,:)     !> (4, slot): the block beyond side W, E, S, N; 0 at the edge
    integer(I4), allocatable :: rpos(:,:)    !> (4, slot): the strip's place among those received, 0 if not received
    integer(I4), allocatable :: send(:,:)    !> (2, strip): the slot that gives it, the side of the block that receives
    integer(I4), allocatable :: scount(:)    !> (0:nranks-1): values sent to each rank
    integer(I4), allocatable :: rcount(:)    !> (0:nranks-1): values received from each rank
  end type morph_ls_grid_t

  public :: Morph_LS_Grid_Init, Morph_LS_Grid_Scatter, Morph_LS_Grid_Gather
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


  !> The number of this rank's blocks; -1 for a grid not initialised, so that
  !> no shape matches it.
  pure integer(I4) function nslots ( grid )
    type(morph_ls_grid_t), intent(in) :: grid
    nslots = -1
    if ( allocated(grid%blk) ) nslots = size(grid%blk)
  end function nslots


  !> The grid, its split among ranks and the list of the strips this rank
  !> sends and receives. Talks to no other rank: exch gives this rank and their
  !> number. Refused, grid not initialised: nbn < NG or fewer than one block
  !> along an axis, or owner without one value per block (MORPH_ERR_SHAPE);
  !> h <= 0 or a block on a rank that does not exist (MORPH_ERR_INPUT); an
  !> exchange that is not rank 0..n-1 of n >= 1 (MORPH_ERR_EXCHANGE).
  subroutine Morph_LS_Grid_Init ( grid, nbn, nblk, h, owner, exch, status )
    type(morph_ls_grid_t),   intent(out) :: grid
    integer(I4),             intent(in)  :: nbn, nblk(2)
    real(R8),                intent(in)  :: h
    integer(I4),             intent(in)  :: owner(:)
    class(morph_exchange_t), intent(in)  :: exch
    type(morph_status_t),    intent(out) :: status
    integer(I4), allocatable :: tmp(:,:)
    integer(I4) :: nb, me, nr, b, s, d, p, k, nmine, nbs(4)
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
    grid%blk    = pack( [ (b, b = 1, nb) ], owner == me )
    nmine = size(grid%blk)
    allocate( grid%slot(nb), grid%nbr(4,nmine), grid%rpos(4,nmine) )
    allocate( grid%scount(0:nr-1), grid%rcount(0:nr-1) )
    grid%slot   = 0
    grid%rpos   = 0
    grid%scount = 0
    grid%rcount = 0
    do s = 1, nmine
      grid%slot(grid%blk(s)) = s
      grid%nbr(:,s) = neighbours ( nblk, grid%blk(s) )
    enddo

    ! received: from each other rank in rank order, my blocks in order, their sides
    k = 0
    do p = 0, nr-1
      if ( p == me ) cycle
      do s = 1, nmine
        do d = 1, 4
          if ( grid%nbr(d,s) == 0 ) cycle
          if ( owner(grid%nbr(d,s)) /= p ) cycle
          k = k + 1
          grid%rpos(d,s) = k
          grid%rcount(p) = grid%rcount(p) + NG*nbn
        enddo
      enddo
    enddo

    ! sent: to each other rank in rank order, its blocks in order, their sides
    ! beyond which the block is mine (at most one per side of my blocks)
    allocate( tmp(2, 4*nmine) )
    k = 0
    do p = 0, nr-1
      if ( p == me ) cycle
      do b = 1, nb
        if ( owner(b) /= p ) cycle
        nbs = neighbours ( nblk, b )
        do d = 1, 4
          if ( nbs(d) == 0 ) cycle
          if ( owner(nbs(d)) /= me ) cycle
          k = k + 1
          tmp(:,k) = [ grid%slot(nbs(d)), d ]
          grid%scount(p) = grid%scount(p) + NG*nbn
        enddo
      enddo
    enddo
    grid%send = tmp(:,1:k)
  end subroutine Morph_LS_Grid_Init


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


  !> This rank's blocks from the whole grid: phi(:,:,s) = the nodes of block
  !> blk(s) of full(nblk(1) nbn, nblk(2) nbn).
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


  !> This rank's blocks into the whole grid. Only the nodes of its blocks are
  !> written, so that every rank can gather into one array shared in memory.
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
  !> blocks are copied; at the edge of the grid, extrapolated. A collective
  !> (see the header); w of another shape is refused on this rank only.
  subroutine Morph_LS_Halo_Fill ( grid, exch, w, status )
    type(morph_ls_grid_t),   intent(in)    :: grid
    class(morph_exchange_t), intent(inout) :: exch
    real(R8),                intent(inout) :: w(1-NG:, 1-NG:, :)
    type(morph_status_t),    intent(out)   :: status
    real(R8), allocatable :: sbuf(:), rbuf(:)
    integer(I4) :: m, k, s, d, nb, g(4), r(4)
    character(len=MORPH_MSGLEN) :: msg

    call status%clear()
    if ( size(w,1) /= grid%nbn + 2*NG .or. size(w,2) /= grid%nbn + 2*NG .or. size(w,3) /= nslots(grid) ) then
      write(msg,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'level-set ghosts: w is ', size(w,1), ' x ', size(w,2), &
            ' x ', size(w,3), ' (the grid: ', grid%nbn, ' nodes per block side, ', nslots(grid), &
            ' blocks here, 3 ghost layers)'
      call status%fail ( MORPH_ERR_SHAPE, msg )
      return
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

    !$omp parallel do private(s,d,nb,k,g,r)
    do s = 1, size(grid%blk)
      do d = 1, 4
        nb = grid%nbr(d,s)
        g  = ghost_box ( grid%nbn, d )
        if ( nb == 0 ) then                               ! the edge of the grid
          call Morph_LS_Ghosts_Extrap ( w(:,:,s), (d+1)/2, 2 - mod(d,2) )
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
  !> ghosts filled from the neighbouring blocks at every stage. A collective.
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
      call Morph_LS_Halo_Fill ( grid, exch, w, status )
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
