!>@brief Unit test of the narrow band of the level-set grid of blocks
!> (Morph_LS_Grid: the active set, Morph_LS_Grid_Update, Morph_LS_Grid_Remap,
!> and the ghosts and the step on the active blocks only). The reference is
!> an oracle written here on the whole grid: for each active block, the
!> ghosts of a side are the whole grid's nodes beyond it when the block there
!> is active; at the edge of the grid the linear extrapolation of the block's
!> own nodes (Morph_LS_Ghosts_Extrap); toward an inactive block, given V (a
!> fill with V, and every stage of a step), Morph_LS_Ghosts_Upwind, else the
!> linear one; then the rate and the stage of Morph_LS_HJ; the other blocks
!> keep their values. Every comparison is bit for bit (the 64 bits of each
!> value, NaN included):
!>   * the ghosts toward an inactive block (Morph_LS_Ghosts_Upwind), on lines
!>     of integers, exact with or without fused multiply-adds: linear where
!>     the characteristics enter, quadratic where they leave, V = 0 and a
!>     flat line quadratic, both sides of both axes, the other ghosts and the
!>     corners untouched; on a field of reals, where every line enters, the
!>     bits of Morph_LS_Ghosts_Extrap; a V of another shape refused by the
!>     fill; the band edges of the steps' cases hold lines of both kinds;
!>   * Init with an active set: the slots and the strips of Init followed by
!>     an update to that set; refusals; no active block at all;
!>   * the update: the set is what at least one rank asks for, the same on
!>     every rank, the schedule of Init with that set; from(s) is the old slot;
!>     a refusal on one rank is taken by all, the grids untouched; grids of
!>     different sizes are refused;
!>   * the remap: the kept slots bit for bit in their new places, NaN in the
!>     new ones, refusals;
!>   * the band: the ghosts of every block after one fill without V and one
!>     with V are the oracle's, the strips counted from the active pairs
!>     alone; the steps give the
!>     oracle's bits for active sets of every kind (all, none, one block, a
!>     checkerboard, a ring, scattered), blocks of 3, 5 and 8 nodes, 1 to 4
!>     ranks; over several steps the band follows a moving circle, blocks
!>     activated ahead (filled with the distance to it) and dropped behind,
!>     with the requests split among the ranks in two ways.
!> Built without OpenMP, the simulated ranks are skipped, and the output says
!> so. The values printed must be the same at 1 and 8 threads.
!> A serial exchange whose allgather number failat gathers, then reports a
!> failure: a caller must stop on the status, whatever the buffer holds.
module test_band_failing_exchange
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m
  use Morph_Exchange,  only: morph_exchange_serial_t
  implicit none
  private

  type, public, extends(morph_exchange_serial_t) :: failing_exchange_t
    integer(I4) :: failat = 0                 !> the allgather that fails, 0 for none
    integer(I4) :: calls  = 0                 !> the allgathers so far
  contains
    procedure :: allgather => failing_allgather
  end type failing_exchange_t

contains

  subroutine failing_allgather ( self, mine, all, counts, status )
    class(failing_exchange_t), intent(inout) :: self
    real(R8),                  intent(in)    :: mine(:)
    real(R8),                  intent(out)   :: all(:)
    integer(I4),               intent(in)    :: counts(0:)
    type(morph_status_t),      intent(out)   :: status
    self%calls = self%calls + 1
    call self%morph_exchange_serial_t%allgather ( mine, all, counts, status )
    if ( self%calls == self%failat ) &
      call status%fail ( MORPH_ERR_EXCHANGE, 'failing exchange: allgather refused on purpose' )
  end subroutine failing_allgather

end module test_band_failing_exchange


program test_morph_ls_band
  use iso_fortran_env, only: I4 => int32, I8 => int64, R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  !$ use omp_lib, only: omp_get_thread_num
  use Morph_Types_m
  use Morph_Exchange
  use Morph_Exchange_Sim
  use Morph_LS_HJ,   only: MORPH_LS_NG, Morph_LS_Ghosts_Extrap, Morph_LS_Ghosts_Upwind, Morph_LS_Rate, &
                           Morph_LS_RK3_Combine, Morph_LS_HJ_Step_Full
  use Morph_LS_Grid
  use test_band_failing_exchange

  implicit none

  integer(I4), parameter :: NG = MORPH_LS_NG
  real(R8),    parameter :: H  = 1.0_R8 / 32.0_R8
  real(R8),    parameter :: DT = 0.4_R8 * H
  real(R8),    parameter :: SENTINEL = -7.0_R8      ! no field value of the test
  integer(I4), parameter :: NSTEP = 4
  integer(I4), parameter :: KFRONT = 6              ! steps of the band that follows a circle

  integer(I4) :: nfail
  logical     :: omp

  nfail = 0
  omp = .false.
  !$ omp = .true.
  write(*,'(A)') '== Morph level-set narrow band on the grid of blocks =='

  call init_active ()
  call update_one_rank ()
  call update_exchange_fails ()
  call remap ()
  call ghosts_upwind ()
  call halo_one_rank ()
  call steps_one_rank ()
  call front_one_rank ()
  if ( omp ) then
    call update_ranks ()
    call update_refusals ()
    call update_two_refusals ()
    call halo_ranks ()
    call steps_ranks ()
    call front_ranks ()
  else
    write(*,'(A)') '   SKIP the simulated ranks: built without OpenMP'
  endif

  if ( nfail > 0 ) then
    write(*,'(I0,A)') nfail, ' NARROW-BAND CHECK(S) FAILED'
    error stop 1
  endif
  write(*,'(A)') 'RESULT: all narrow-band checks PASS'

contains

  subroutine check ( cond, what )
    logical,          intent(in) :: cond
    character(len=*), intent(in) :: what
    if ( cond ) then
      write(*,'(A,A,A)') '   ', what, ': PASS'
    else
      write(*,'(A,A)') '   FAIL ', what
      nfail = nfail + 1
    endif
  end subroutine check


  !> The 64 bits of every value are the same.
  logical function same ( a, b )
    real(R8), intent(in) :: a(:,:), b(:,:)
    same = size(a,1) == size(b,1) .and. size(a,2) == size(b,2)
    if ( same ) same = all( transfer(a, 0_I8, size(a)) == transfer(b, 0_I8, size(b)) )
  end function same


  logical function same3 ( a, b )
    real(R8), intent(in) :: a(:,:,:), b(:,:,:)
    same3 = size(a,1) == size(b,1) .and. size(a,2) == size(b,2) .and. size(a,3) == size(b,3)
    if ( same3 ) same3 = all( transfer(a, 0_I8, size(a)) == transfer(b, 0_I8, size(b)) )
  end function same3


  real(R8) function qnan ()
    qnan = ieee_value( 0.0_R8, ieee_quiet_nan )
  end function qnan


  !> The test's fields on an nx x ny grid: the distance to a circle away from
  !> the grid's centre, and V from -1 to 1/2 across x (both Godunov branches).
  subroutine fields ( nx, ny, phi, vls )
    integer(I4),           intent(in)  :: nx, ny
    real(R8), allocatable, intent(out) :: phi(:,:), vls(:,:)
    real(R8)    :: xc, yc, rad
    integer(I4) :: i, j
    allocate( phi(nx,ny), vls(nx,ny) )
    xc  = 0.37_R8 * real(nx,R8) * H
    yc  = 0.55_R8 * real(ny,R8) * H
    rad = 0.3_R8 * real(min(nx,ny),R8) * H
    do j = 1, ny
      do i = 1, nx
        phi(i,j) = sqrt( (real(i-1,R8)*H - xc)**2 + (real(j-1,R8)*H - yc)**2 ) - rad
        vls(i,j) = -1.0_R8 + 1.5_R8 * real(i-1,R8) / real(nx-1,R8)
      enddo
    enddo
  end subroutine fields


  !> The whole grid with NaN on the nodes of the inactive blocks.
  function masked ( nbn, nblk, act, full ) result ( u )
    integer(I4), intent(in) :: nbn, nblk(2)
    logical,     intent(in) :: act(:)
    real(R8),    intent(in) :: full(:,:)
    real(R8), allocatable   :: u(:,:)
    integer(I4) :: b, i0, j0
    u = full
    do b = 1, nblk(1)*nblk(2)
      if ( act(b) ) cycle
      i0 = mod(b-1, nblk(1)) * nbn
      j0 = ((b-1)/nblk(1)) * nbn
      u(i0+1:i0+nbn, j0+1:j0+nbn) = qnan()
    enddo
  end function masked


  !> Block b of the whole grid cur with its ghosts as the band must read them:
  !> beyond a side, the whole grid's nodes when the block there is active; at
  !> the edge of the grid the linear extrapolation of the block's own nodes;
  !> toward an inactive block, Morph_LS_Ghosts_Upwind from vfull (the whole
  !> grid's V) when it is given, else the linear one; the corners SENTINEL.
  !> Written here again from the layout of the kernels.
  function oracle_block ( nbn, nblk, act, cur, b, vfull ) result ( wb )
    integer(I4), intent(in) :: nbn, nblk(2), b
    logical,     intent(in) :: act(:)
    real(R8),    intent(in) :: cur(:,:)
    real(R8),    intent(in), optional :: vfull(:,:)
    real(R8), allocatable   :: wb(:,:)
    integer(I4) :: ib, jb, i0, j0, d
    logical :: inside(4), there(4)
    allocate( wb(1-NG:nbn+NG, 1-NG:nbn+NG) )
    wb = SENTINEL
    ib = mod(b-1, nblk(1)) + 1
    jb = (b-1)/nblk(1) + 1
    i0 = (ib-1)*nbn
    j0 = (jb-1)*nbn
    wb(1:nbn, 1:nbn) = cur(i0+1:i0+nbn, j0+1:j0+nbn)
    inside = [ ib > 1, ib < nblk(1), jb > 1, jb < nblk(2) ]
    there  = .false.
    if ( inside(1) ) there(1) = act(b-1)
    if ( inside(2) ) there(2) = act(b+1)
    if ( inside(3) ) there(3) = act(b-nblk(1))
    if ( inside(4) ) there(4) = act(b+nblk(1))
    if ( there(1) ) wb(1-NG:0, 1:nbn)       = cur(i0+1-NG:i0, j0+1:j0+nbn)
    if ( there(2) ) wb(nbn+1:nbn+NG, 1:nbn) = cur(i0+nbn+1:i0+nbn+NG, j0+1:j0+nbn)
    if ( there(3) ) wb(1:nbn, 1-NG:0)       = cur(i0+1:i0+nbn, j0+1-NG:j0)
    if ( there(4) ) wb(1:nbn, nbn+1:nbn+NG) = cur(i0+1:i0+nbn, j0+nbn+1:j0+nbn+NG)
    do d = 1, 4
      if ( there(d) ) cycle
      if ( inside(d) .and. present(vfull) ) then
        call Morph_LS_Ghosts_Upwind ( wb, vfull(i0+1:i0+nbn, j0+1:j0+nbn), (d+1)/2, 2 - mod(d,2) )
      else
        call Morph_LS_Ghosts_Extrap ( wb, (d+1)/2, 2 - mod(d,2) )
      endif
    enddo
  end function oracle_block


  !> nsteps steps of the band on the whole grid u, written here again: each
  !> active block moves with the ghosts of oracle_block; the others keep their
  !> values.
  subroutine oracle_steps ( nbn, nblk, act, vls, u, nsteps )
    integer(I4), intent(in)    :: nbn, nblk(2), nsteps
    logical,     intent(in)    :: act(:)
    real(R8),    intent(in)    :: vls(:,:)
    real(R8),    intent(inout) :: u(:,:)
    real(R8), allocatable :: u0(:,:), cur(:,:), wb(:,:), rate(:,:)
    integer(I4) :: k, stage, b, i0, j0
    allocate( rate(nbn,nbn) )
    do k = 1, nsteps
      u0 = u
      do stage = 1, 3
        cur = u
        do b = 1, nblk(1)*nblk(2)
          if ( .not. act(b) ) cycle
          wb = oracle_block ( nbn, nblk, act, cur, b, vls )
          i0 = mod(b-1, nblk(1)) * nbn
          j0 = ((b-1)/nblk(1)) * nbn
          call Morph_LS_Rate ( wb, vls(i0+1:i0+nbn, j0+1:j0+nbn), H, rate )
          call Morph_LS_RK3_Combine ( stage, u0(i0+1:i0+nbn, j0+1:j0+nbn), rate, DT, &
                                      u(i0+1:i0+nbn, j0+1:j0+nbn) )
        enddo
      enddo
    enddo
  end subroutine oracle_steps


  !> The slots, the strips and the active set of two grids are the same.
  logical function schedule_same ( a, b ) result ( ok )
    type(morph_ls_grid_t), intent(in) :: a, b
    ok = allocated(a%blk) .and. allocated(b%blk)
    if ( .not. ok ) return
    ok = size(a%active) == size(b%active) .and. size(a%blk) == size(b%blk) .and. &
         size(a%slot) == size(b%slot) .and. size(a%send,2) == size(b%send,2) .and. &
         size(a%scount) == size(b%scount) .and. size(a%rcount) == size(b%rcount)
    if ( .not. ok ) return
    ok = all( a%active .eqv. b%active ) .and. all( a%blk == b%blk ) .and. all( a%slot == b%slot ) .and. &
         all( a%nbr == b%nbr ) .and. all( a%rpos == b%rpos ) .and. all( a%send == b%send ) .and. &
         all( a%scount == b%scount ) .and. all( a%rcount == b%rcount )
  end function schedule_same


  !> The ghosts of this rank's slots after one fill are the oracle's for their
  !> blocks: the neighbours' nodes where the neighbour is active, the
  !> extrapolation elsewhere (from vfull, the whole grid's V, when the fill had
  !> V), the corners untouched.
  logical function ghosts_ok ( grid, full, w, vfull ) result ( ok )
    type(morph_ls_grid_t), intent(in) :: grid
    real(R8),              intent(in) :: full(:,:), w(1-NG:,1-NG:,:)
    real(R8), optional,    intent(in) :: vfull(:,:)
    integer(I4) :: s
    ok = .true.
    do s = 1, size(grid%blk)
      if ( .not. same( w(:,:,s), oracle_block ( grid%nbn, grid%nblk, grid%active, full, grid%blk(s), vfull ) ) ) &
        ok = .false.
    enddo
  end function ghosts_ok


  !> The values each rank sends to and receives from each other rank, counted
  !> here from the pairs of active blocks across a side alone.
  logical function counts_ok ( grid ) result ( ok )
    type(morph_ls_grid_t), intent(in) :: grid
    integer(I4) :: snd(0:grid%nranks-1), rcv(0:grid%nranks-1), b, c, ib, jb, ic, jc, nb
    snd = 0
    rcv = 0
    nb = grid%nblk(1)*grid%nblk(2)
    do b = 1, nb
      do c = 1, nb
        if ( .not. ( grid%active(b) .and. grid%active(c) ) ) cycle
        ib = mod(b-1, grid%nblk(1)); jb = (b-1)/grid%nblk(1)
        ic = mod(c-1, grid%nblk(1)); jc = (c-1)/grid%nblk(1)
        if ( abs(ib-ic) + abs(jb-jc) /= 1 ) cycle                 ! c beyond one side of b
        if ( grid%owner(b) == grid%me .and. grid%owner(c) /= grid%me ) &
          rcv(grid%owner(c)) = rcv(grid%owner(c)) + NG*grid%nbn   ! b's ghosts from c
        if ( grid%owner(c) == grid%me .and. grid%owner(b) /= grid%me ) &
          snd(grid%owner(b)) = snd(grid%owner(b)) + NG*grid%nbn   ! c's nodes to b
      enddo
    enddo
    ok = all( grid%scount == snd ) .and. all( grid%rcount == rcv )
  end function counts_ok


  !> One rank's part of a band of fixed active set: its grid, its blocks of
  !> phi0 and v0, one fill of the ghosts without V and one with V (halo: both
  !> are the oracle's), the strips counted (cnt), nsteps steps, its blocks
  !> gathered into out.
  subroutine band_advance ( exch, nbn, nblk, owner, act, phi0, v0, nsteps, out, halo, cnt, st )
    class(morph_exchange_t), intent(inout) :: exch
    integer(I4),             intent(in)    :: nbn, nblk(2), owner(:), nsteps
    logical,                 intent(in)    :: act(:)
    real(R8),                intent(in)    :: phi0(:,:), v0(:,:)
    real(R8),                intent(inout) :: out(:,:)
    logical,                 intent(out)   :: halo, cnt
    type(morph_status_t),    intent(out)   :: st
    type(morph_ls_grid_t) :: grid
    real(R8), allocatable :: phi(:,:,:), vb(:,:,:), w(:,:,:)
    integer(I4) :: k, ns

    halo = .false.
    cnt  = .false.
    call Morph_LS_Grid_Init ( grid, nbn, nblk, H, owner, exch, st, act )
    if ( .not. st%ok() ) return
    ns = size(grid%blk)
    allocate( phi(nbn,nbn,ns), vb(nbn,nbn,ns), w(1-NG:nbn+NG, 1-NG:nbn+NG, ns) )
    call Morph_LS_Grid_Scatter ( grid, phi0, phi, st )
    if ( .not. st%ok() ) return
    call Morph_LS_Grid_Scatter ( grid, v0, vb, st )
    if ( .not. st%ok() ) return
    w = SENTINEL
    w(1:nbn,1:nbn,:) = phi
    call Morph_LS_Halo_Fill ( grid, exch, w, st )
    if ( .not. st%ok() ) return
    halo = ghosts_ok ( grid, phi0, w )
    w = SENTINEL
    w(1:nbn,1:nbn,:) = phi
    call Morph_LS_Halo_Fill ( grid, exch, w, st, vb )
    if ( .not. st%ok() ) return
    halo = halo .and. ghosts_ok ( grid, phi0, w, v0 )
    cnt  = counts_ok ( grid )
    do k = 1, nsteps
      call Morph_LS_HJ_Step ( grid, exch, phi, vb, DT, st )
      if ( .not. st%ok() ) return
    enddo
    call Morph_LS_Grid_Gather ( grid, phi, out, st )
  end subroutine band_advance


  ! --- the band that follows a moving circle --------------------------------

  !> The signed distance to the circle of step k at node (i,j): the circle
  !> grows by .6 block per step from .25 of the grid's smaller side.
  real(R8) function dist_k ( nbn, nx, ny, k, i, j )
    integer(I4), intent(in) :: nbn, nx, ny, k, i, j
    real(R8) :: xc, yc, rad
    xc  = 0.45_R8 * real(nx,R8) * H
    yc  = 0.40_R8 * real(ny,R8) * H
    rad = 0.25_R8 * real(min(nx,ny),R8) * H + 0.6_R8 * real(k*nbn,R8) * H
    dist_k = sqrt( (real(i-1,R8)*H - xc)**2 + (real(j-1,R8)*H - yc)**2 ) - rad
  end function dist_k


  !> The blocks within 1.5 blocks of the circle of step k (min |distance| over
  !> their nodes).
  function wanted_k ( nbn, nblk, k ) result ( w )
    integer(I4), intent(in) :: nbn, nblk(2), k
    logical :: w(nblk(1)*nblk(2))
    integer(I4) :: b, i0, j0, i, j
    real(R8) :: dmin
    do b = 1, nblk(1)*nblk(2)
      i0 = mod(b-1, nblk(1)) * nbn
      j0 = ((b-1)/nblk(1)) * nbn
      dmin = huge(1.0_R8)
      do j = j0+1, j0+nbn
        do i = i0+1, i0+nbn
          dmin = min( dmin, abs( dist_k ( nbn, nblk(1)*nbn, nblk(2)*nbn, k, i, j ) ) )
        enddo
      enddo
      w(b) = dmin < 1.5_R8 * real(nbn,R8) * H
    enddo
  end function wanted_k


  !> The oracle of the band that follows the circle: at step k the blocks of
  !> wanted_k move; a block that joins takes the distance of step k, one that
  !> leaves becomes NaN; then one step.
  function oracle_front ( nbn, nblk, v0 ) result ( u )
    integer(I4), intent(in) :: nbn, nblk(2)
    real(R8),    intent(in) :: v0(:,:)
    real(R8), allocatable   :: u(:,:)
    logical :: act(nblk(1)*nblk(2)), nxt(nblk(1)*nblk(2))
    integer(I4) :: k, b, i0, j0, i, j, nx, ny
    nx = nblk(1)*nbn
    ny = nblk(2)*nbn
    allocate( u(nx,ny) )
    u   = qnan()
    act = .false.
    do k = 0, KFRONT-1
      nxt = wanted_k ( nbn, nblk, k )
      do b = 1, nblk(1)*nblk(2)
        i0 = mod(b-1, nblk(1)) * nbn
        j0 = ((b-1)/nblk(1)) * nbn
        if ( nxt(b) .and. .not. act(b) ) then
          do j = j0+1, j0+nbn
            do i = i0+1, i0+nbn
              u(i,j) = dist_k ( nbn, nx, ny, k, i, j )
            enddo
          enddo
        else if ( act(b) .and. .not. nxt(b) ) then
          u(i0+1:i0+nbn, j0+1:j0+nbn) = qnan()
        endif
      enddo
      act = nxt
      call oracle_steps ( nbn, nblk, act, v0, u, 1 )
    enddo
  end function oracle_front


  !> One rank's part of the band that follows the circle: a grid with no
  !> active block, then at each step its requests (mode 1: the wanted blocks
  !> it owns and the wanted neighbours of its active blocks; mode 2: the
  !> wanted blocks b with mod(b + rank, n) = 0, whoever owns them), the
  !> update, the remap of phi, the new slots filled with the distance of the
  !> step, V scattered again, one step; at the end its blocks into out.
  subroutine front_advance ( exch, nbn, nblk, owner, mode, v0, out, filled, st )
    class(morph_exchange_t), intent(inout) :: exch
    integer(I4),             intent(in)    :: nbn, nblk(2), owner(:), mode
    real(R8),                intent(in)    :: v0(:,:)
    real(R8),                intent(inout) :: out(:,:)
    logical,                 intent(out)   :: filled
    type(morph_status_t),    intent(out)   :: st
    type(morph_ls_grid_t) :: grid
    real(R8),    allocatable :: phi(:,:,:), vb(:,:,:)
    integer(I4), allocatable :: from(:)
    logical :: wk(nblk(1)*nblk(2)), want(nblk(1)*nblk(2)), none(nblk(1)*nblk(2))
    integer(I4) :: k, b, s, d, i0, j0, i, j, nx, ny, me, nr, nbs(4)

    filled = .true.
    nx = nblk(1)*nbn
    ny = nblk(2)*nbn
    none = .false.
    call Morph_LS_Grid_Init ( grid, nbn, nblk, H, owner, exch, st, none )
    if ( .not. st%ok() ) return
    me = exch%rank()
    nr = exch%nranks()
    allocate( phi(nbn,nbn,0) )
    do k = 0, KFRONT-1
      wk = wanted_k ( nbn, nblk, k )
      want = .false.
      do b = 1, nblk(1)*nblk(2)
        if ( .not. wk(b) ) cycle
        if ( mode == 1 ) then
          if ( owner(b) == me ) want(b) = .true.
          nbs = 0
          if ( mod(b-1, nblk(1)) > 0 )       nbs(1) = b - 1
          if ( mod(b-1, nblk(1)) < nblk(1)-1 ) nbs(2) = b + 1
          if ( (b-1)/nblk(1) > 0 )           nbs(3) = b - nblk(1)
          if ( (b-1)/nblk(1) < nblk(2)-1 )   nbs(4) = b + nblk(1)
          do d = 1, 4
            if ( nbs(d) == 0 ) cycle
            if ( owner(nbs(d)) == me .and. grid%active(nbs(d)) ) want(b) = .true.
          enddo
        else
          if ( mod(b + me, nr) == 0 ) want(b) = .true.
        endif
      enddo
      call Morph_LS_Grid_Update ( grid, exch, want, from, st )
      if ( .not. st%ok() ) return
      call Morph_LS_Grid_Remap ( from, phi, st )
      if ( .not. st%ok() ) return
      do s = 1, size(grid%blk)
        if ( from(s) /= 0 ) cycle
        if ( .not. all( ieee_is_nan( phi(:,:,s) ) ) ) filled = .false.
        b  = grid%blk(s)
        i0 = mod(b-1, nblk(1)) * nbn
        j0 = ((b-1)/nblk(1)) * nbn
        do j = 1, nbn
          do i = 1, nbn
            phi(i,j,s) = dist_k ( nbn, nx, ny, k, i0+i, j0+j )
          enddo
        enddo
      enddo
      if ( allocated(vb) ) deallocate( vb )
      allocate( vb(nbn,nbn,size(grid%blk)) )
      call Morph_LS_Grid_Scatter ( grid, v0, vb, st )
      if ( .not. st%ok() ) return
      call Morph_LS_HJ_Step ( grid, exch, phi, vb, DT, st )
      if ( .not. st%ok() ) return
    enddo
    call Morph_LS_Grid_Gather ( grid, phi, out, st )
  end subroutine front_advance


  ! ===========================================================================
  subroutine init_active ()
    type(morph_exchange_serial_t) :: ex
    type(morph_ls_grid_t) :: ga, gb
    type(morph_status_t)  :: st, st2
    integer(I4), allocatable :: from(:)
    logical  :: act(12), none(12)
    real(R8) :: phi(5,5,0), vb(5,5,0), full(20,15)
    integer(I4) :: k

    do k = 1, 12
      act(k) = mod(k, 3) /= 1
    enddo
    call Morph_LS_Grid_Init ( ga, 5, [4,3], H, [ (0, k = 1, 12) ], ex, st, act )
    call Morph_LS_Grid_Init ( gb, 5, [4,3], H, [ (0, k = 1, 12) ], ex, st2 )
    call check ( st%ok() .and. st2%ok() .and. size(gb%blk) == 12, 'Init without an active set: every block active' )
    call Morph_LS_Grid_Update ( gb, ex, act, from, st2 )
    call check ( st%ok() .and. st2%ok() .and. schedule_same ( ga, gb ) .and. size(ga%blk) == count(act), &
                 'Init with an active set: the slots and strips of Init then an update to it' )

    call Morph_LS_Grid_Init ( ga, 5, [4,3], H, [ (0, k = 1, 12) ], ex, st, act(1:11) )
    call check ( st%code == MORPH_ERR_SHAPE .and. index(st%message, '11 active flags for 12 blocks') > 0 .and. &
                 .not. allocated(ga%blk), 'Init with 11 active flags for 12 blocks: refused, no grid' )

    none = .false.
    full = 1.0_R8
    call Morph_LS_Grid_Init ( ga, 5, [4,3], H, [ (0, k = 1, 12) ], ex, st, none )
    call check ( st%ok() .and. size(ga%blk) == 0 .and. all( ga%slot == 0 ) .and. size(ga%send,2) == 0, &
                 'Init with no active block: no slot' )
    call Morph_LS_Grid_Scatter ( ga, full, phi, st )
    vb = 0.0_R8
    if ( st%ok() ) call Morph_LS_HJ_Step ( ga, ex, phi, vb, DT, st )
    if ( st%ok() ) call Morph_LS_Grid_Gather ( ga, phi, full, st )
    call check ( st%ok() .and. all( full == 1.0_R8 ), &
                 'no active block: scatter, the step and gather accept no slot and write nothing' )
  end subroutine init_active


  ! ===========================================================================
  subroutine update_one_rank ()
    type(morph_exchange_serial_t) :: ex
    type(morph_sim_world_t), target :: world
    type(morph_exchange_sim_t) :: other
    type(morph_ls_grid_t) :: g, keep, gi, gnot
    type(morph_status_t)  :: st
    integer(I4), allocatable :: from(:), expect(:), oldblk(:)
    logical :: w1(20), w2(20), none(20)
    integer(I4) :: k, s, b

    call Morph_LS_Grid_Init ( g, 4, [5,4], H, [ (0, k = 1, 20) ], ex, st )
    do k = 1, 20
      w1(k) = mod(3*k, 7) < 4
      w2(k) = mod(5*k, 6) < 3
    enddo
    call Morph_LS_Grid_Update ( g, ex, w1, from, st )
    call check ( st%ok() .and. all( g%active .eqv. w1 ) .and. size(from) == count(w1) .and. &
                 all( from == g%blk ), 'update, 1 rank: the set asked for; from(s) is the block''s old slot' )
    oldblk = g%blk
    call Morph_LS_Grid_Update ( g, ex, w2, from, st )
    allocate( expect(size(g%blk)) )
    expect = 0
    do s = 1, size(g%blk)
      do k = 1, size(oldblk)
        if ( oldblk(k) == g%blk(s) ) expect(s) = k
      enddo
    enddo
    call check ( st%ok() .and. all( g%active .eqv. w2 ) .and. size(from) == size(expect) .and. &
                 all( from == expect ) .and. any( expect == 0 ) .and. any( expect > 0 ) .and. &
                 count( w1 .and. .not. w2 ) > 0, &
                 'update, 1 rank: blocks kept (their old slots), added (0) and dropped' )
    call Morph_LS_Grid_Init ( gi, 4, [5,4], H, [ (0, k = 1, 20) ], ex, st, w2 )
    call check ( st%ok() .and. schedule_same ( g, gi ), 'update, 1 rank: the schedule of Init with that set' )
    keep = g
    call Morph_LS_Grid_Update ( g, ex, w2, from, st )
    call check ( st%ok() .and. schedule_same ( g, keep ) .and. all( from == [ (s, s = 1, size(g%blk)) ] ), &
                 'update, 1 rank: the same set again, the same schedule and from(s) = s' )

    call Morph_LS_Grid_Update ( g, ex, w2(1:19), from, st )
    call check ( st%code == MORPH_ERR_SHAPE .and. index(st%message, '19 requests for 20 blocks') > 0 .and. &
                 schedule_same ( g, keep ) .and. all( from == [ (s, s = 1, size(g%blk)) ] ), &
                 'update refused, 19 requests for 20 blocks: grid untouched, from(s) = s' )
    call Morph_LS_Grid_Update ( gnot, ex, w2, from, st )
    call check ( st%code == MORPH_ERR_SHAPE .and. index(st%message, 'not initialised') > 0 .and. size(from) == 0, &
                 'update refused, a grid not initialised: no slot in from' )
    call world%init ( 2, st )
    other%me = 1
    other%world => world
    call Morph_LS_Grid_Init ( gnot, 4, [5,4], H, [ (mod(k,2), k = 1, 20) ], other, st )
    keep = gnot
    call Morph_LS_Grid_Update ( gnot, ex, w2, from, st )
    call check ( st%code == MORPH_ERR_EXCHANGE .and. index(st%message, 'the grid rank 1 of 2') > 0 .and. &
                 schedule_same ( gnot, keep ) .and. all( from == [ (s, s = 1, size(gnot%blk)) ] ), &
                 'update refused, the exchange of another rank: grid untouched' )
    other%me = 0
    call Morph_LS_Grid_Init ( gnot, 4, [5,4], H, [ (mod(k,2), k = 1, 20) ], other, st )
    keep = gnot
    call Morph_LS_Grid_Update ( gnot, ex, w2, from, st )
    call check ( st%code == MORPH_ERR_EXCHANGE .and. index(st%message, 'rank 0 of 1, the grid rank 0 of 2') > 0 &
                 .and. schedule_same ( gnot, keep ), &
                 'update refused, an exchange of the grid''s rank among another number of ranks' )

    none = .false.
    call Morph_LS_Grid_Update ( g, ex, none, from, st )
    call check ( st%ok() .and. size(g%blk) == 0 .and. size(from) == 0 .and. .not. any( g%active ), &
                 'update, 1 rank: nothing asked, no block active' )
    call Morph_LS_Grid_Update ( g, ex, w1, from, st )
    b = 0
    if ( size(from) > 0 ) b = maxval(from)
    call check ( st%ok() .and. all( g%active .eqv. w1 ) .and. size(from) == count(w1) .and. b == 0, &
                 'update, 1 rank: from no block to some, every slot new' )
  end subroutine update_one_rank


  ! ===========================================================================
  !> The exchange fails in the update's first allgather, then in its second:
  !> the failure is returned, the grid untouched, from(s) = s.
  subroutine update_exchange_fails ()
    type(failing_exchange_t) :: fx
    type(morph_ls_grid_t) :: g, keep
    type(morph_status_t)  :: st
    integer(I4), allocatable :: from(:)
    logical :: w(12)
    integer(I4) :: k, s
    character(len=6), parameter :: nth(2) = [ 'first ', 'second' ]

    w = .true.
    w(3) = .false.
    do k = 1, 2
      fx%failat = 0
      call Morph_LS_Grid_Init ( g, 4, [4,3], H, [ (0, s = 1, 12) ], fx, st )
      keep = g
      fx%calls  = 0
      fx%failat = k
      call Morph_LS_Grid_Update ( g, fx, w, from, st )
      call check ( st%code == MORPH_ERR_EXCHANGE .and. index(st%message, 'on purpose') > 0 .and. fx%calls == k &
                   .and. schedule_same ( g, keep ) .and. all( from == [ (s, s = 1, size(g%blk)) ] ), &
                   'update, the exchange fails in its '//trim(nth(k))//' allgather: returned, grid untouched' )
    enddo
  end subroutine update_exchange_fails


  !> Ranks 0 and 2 of 3 refuse with different codes: each reports its own,
  !> rank 1 the first refusing rank's, rank 0's.
  subroutine update_two_refusals ()
    type(morph_sim_world_t), target :: world
    type(morph_status_t) :: st
    integer(I4) :: codes(0:2), owner(12), k
    logical :: named(0:2)

    do k = 1, 12
      owner(k) = mod(k, 3)
    enddo
    call world%init ( 3, st )
    codes = -1
    named = .false.
    !$omp parallel num_threads(3)
    block
      type(morph_exchange_sim_t) :: ex, oth
      type(morph_ls_grid_t) :: g
      type(morph_status_t) :: s
      integer(I4), allocatable :: from(:)
      logical :: want(12)
      integer(I4) :: t
      t = 0
      !$ t = omp_get_thread_num()
      ex%me = t
      ex%world => world
      oth%me = 1
      oth%world => world
      want = .true.
      if ( t == 2 ) then
        call Morph_LS_Grid_Init ( g, 4, [4,3], H, owner, oth, s )
      else
        call Morph_LS_Grid_Init ( g, 4, [4,3], H, owner, ex, s )
      endif
      if ( t == 0 ) then
        call Morph_LS_Grid_Update ( g, ex, want(1:10), from, s )
      else
        call Morph_LS_Grid_Update ( g, ex, want, from, s )
      endif
      codes(t) = s%code
      select case ( t )
      case ( 0 )
        named(t) = index(s%message, '10 requests for 12 blocks') > 0
      case ( 1 )
        named(t) = index(s%message, 'refused on rank 0 (status') > 0
      case default
        named(t) = index(s%message, 'the grid rank 1 of 3') > 0
      end select
    end block
    !$omp end parallel
    call check ( st%ok() .and. all( codes == [ MORPH_ERR_SHAPE, MORPH_ERR_SHAPE, MORPH_ERR_EXCHANGE ] ) .and. &
                 all( named ), 'update refused on ranks 0 and 2 with different codes: rank 1 reports rank 0''s' )
  end subroutine update_two_refusals


  ! ===========================================================================
  subroutine remap ()
    type(morph_status_t)  :: st
    real(R8), allocatable :: a(:,:,:), a0(:,:,:), b(:,:,:)
    integer(I4) :: i, j, s
    logical :: kept, nans

    allocate( a(3,2,4) )
    do s = 1, 4
      do j = 1, 2
        do i = 1, 3
          a(i,j,s) = real(100*s + 10*j + i, R8) / 7.0_R8
        enddo
      enddo
    enddo
    a(2,1,3) = -0.0_R8
    a0 = a
    call Morph_LS_Grid_Remap ( [ 3, 0, 1, 4, 0 ], a, st )
    kept = st%ok() .and. size(a,3) == 5
    if ( kept ) kept = same3( a(:,:,1:1), a0(:,:,3:3) ) .and. same3( a(:,:,3:3), a0(:,:,1:1) ) .and. &
                       same3( a(:,:,4:4), a0(:,:,4:4) )
    call check ( kept, 'remap: the kept slots bit for bit in their new places, a -0 too' )
    nans = st%ok() .and. size(a,3) == 5
    if ( nans ) nans = all( ieee_is_nan( a(:,:,2) ) ) .and. all( ieee_is_nan( a(:,:,5) ) )
    call check ( nans, 'remap: a NaN in every node of the new slots' )

    b = a0
    call Morph_LS_Grid_Remap ( [ 1, 5 ], b, st )
    call check ( st%code == MORPH_ERR_SHAPE .and. index(st%message, 'slots 1 .. 5 of a field of 4 slots') > 0 .and. &
                 same3( b, a0 ), 'remap refused, slot 5 of a field of 4: the field untouched' )
    call Morph_LS_Grid_Remap ( [ 2, -1 ], b, st )
    call check ( st%code == MORPH_ERR_SHAPE .and. same3( b, a0 ), 'remap refused, slot -1: the field untouched' )
    deallocate( b )
    call Morph_LS_Grid_Remap ( [ 1 ], b, st )
    call check ( st%code == MORPH_ERR_SHAPE .and. index(st%message, 'not allocated') > 0, &
                 'remap refused, a field not allocated' )
    b = a0
    call Morph_LS_Grid_Remap ( [ integer(I4) :: ], b, st )
    call check ( st%ok() .and. size(b,1) == 3 .and. size(b,2) == 2 .and. size(b,3) == 0, 'remap to no slot' )
  end subroutine remap


  ! ===========================================================================
  !> Morph_LS_Ghosts_Upwind on its own; the fill's refusal of a V of another
  !> shape; the kinds of the band edges in the steps' cases.
  subroutine ghosts_upwind ()
    type(morph_exchange_serial_t) :: ex
    type(morph_ls_grid_t) :: g
    type(morph_status_t)  :: st
    real(R8) :: base(1-NG:4+NG, 1-NG:4+NG), baset(1-NG:4+NG, 1-NG:4+NG), a(1-NG:4+NG, 1-NG:4+NG)
    real(R8) :: r(1-NG:4+NG, 1-NG:4+NG), v(4,4), lo(NG,4), hi(NG,4)
    real(R8), allocatable :: phi0(:,:), v0(:,:), u8(:,:), e8(:,:), x8(:,:), v8(:,:), w(:,:,:), vbad(:,:,:)
    integer(I4) :: i, j, k, axis, side, nin, nout, nin2, nout2
    logical :: ok

    ! lines 1-3 (i-1)^2, line 4 a flat start; V at the two edge nodes -1, +1,
    ! 0, -1, the opposite sign at the nodes next to them (V is read at the
    ! edge node only)
    base = SENTINEL
    do j = 1, 4
      do i = 1, 4
        base(i,j) = real((i-1)**2, R8)
      enddo
    enddo
    base(1:4,4) = [ 5.0_R8, 5.0_R8, 7.0_R8, 12.0_R8 ]
    v(:,1) = [ -1.0_R8,  1.0_R8,  1.0_R8, -1.0_R8 ]
    v(:,2) = [  1.0_R8, -1.0_R8, -1.0_R8,  1.0_R8 ]
    v(:,3) =    0.0_R8
    v(:,4) = [ -1.0_R8,  1.0_R8,  1.0_R8, -1.0_R8 ]
    ! low side: leaves, enters, V = 0, flat (V d = -0); high side: enters, leaves, V = 0, enters
    lo(:,1) = [  1.0_R8,  4.0_R8,  9.0_R8 ]
    lo(:,2) = [ -1.0_R8, -2.0_R8, -3.0_R8 ]
    lo(:,3) = [  1.0_R8,  4.0_R8,  9.0_R8 ]
    lo(:,4) = [  7.0_R8, 11.0_R8, 17.0_R8 ]
    hi(:,1) = [ 14.0_R8, 19.0_R8, 24.0_R8 ]
    hi(:,2) = [ 16.0_R8, 25.0_R8, 36.0_R8 ]
    hi(:,3) = [ 16.0_R8, 25.0_R8, 36.0_R8 ]
    hi(:,4) = [ 17.0_R8, 22.0_R8, 27.0_R8 ]
    baset = SENTINEL
    baset(1:4,1:4) = transpose( base(1:4,1:4) )
    do axis = 1, 2
      do side = 1, 2
        if ( axis == 1 ) then
          a = base
          r = base
          call Morph_LS_Ghosts_Upwind ( a, v, 1, side )
          do k = 1, NG
            if ( side == 1 ) then
              r(1-k,1:4) = lo(k,:)
            else
              r(4+k,1:4) = hi(k,:)
            endif
          enddo
        else
          a = baset
          r = baset
          call Morph_LS_Ghosts_Upwind ( a, transpose(v), 2, side )
          do k = 1, NG
            if ( side == 1 ) then
              r(1:4,1-k) = lo(k,:)
            else
              r(1:4,4+k) = hi(k,:)
            endif
          enddo
        endif
        call check ( same( a, r ), 'upwind ghosts on integers, axis ' // achar(48+axis) // ' side ' // &
                     achar(48+side) // ': linear where the line enters, quadratic where it leaves, V = 0 and ' // &
                     'a flat line quadratic, the rest untouched' )
      enddo
    enddo

    ! on reals, every line entering: the bits of the linear extrapolation
    call fields ( 8, 8, phi0, v0 )
    allocate( u8(1-NG:8+NG, 1-NG:8+NG), v8(8,8) )
    u8 = SENTINEL
    u8(1:8,1:8) = phi0
    ok = .true.
    do axis = 1, 2
      do side = 1, 2
        v8 = 1.0_R8
        do j = 1, 8
          if ( axis == 1 .and. side == 1 ) v8(1,j) =  sign( 1.0_R8, phi0(2,j) - phi0(1,j) )
          if ( axis == 1 .and. side == 2 ) v8(8,j) = -sign( 1.0_R8, phi0(8,j) - phi0(7,j) )
          if ( axis == 2 .and. side == 1 ) v8(j,1) =  sign( 1.0_R8, phi0(j,2) - phi0(j,1) )
          if ( axis == 2 .and. side == 2 ) v8(j,8) = -sign( 1.0_R8, phi0(j,8) - phi0(j,7) )
        enddo
        e8 = u8
        x8 = u8
        call Morph_LS_Ghosts_Upwind ( e8, v8, axis, side )
        call Morph_LS_Ghosts_Extrap ( x8, axis, side )
        if ( .not. same( e8, x8 ) ) ok = .false.
      enddo
    enddo
    call check ( ok, 'upwind ghosts on reals, every line entering: the bits of the linear extrapolation' )

    ! the fill refuses a V of another shape along each of its axes, w untouched
    call Morph_LS_Grid_Init ( g, 3, [2,2], H, [0,0,0,0], ex, st, [ .true., .false., .true., .true. ] )
    allocate( w(1-NG:3+NG, 1-NG:3+NG, 3) )
    w = SENTINEL
    ok = .true.
    do k = 1, 3
      if ( allocated(vbad) ) deallocate( vbad )
      if ( k == 1 ) allocate( vbad(2,3,3) )
      if ( k == 2 ) allocate( vbad(3,2,3) )
      if ( k == 3 ) allocate( vbad(3,3,2) )
      vbad = -1.0_R8
      call Morph_LS_Halo_Fill ( g, ex, w, st, vbad )
      if ( .not. ( st%code == MORPH_ERR_SHAPE .and. index(st%message, 'level-set ghosts: V is') > 0 ) ) ok = .false.
    enddo
    call check ( ok .and. index(st%message, 'V is 3 x 3 x 2') > 0 .and. all( w == SENTINEL ), &
                 'ghosts with a V of 2 x 3 x 3, 3 x 2 x 3 or 3 x 3 x 2 for 3 x 3 x 3: refused, w untouched' )

    ! the band edges of two steps' cases hold lines that enter and lines that leave
    call fields ( 30, 25, phi0, v0 )
    call edge_kinds ( 5, [6,5], pattern(6, [6,5]), phi0, v0, nin, nout )
    call fields ( 18, 15, phi0, v0 )
    call edge_kinds ( 3, [6,5], pattern(4, [6,5]), phi0, v0, nin2, nout2 )
    call check ( nin > 0 .and. nout > 0 .and. nin2 > 0 .and. nout2 > 0, &
                 'the band edges of the scattered blocks of 5 and of the checkerboard of 3 hold lines of both kinds' )
  end subroutine ghosts_upwind


  !> The lines across the sides of the active blocks that face an inactive
  !> block inside the grid: nin where the characteristics enter (the rule of
  !> Morph_LS_Ghosts_Upwind, V at the edge node), nout where they leave.
  subroutine edge_kinds ( nbn, nblk, act, phi, v, nin, nout )
    integer(I4), intent(in)  :: nbn, nblk(2)
    logical,     intent(in)  :: act(:)
    real(R8),    intent(in)  :: phi(:,:), v(:,:)
    integer(I4), intent(out) :: nin, nout
    integer(I4) :: b, ib, jb, i0, j0, t, e, q
    logical :: inflow
    nin  = 0
    nout = 0
    do b = 1, nblk(1)*nblk(2)
      if ( .not. act(b) ) cycle
      ib = mod(b-1, nblk(1)) + 1
      jb = (b-1)/nblk(1) + 1
      i0 = (ib-1)*nbn
      j0 = (jb-1)*nbn
      do t = 1, nbn
        if ( ib > 1 ) then
          if ( .not. act(b-1) ) then
            e = i0+1; q = j0+t
            inflow = v(e,q) * (phi(e+1,q) - phi(e,q)) > 0.0_R8
            if ( inflow ) then; nin = nin + 1; else; nout = nout + 1; endif
          endif
        endif
        if ( ib < nblk(1) ) then
          if ( .not. act(b+1) ) then
            e = i0+nbn; q = j0+t
            inflow = v(e,q) * (phi(e,q) - phi(e-1,q)) < 0.0_R8
            if ( inflow ) then; nin = nin + 1; else; nout = nout + 1; endif
          endif
        endif
        if ( jb > 1 ) then
          if ( .not. act(b-nblk(1)) ) then
            e = j0+1; q = i0+t
            inflow = v(q,e) * (phi(q,e+1) - phi(q,e)) > 0.0_R8
            if ( inflow ) then; nin = nin + 1; else; nout = nout + 1; endif
          endif
        endif
        if ( jb < nblk(2) ) then
          if ( .not. act(b+nblk(1)) ) then
            e = j0+nbn; q = i0+t
            inflow = v(q,e) * (phi(q,e) - phi(q,e-1)) < 0.0_R8
            if ( inflow ) then; nin = nin + 1; else; nout = nout + 1; endif
          endif
        endif
      enddo
    enddo
  end subroutine edge_kinds


  !> Active sets on a grid of nblk blocks: all, none, one, a checkerboard, a
  !> ring, scattered.
  function pattern ( kind, nblk ) result ( act )
    integer(I4), intent(in) :: kind, nblk(2)
    logical :: act(nblk(1)*nblk(2))
    integer(I4) :: b, ib, jb
    do b = 1, nblk(1)*nblk(2)
      ib = mod(b-1, nblk(1)) + 1
      jb = (b-1)/nblk(1) + 1
      select case ( kind )
      case ( 1 ); act(b) = .true.
      case ( 2 ); act(b) = .false.
      case ( 3 ); act(b) = b == (nblk(1)*nblk(2) + 1)/2
      case ( 4 ); act(b) = mod(ib + jb, 2) == 0
      case ( 5 ); act(b) = ib == 1 .or. ib == nblk(1) .or. jb == 1 .or. jb == nblk(2)
      case default; act(b) = mod(7*b, 5) < 3
      end select
    enddo
  end function pattern


  subroutine halo_one_rank ()
    type(morph_exchange_serial_t) :: ex
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:), v0(:,:), out(:,:)
    logical :: halo, cnt
    integer(I4) :: k

    call fields ( 30, 25, phi0, v0 )
    out = SENTINEL + 0.0_R8 * phi0
    call band_advance ( ex, 5, [6,5], [ (0, k = 1, 30) ], pattern(6, [6,5]), phi0, v0, 0, out, halo, cnt, st )
    call check ( st%ok() .and. halo .and. cnt, &
                 'ghosts, 1 rank, scattered band: the active neighbours'' nodes, extrapolated elsewhere' )
  end subroutine halo_one_rank


  !> One case of the band on 1 rank: NSTEP steps on the active set of kind,
  !> against the oracle; kind 1 also against the whole grid's step.
  subroutine one_rank_case ( nbn, nblk, kind, name )
    integer(I4),      intent(in) :: nbn, nblk(2), kind
    character(len=*), intent(in) :: name
    type(morph_exchange_serial_t) :: ex
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:), v0(:,:), out(:,:), ref(:,:), full(:,:)
    logical, allocatable  :: act(:)
    logical :: halo, cnt, whole
    integer(I4) :: k

    call fields ( nbn*nblk(1), nbn*nblk(2), phi0, v0 )
    act = pattern ( kind, nblk )
    ref = masked ( nbn, nblk, act, phi0 )
    call oracle_steps ( nbn, nblk, act, v0, ref, NSTEP )
    whole = .true.
    if ( kind == 1 ) then
      full = phi0
      do k = 1, NSTEP
        call Morph_LS_HJ_Step_Full ( full, v0, H, DT, st )
      enddo
      whole = same( full, ref )
    endif
    out = phi0
    out = qnan()
    call band_advance ( ex, nbn, nblk, [ (0, k = 1, nblk(1)*nblk(2)) ], act, phi0, v0, NSTEP, out, halo, cnt, st )
    call check ( st%ok() .and. halo .and. same( out, ref ) .and. whole, name )
  end subroutine one_rank_case


  subroutine steps_one_rank ()
    call one_rank_case ( 8, [5,4], 1, 'band, 1 rank, every block active: the oracle''s bits, the whole grid''s too' )
    call one_rank_case ( 5, [6,5], 2, 'band, 1 rank, no block active: nothing moves' )
    call one_rank_case ( 5, [5,5], 3, 'band, 1 rank, one block: the oracle''s bits' )
    call one_rank_case ( 3, [6,5], 4, 'band, 1 rank, a checkerboard of blocks of 3: the oracle''s bits' )
    call one_rank_case ( 8, [5,4], 5, 'band, 1 rank, a ring of blocks of 8: the oracle''s bits' )
    call one_rank_case ( 5, [6,5], 6, 'band, 1 rank, scattered blocks of 5: the oracle''s bits' )
  end subroutine steps_one_rank


  subroutine front_one_rank ()
    type(morph_exchange_serial_t) :: ex
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:), v0(:,:), out(:,:), ref(:,:)
    logical :: filled, moved
    integer(I4) :: k

    call fields ( 40, 35, phi0, v0 )
    ref = oracle_front ( 5, [8,7], v0 )
    out = phi0
    out = qnan()
    call front_advance ( ex, 5, [8,7], [ (0, k = 1, 56) ], 1, v0, out, filled, st )
    moved = count( .not. ieee_is_nan(ref) ) > 0 .and. &
            any( wanted_k(5, [8,7], 0) .and. .not. wanted_k(5, [8,7], KFRONT-1) ) .and. &
            any( wanted_k(5, [8,7], KFRONT-1) .and. .not. wanted_k(5, [8,7], 0) )
    call check ( st%ok() .and. filled .and. moved .and. same( out, ref ), &
                 'band following a circle, 1 rank: blocks added and dropped, the oracle''s bits' )
  end subroutine front_one_rank


  ! ===========================================================================
  !> n ranks: each asks for a different subset; the set is their union, the
  !> same on every rank, with the schedule of Init with that set; from(s) is
  !> the old slot on each rank.
  subroutine update_ranks ()
    call update_ranks_n ( 2 )
    call update_ranks_n ( 3 )
    call update_ranks_n ( 4 )
  end subroutine update_ranks


  subroutine update_ranks_n ( n )
    integer(I4), intent(in) :: n
    type(morph_sim_world_t), target :: world
    type(morph_status_t) :: st
    integer(I4) :: owner(20), codes(0:n-1), k
    logical :: sameset(0:n-1), sched(0:n-1), froms(0:n-1), empty(0:n-1), union(20)
    character(len=8) :: tag

    do k = 1, 20
      owner(k) = mod(k + k/4, n)       ! every rank owns blocks
      union(k) = mod(7*k, 5) /= 0      ! blocks 5, 10, 15, 20 dropped, inside the ranks' lists
    enddo
    call world%init ( n, st )
    codes = -1
    !$omp parallel num_threads(n)
    block
      type(morph_exchange_sim_t) :: ex
      type(morph_ls_grid_t) :: g, gi
      type(morph_status_t) :: s
      integer(I4), allocatable :: from(:), oldblk(:)
      logical :: want(20), none(20)
      integer(I4) :: t, b, q, expect
      t = 0
      !$ t = omp_get_thread_num()
      ex%me = t
      ex%world => world
      do b = 1, 20
        want(b) = union(b) .and. ( mod(b, n) == t .or. mod(b+1, n) == t )    ! two ranks ask for each
      enddo
      call Morph_LS_Grid_Init ( g, 4, [5,4], H, owner, ex, s )
      oldblk = g%blk
      call Morph_LS_Grid_Update ( g, ex, want, from, s )
      codes(t) = s%code
      sameset(t) = all( g%active .eqv. union )
      call Morph_LS_Grid_Init ( gi, 4, [5,4], H, owner, ex, s, union )
      sched(t) = s%ok() .and. schedule_same ( g, gi )
      froms(t) = size(from) == size(g%blk)
      if ( froms(t) ) then
        do q = 1, size(g%blk)
          expect = 0
          do b = 1, size(oldblk)
            if ( oldblk(b) == g%blk(q) ) expect = b
          enddo
          if ( from(q) /= expect ) froms(t) = .false.
        enddo
      endif
      none = .false.
      call Morph_LS_Grid_Update ( g, ex, none, from, s )
      empty(t) = s%ok() .and. size(g%blk) == 0 .and. .not. any( g%active ) .and. size(from) == 0
    end block
    !$omp end parallel
    write(tag,'(I0,A)') n, ' ranks'
    call check ( st%ok() .and. all( codes == MORPH_OK ) .and. all( sameset ) .and. all( sched ), &
                 'update, '//trim(tag)//': the union of the requests on every rank, the schedule of Init with it' )
    call check ( all( froms ), 'update, '//trim(tag)//': from(s) is the old slot on every rank' )
    call check ( all( empty ), 'update, '//trim(tag)//': nothing asked anywhere, no block on any rank' )
  end subroutine update_ranks_n


  !> A refusal on one rank is taken by every rank, the grids untouched.
  subroutine update_refusals ()
    type(morph_sim_world_t), target :: world
    type(morph_status_t) :: st
    integer(I4) :: codes(0:1,4), owner(12), k
    logical :: kept(0:1,4), named(0:1,4)

    do k = 1, 12
      owner(k) = mod(k, 2)
    enddo
    call world%init ( 2, st )
    codes = -1
    kept  = .false.
    named = .false.
    !$omp parallel num_threads(2)
    block
      type(morph_exchange_sim_t) :: ex, oth
      type(morph_ls_grid_t) :: g, keep, og, gnot
      type(morph_status_t) :: s
      integer(I4), allocatable :: from(:)
      logical :: want(12)
      integer(I4) :: t, q
      t = 0
      !$ t = omp_get_thread_num()
      ex%me = t
      ex%world => world
      want = .true.
      want(5) = .false.
      call Morph_LS_Grid_Init ( g, 4, [4,3], H, owner, ex, s, [ (mod(q,3) /= 0, q = 1, 12) ] )
      keep = g
      ! 1: rank 1 passes 11 requests
      if ( t == 1 ) then
        call Morph_LS_Grid_Update ( g, ex, want(1:11), from, s )
      else
        call Morph_LS_Grid_Update ( g, ex, want, from, s )
      endif
      codes(t,1) = s%code
      kept(t,1)  = schedule_same ( g, keep ) .and. all( from == [ (q, q = 1, size(g%blk)) ] )
      named(t,1) = merge( index(s%message, '11 requests for 12 blocks') > 0, &
                          index(s%message, 'refused on rank 1 (status') > 0, t == 1 )
      ! 2: rank 0 updates a grid made for rank 1
      oth%me = 1
      oth%world => world
      if ( t == 0 ) then
        call Morph_LS_Grid_Init ( og, 4, [4,3], H, owner, oth, s )
        keep = og
        call Morph_LS_Grid_Update ( og, ex, want, from, s )
        kept(t,2) = schedule_same ( og, keep )
      else
        call Morph_LS_Grid_Init ( g, 4, [4,3], H, owner, ex, s )
        keep = g
        call Morph_LS_Grid_Update ( g, ex, want, from, s )
        kept(t,2) = schedule_same ( g, keep )
      endif
      codes(t,2) = s%code
      named(t,2) = merge( index(s%message, 'the grid rank 1 of 2') > 0, &
                          index(s%message, 'refused on rank 0 (status') > 0, t == 0 )
      ! 3: grids of different sizes, the requests within the smaller one (so that only the check of the
      ! sizes stops them)
      if ( t == 0 ) then
        call Morph_LS_Grid_Init ( g, 4, [4,3], H, owner, ex, s )
        keep = g
        call Morph_LS_Grid_Update ( g, ex, [ want(1:9), .false., .false., .false. ], from, s )
      else
        call Morph_LS_Grid_Init ( g, 4, [3,3], H, owner(1:9), ex, s )
        keep = g
        call Morph_LS_Grid_Update ( g, ex, want(1:9), from, s )
      endif
      codes(t,3) = s%code
      kept(t,3)  = schedule_same ( g, keep )
      named(t,3) = index(s%message, 'grids of different sizes') > 0
      ! 4: rank 0's grid not initialised
      if ( t == 0 ) then
        call Morph_LS_Grid_Update ( gnot, ex, want, from, s )
        kept(t,4) = size(from) == 0
      else
        call Morph_LS_Grid_Init ( g, 4, [4,3], H, owner, ex, s )
        keep = g
        call Morph_LS_Grid_Update ( g, ex, want, from, s )
        kept(t,4) = schedule_same ( g, keep )
      endif
      codes(t,4) = s%code
      named(t,4) = merge( index(s%message, 'not initialised') > 0, &
                          index(s%message, 'refused on rank 0 (status') > 0, t == 0 )
    end block
    !$omp end parallel
    call check ( st%ok() .and. all( codes(:,1) == MORPH_ERR_SHAPE ) .and. all( kept(:,1) ) .and. all( named(:,1) ), &
                 'update refused on rank 1 (11 requests for 12 blocks): by both, grids untouched, from(s) = s' )
    call check ( all( codes(:,2) == MORPH_ERR_EXCHANGE ) .and. all( kept(:,2) ) .and. all( named(:,2) ), &
                 'update refused, rank 0 with rank 1''s grid: by both, the grid untouched' )
    call check ( all( codes(:,3) == MORPH_ERR_EXCHANGE ) .and. all( kept(:,3) ) .and. all( named(:,3) ), &
                 'update refused, grids of 12 and 9 blocks on the two ranks: by both, untouched' )
    call check ( all( codes(:,4) == MORPH_ERR_SHAPE ) .and. all( kept(:,4) ) .and. all( named(:,4) ), &
                 'update refused, rank 0''s grid not initialised: by both' )
  end subroutine update_refusals


  ! ===========================================================================
  !> n simulated ranks, a band of fixed active set: the oracle's ghosts and
  !> strips on every rank, and the oracle's bits after NSTEP steps.
  subroutine ranks_case ( n, nbn, nblk, owner, kind, nsteps, name )
    integer(I4),      intent(in) :: n, nbn, nblk(2), owner(:), kind, nsteps
    character(len=*), intent(in) :: name
    type(morph_sim_world_t), target :: world
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:), v0(:,:), ref(:,:), out(:,:)
    logical, allocatable  :: act(:)
    integer(I4) :: codes(0:n-1)
    logical :: halos(0:n-1), cnts(0:n-1)

    call fields ( nbn*nblk(1), nbn*nblk(2), phi0, v0 )
    act = pattern ( kind, nblk )
    ref = masked ( nbn, nblk, act, phi0 )
    call oracle_steps ( nbn, nblk, act, v0, ref, nsteps )
    out = phi0
    out = qnan()
    call world%init ( n, st )
    codes = -1
    halos = .false.
    cnts  = .false.
    !$omp parallel num_threads(n)
    block
      type(morph_exchange_sim_t) :: ex
      type(morph_status_t) :: s
      integer(I4) :: t
      logical :: hk, ck
      t = 0
      !$ t = omp_get_thread_num()
      ex%me = t
      ex%world => world
      call band_advance ( ex, nbn, nblk, owner, act, phi0, v0, nsteps, out, hk, ck, s )
      codes(t) = s%code
      halos(t) = hk
      cnts(t)  = ck
    end block
    !$omp end parallel
    call check ( st%ok() .and. all( codes == MORPH_OK ) .and. all( halos ) .and. all( cnts ) .and. &
                 same( out, ref ), name )
  end subroutine ranks_case


  subroutine halo_ranks ()
    integer(I4) :: own(30), k
    do k = 1, 30
      own(k) = mod(5*k, 3)
    enddo
    call ranks_case ( 3, 5, [6,5], own, 6, 0, &
                      'ghosts, 3 ranks, scattered band: the oracle''s on every rank, strips of active pairs only' )
  end subroutine halo_ranks


  subroutine steps_ranks ()
    integer(I4) :: own20(20), own30(30), k, ib, jb
    do k = 1, 20
      ib = mod(k-1, 5) + 1
      jb = (k-1)/5 + 1
      own20(k) = mod(ib + jb, 2)
    enddo
    call ranks_case ( 2, 8, [5,4], own20, 6, NSTEP, 'band, 2 ranks, checkerboard owners, scattered band' )
    call ranks_case ( 2, 8, [5,4], own20, 4, NSTEP, 'band, 2 ranks, the band a checkerboard too: no strip' )
    do k = 1, 30
      own30(k) = mod(k-1, 3)
    enddo
    call ranks_case ( 3, 5, [6,5], own30, 5, NSTEP, 'band, 3 ranks, blocks dealt in turn, a ring' )
    do k = 1, 30
      own30(k) = 2*mod(k-1, 2)
    enddo
    call ranks_case ( 3, 5, [6,5], own30, 6, NSTEP, 'band, 3 ranks, rank 1 without blocks, scattered band' )
    do k = 1, 20
      own20(k) = (k-1)/5
    enddo
    call ranks_case ( 4, 8, [5,4], own20, 3, NSTEP, 'band, 4 ranks, one block active on one rank' )
    call ranks_case ( 4, 8, [5,4], own20, 2, NSTEP, 'band, 4 ranks, no block active' )
    do k = 1, 30
      own30(k) = mod(5*k, 4)
    enddo
    call ranks_case ( 4, 3, [6,5], own30, 6, NSTEP, 'band, 4 ranks, blocks of 3, scattered owners and band' )
  end subroutine steps_ranks


  subroutine front_ranks ()
    integer(I4) :: own(56), k
    do k = 1, 56
      own(k) = mod(7*k, 3)
    enddo
    call front_case ( 3, own, 1, 'band following a circle, 3 ranks, each asking for its own and its neighbours' )
    do k = 1, 56
      own(k) = (k-1)/14
    enddo
    call front_case ( 4, own, 2, 'band following a circle, 4 ranks, requests dealt whoever owns the block' )
  end subroutine front_ranks


  subroutine front_case ( n, owner, mode, name )
    integer(I4),      intent(in) :: n, owner(:), mode
    character(len=*), intent(in) :: name
    type(morph_sim_world_t), target :: world
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:), v0(:,:), out(:,:), ref(:,:)
    integer(I4) :: codes(0:n-1)
    logical :: fills(0:n-1)

    call fields ( 40, 35, phi0, v0 )
    ref = oracle_front ( 5, [8,7], v0 )
    out = phi0
    out = qnan()
    call world%init ( n, st )
    codes = -1
    fills = .false.
    !$omp parallel num_threads(n)
    block
      type(morph_exchange_sim_t) :: ex
      type(morph_status_t) :: s
      integer(I4) :: t
      logical :: f
      t = 0
      !$ t = omp_get_thread_num()
      ex%me = t
      ex%world => world
      call front_advance ( ex, 5, [8,7], owner, mode, v0, out, f, s )
      codes(t) = s%code
      fills(t) = f
    end block
    !$omp end parallel
    call check ( st%ok() .and. all( codes == MORPH_OK ) .and. all( fills ) .and. same( out, ref ), name )
  end subroutine front_case

end program test_morph_ls_band
