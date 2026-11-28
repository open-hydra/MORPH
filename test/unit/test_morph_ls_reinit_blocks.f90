!> Unit test: the reinitialisation and the eikonal error of Morph_LS_Reinit on
!> the grid of blocks of Morph_LS_Grid.
!>
!> Every comparison of fields is bit for bit (the 64 bits of each value):
!>   * every block active, blocks of 3, 5 and 8 nodes: Morph_LS_Reinit_Grid
!>     gives the bits of Morph_LS_Reinit_Full on the whole grid (4 steps with
!>     the correction, 3 without; 0 steps leave phi untouched);
!>   * a band (one block alone, a checkerboard, a frame, scattered blocks, a
!>     ring): the bits of an oracle written here on the whole grid, which
!>     gives each active block the ghosts the band reads (the whole grid's
!>     nodes beyond a side whose block is active, else the linear
!>     extrapolation of the block's own nodes) and calls the public kernels
!>     of Morph_LS_Reinit; the inactive blocks untouched;
!>   * the reach: on a ring of blocks, from the same phi, every node farther
!>     than 4 iters nodes (L1) from all inactive nodes has the bits of the
!>     whole grid (the stencils: 3 nodes for a step, 1 for the correction);
!>   * the measure: on one block covering the grid, the bits of
!>     Morph_LS_Eikonal_Error; on blocks, the bits of the oracle's sum of
!>     Morph_LS_Eikonal_Terms in block order, and within 2 (n+1) u E of the
!>     whole grid's (u = 2^-53; two recursive sums of the same n non-negative
!>     terms, each within (n-1) u, and two divisions, each within u: 2 n u E
!>     at first order, the extra 2 u E covering the higher orders while
!>     n < 6.7e7); +0 with no node in the band or no active block;
!>   * 2 to 4 simulated ranks, blocks of 3, 5 and 8, one rank with no active
!>     block: the bits of one rank, E the same bits on every rank;
!>     Morph_LS_Reinit_If_Needed_Grid below, at and above the threshold;
!>   * refusals, phi untouched: a grid whose Init was refused (on one rank
!>     and on one of two), shapes, iters < 0, beta < 0 or NaN, shape
!>     and iters both wrong (the shape reported), an exchange that is not
!>     the grid's, one rank's refusal taken by every rank (two refusing: rank
!>     0 takes the first one's code), iters, beta or a threshold that differ
!>     between ranks; an exchange that fails at each of the collectives of
!>     each kind that the calls make (allgather, swap, ordered sum), and one
!>     call past the last of each kind that succeeds (the counts are those of
!>     the calls, also at 0 steps and without the correction).
!> Built without OpenMP, the simulated ranks are skipped (said so).
module test_reinit_failing_exchange
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m
  use Morph_Exchange,  only: morph_exchange_serial_t
  implicit none
  private

  !> One rank; the call of a kind whose count reaches fail(kind) fails (kind
  !> 1 allgather, 2 swap, 3 ordered sum; 0 for none). calls(kind): the calls.
  type, public, extends(morph_exchange_serial_t) :: failing_exchange_t
    integer(I4) :: fail(3)  = 0
    integer(I4) :: calls(3) = 0
  contains
    procedure :: allgather   => f_allgather
    procedure :: swap        => f_swap
    procedure :: ordered_sum => f_ordered_sum
  end type failing_exchange_t

contains

  subroutine f_allgather ( self, mine, all, counts, status )
    class(failing_exchange_t), intent(inout) :: self
    real(R8),                  intent(in)    :: mine(:)
    real(R8),                  intent(out)   :: all(:)
    integer(I4),               intent(in)    :: counts(0:)
    type(morph_status_t),      intent(out)   :: status
    self%calls(1) = self%calls(1) + 1
    call self%morph_exchange_serial_t%allgather ( mine, all, counts, status )
    if ( self%calls(1) == self%fail(1) ) call status%fail ( MORPH_ERR_EXCHANGE, 'failing exchange: allgather' )
  end subroutine f_allgather

  subroutine f_swap ( self, scount, sbuf, rcount, rbuf, status )
    class(failing_exchange_t), intent(inout) :: self
    integer(I4),               intent(in)    :: scount(0:)
    real(R8),                  intent(in)    :: sbuf(:)
    integer(I4),               intent(in)    :: rcount(0:)
    real(R8),                  intent(out)   :: rbuf(:)
    type(morph_status_t),      intent(out)   :: status
    self%calls(2) = self%calls(2) + 1
    call self%morph_exchange_serial_t%swap ( scount, sbuf, rcount, rbuf, status )
    if ( self%calls(2) == self%fail(2) ) call status%fail ( MORPH_ERR_EXCHANGE, 'failing exchange: swap' )
  end subroutine f_swap

  subroutine f_ordered_sum ( self, nblk, owner, partial, total, status )
    class(failing_exchange_t), intent(inout) :: self
    integer(I4),               intent(in)    :: nblk
    integer(I4),               intent(in)    :: owner(:)
    real(R8),                  intent(in)    :: partial(:)
    real(R8),                  intent(out)   :: total
    type(morph_status_t),      intent(out)   :: status
    self%calls(3) = self%calls(3) + 1
    call self%morph_exchange_serial_t%ordered_sum ( nblk, owner, partial, total, status )
    if ( self%calls(3) == self%fail(3) ) call status%fail ( MORPH_ERR_EXCHANGE, 'failing exchange: ordered sum' )
  end subroutine f_ordered_sum

end module test_reinit_failing_exchange


program test_morph_ls_reinit_blocks
  use iso_fortran_env, only: I4 => int32, I8 => int64, R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  !$ use omp_lib, only: omp_get_thread_num
  use Morph_Types_m
  use Morph_Exchange
  use Morph_Exchange_Sim
  use Morph_LS_HJ,     only: MORPH_LS_NG, Morph_LS_Ghosts_Extrap
  use Morph_LS_Grid
  use Morph_LS_Reinit
  use test_reinit_failing_exchange

  implicit none

  integer(I4), parameter :: NG = MORPH_LS_NG
  real(R8),    parameter :: H = 1.0_R8 / 32.0_R8
  real(R8),    parameter :: BAND = 2.0_R8 * H
  real(R8),    parameter :: BETA = 0.5_R8
  real(R8),    parameter :: SENTINEL = -7.0_R8               ! no field value of the test
  real(R8),    parameter :: U = 0.5_R8 * epsilon(1.0_R8)     ! 2^-53
  integer(I4), parameter :: NBNS(3) = [3, 5, 8]              ! the block sizes

  integer(I4) :: nfail
  logical     :: omp

  nfail = 0
  omp = .false.
  !$ omp = .true.
  write(*,'(A)') '== Morph level-set reinitialisation on the grid of blocks =='

  call measure_one_rank ()
  call reinit_every_block ()
  call reinit_band ()
  call reach ()
  call if_needed_one_rank ()
  call refusals_one_rank ()
  call failing_exchange ()
  if ( omp ) then
    call ranks ()
    call if_needed_ranks ()
    call refusals_ranks ()
    call two_refusals_ranks ()
    call uninitialised_rank ()
  else
    write(*,'(A)') '   SKIP the simulated ranks: built without OpenMP'
  endif

  if ( nfail > 0 ) then
    write(*,'(I0,A)') nfail, ' REINITIALISATION-ON-BLOCKS CHECK(S) FAILED'
    error stop 1
  endif
  write(*,'(A)') 'RESULT: all reinitialisation-on-blocks checks PASS'

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


  logical function same1 ( a, b )
    real(R8), intent(in) :: a, b
    same1 = transfer(a, 0_I8) == transfer(b, 0_I8)
  end function same1


  !> phi on an nx x ny grid: the distance to a circle away from the grid's
  !> centre times 1.5 + 0.5 x/lx, so |grad phi| runs from about 1.2 to 1.8
  !> and the reinitialisation works; its zero level crosses blocks.
  function field ( nx, ny ) result ( phi )
    integer(I4), intent(in) :: nx, ny
    real(R8), allocatable   :: phi(:,:)
    real(R8)    :: x, y, lx, ly
    integer(I4) :: i, j
    allocate( phi(nx,ny) )
    lx = H*real(nx-1, R8)
    ly = H*real(ny-1, R8)
    do j = 1, ny
      do i = 1, nx
        x = H*real(i-1, R8) - 0.53_R8*lx
        y = H*real(j-1, R8) - 0.47_R8*ly
        phi(i,j) = (0.31_R8*min(lx, ly) - hypot(x, y)) * (1.5_R8 + 0.5_R8*x/lx)
      enddo
    enddo
  end function field


  !> The active blocks of a kind: 1 all, 2 none, 3 the middle one, 4 a
  !> checkerboard, 5 the frame, 6 scattered, 7 a ring along the zero level
  !> of field (the blocks whose centre lies within 2.4 block sides of the
  !> circle, the one round its centre left out).
  function pattern ( kind, nbn, nblk ) result ( act )
    integer(I4), intent(in) :: kind, nbn, nblk(2)
    logical :: act(nblk(1)*nblk(2))
    integer(I4) :: b, ib, jb
    real(R8) :: lx, ly, x, y
    lx = H*real(nbn*nblk(1)-1, R8)
    ly = H*real(nbn*nblk(2)-1, R8)
    do b = 1, nblk(1)*nblk(2)
      ib = mod(b-1, nblk(1)) + 1
      jb = (b-1)/nblk(1) + 1
      select case ( kind )
      case ( 1 ); act(b) = .true.
      case ( 2 ); act(b) = .false.
      case ( 3 ); act(b) = b == (nblk(1)*nblk(2) + 1)/2
      case ( 4 ); act(b) = mod(ib + jb, 2) == 0
      case ( 5 ); act(b) = ib == 1 .or. ib == nblk(1) .or. jb == 1 .or. jb == nblk(2)
      case ( 6 ); act(b) = mod(7*b, 5) < 3
      case default
        x = H*(real((ib-1)*nbn, R8) + 0.5_R8*real(nbn-1, R8)) - 0.53_R8*lx
        y = H*(real((jb-1)*nbn, R8) + 0.5_R8*real(nbn-1, R8)) - 0.47_R8*ly
        act(b) = abs( hypot(x, y) - 0.31_R8*min(lx, ly) ) < 2.4_R8*H*real(nbn, R8)
      end select
    enddo
  end function pattern


  !> Block b on rank mod(m*b, n): owners for n ranks.
  function owners ( n, nb, m ) result ( own )
    integer(I4), intent(in) :: n, nb, m
    integer(I4) :: own(nb)
    integer(I4) :: b
    do b = 1, nb
      own(b) = mod(m*b, n)
    enddo
  end function owners


  !> Owners on n ranks such that rank n-1 owns only inactive blocks.
  function owners_no_active ( n, act ) result ( own )
    integer(I4), intent(in) :: n
    logical,     intent(in) :: act(:)
    integer(I4) :: own(size(act))
    integer(I4) :: b
    do b = 1, size(act)
      own(b) = n-1
      if ( act(b) ) own(b) = mod(b, n-1)
    enddo
  end function owners_no_active


  !> The edge flags (W, E, S, N) of block b: the side faces the edge of the
  !> grid or an inactive block.
  function edge_of ( nblk, act, b ) result ( e )
    integer(I4), intent(in) :: nblk(2), b
    logical,     intent(in) :: act(:)
    logical :: e(4)
    integer(I4) :: ib, jb
    ib = mod(b-1, nblk(1)) + 1
    jb = (b-1)/nblk(1) + 1
    e(1) = ib == 1
    if ( .not. e(1) ) e(1) = .not. act(b-1)
    e(2) = ib == nblk(1)
    if ( .not. e(2) ) e(2) = .not. act(b+1)
    e(3) = jb == 1
    if ( .not. e(3) ) e(3) = .not. act(b-nblk(1))
    e(4) = jb == nblk(2)
    if ( .not. e(4) ) e(4) = .not. act(b+nblk(1))
  end function edge_of


  !> Block b of the whole grid f with the ghosts the band reads: the whole
  !> grid's nodes beyond a side whose block is active, else the linear
  !> extrapolation of the block's own nodes; the corners SENTINEL. Written
  !> here again from the layout of Morph_LS_Grid.
  function with_ghosts ( nbn, nblk, act, f, b ) result ( w )
    integer(I4), intent(in) :: nbn, nblk(2), b
    logical,     intent(in) :: act(:)
    real(R8),    intent(in) :: f(:,:)
    real(R8), allocatable   :: w(:,:)
    logical     :: e(4)
    integer(I4) :: i0, j0
    allocate( w(1-NG:nbn+NG, 1-NG:nbn+NG) )
    w = SENTINEL
    i0 = mod(b-1, nblk(1))*nbn
    j0 = ((b-1)/nblk(1))*nbn
    w(1:nbn,1:nbn) = f(i0+1:i0+nbn, j0+1:j0+nbn)
    e = edge_of ( nblk, act, b )
    if ( e(1) ) then
      call Morph_LS_Ghosts_Extrap ( w, 1, 1 )
    else
      w(1-NG:0, 1:nbn) = f(i0-NG+1:i0, j0+1:j0+nbn)
    endif
    if ( e(2) ) then
      call Morph_LS_Ghosts_Extrap ( w, 1, 2 )
    else
      w(nbn+1:nbn+NG, 1:nbn) = f(i0+nbn+1:i0+nbn+NG, j0+1:j0+nbn)
    endif
    if ( e(3) ) then
      call Morph_LS_Ghosts_Extrap ( w, 2, 1 )
    else
      w(1:nbn, 1-NG:0) = f(i0+1:i0+nbn, j0-NG+1:j0)
    endif
    if ( e(4) ) then
      call Morph_LS_Ghosts_Extrap ( w, 2, 2 )
    else
      w(1:nbn, nbn+1:nbn+NG) = f(i0+1:i0+nbn, j0+nbn+1:j0+nbn+NG)
    endif
  end function with_ghosts


  !> The oracle of the reinitialisation on a band, on the whole grid u: the
  !> public kernels on each active block with the ghosts of with_ghosts, the
  !> exchanges replaced by reading the whole grid between the kernels.
  subroutine oracle_reinit ( nbn, nblk, act, u, iters, beta )
    integer(I4), intent(in)    :: nbn, nblk(2), iters
    logical,     intent(in)    :: act(:)
    real(R8),    intent(inout) :: u(:,:)
    real(R8),    intent(in)    :: beta
    real(R8), allocatable :: p0(:,:), st(:,:), s(:,:,:), dd(:,:,:), stb(:,:)
    logical,  allocatable :: pos(:,:,:), cut(:,:,:)
    integer(I4) :: b, it, i0, j0, nb
    nb = nblk(1)*nblk(2)
    p0 = u
    st = u
    allocate( s(nbn,nbn,nb), dd(nbn,nbn,nb), pos(nbn,nbn,nb), cut(nbn,nbn,nb), stb(nbn,nbn) )
    do b = 1, nb
      if ( .not. act(b) ) cycle
      call Morph_LS_Reinit_Setup ( with_ghosts(nbn, nblk, act, p0, b), H, edge_of(nblk, act, b), &
                                   s(:,:,b), pos(:,:,b), cut(:,:,b), dd(:,:,b) )
    enddo
    do it = 1, iters
      do b = 1, nb
        if ( .not. act(b) ) cycle
        i0 = mod(b-1, nblk(1))*nbn
        j0 = ((b-1)/nblk(1))*nbn
        call Morph_LS_Reinit_Advance ( with_ghosts(nbn, nblk, act, u, b), pos(:,:,b), s(:,:,b), cut(:,:,b), &
                                       dd(:,:,b), H, stb )
        st(i0+1:i0+nbn, j0+1:j0+nbn) = stb
      enddo
      do b = 1, nb
        if ( .not. act(b) ) cycle
        i0 = mod(b-1, nblk(1))*nbn
        j0 = ((b-1)/nblk(1))*nbn
        if ( beta > 0.0_R8 ) then
          call Morph_LS_Reinit_Correct ( with_ghosts(nbn, nblk, act, st, b), with_ghosts(nbn, nblk, act, p0, b), &
                                         cut(:,:,b), H, beta, edge_of(nblk, act, b), stb )
          u(i0+1:i0+nbn, j0+1:j0+nbn) = stb
        else
          u(i0+1:i0+nbn, j0+1:j0+nbn) = st(i0+1:i0+nbn, j0+1:j0+nbn)
        endif
      enddo
    enddo
  end subroutine oracle_reinit


  !> The oracle of the measure: Morph_LS_Eikonal_Terms of each active block
  !> (with_ghosts), summed in block order, the counts likewise.
  real(R8) function oracle_measure ( nbn, nblk, act, f ) result ( e )
    integer(I4), intent(in) :: nbn, nblk(2)
    logical,     intent(in) :: act(:)
    real(R8),    intent(in) :: f(:,:)
    real(R8)    :: s, tot, cnt
    integer(I4) :: b, n
    tot = 0.0_R8
    cnt = 0.0_R8
    do b = 1, nblk(1)*nblk(2)
      if ( .not. act(b) ) cycle
      call Morph_LS_Eikonal_Terms ( with_ghosts(nbn, nblk, act, f, b), H, BAND, edge_of(nblk, act, b), s, n )
      tot = tot + s
      cnt = cnt + real(n, R8)
    enddo
    e = 0.0_R8
    if ( cnt > 0.0_R8 ) e = tot / cnt
  end function oracle_measure


  !> The whole grid with SENTINEL on the nodes of the inactive blocks.
  function masked ( nbn, nblk, act, full ) result ( v )
    integer(I4), intent(in) :: nbn, nblk(2)
    logical,     intent(in) :: act(:)
    real(R8),    intent(in) :: full(:,:)
    real(R8), allocatable   :: v(:,:)
    integer(I4) :: b, i0, j0
    v = full
    do b = 1, nblk(1)*nblk(2)
      if ( act(b) ) cycle
      i0 = mod(b-1, nblk(1))*nbn
      j0 = ((b-1)/nblk(1))*nbn
      v(i0+1:i0+nbn, j0+1:j0+nbn) = SENTINEL
    enddo
  end function masked


  !> One rank's share of a band: Init with owner and act, scatter phi0, E
  !> (Morph_LS_Eikonal_Error_Grid), iters steps (Morph_LS_Reinit_Grid) and
  !> gather into out (this rank's blocks only). code: the first failure.
  subroutine band_run ( exch, nbn, nblk, owner, act, phi0, iters, beta, out, e, code )
    class(morph_exchange_t), intent(inout) :: exch
    integer(I4),             intent(in)    :: nbn, nblk(2), owner(:), iters
    logical,                 intent(in)    :: act(:)
    real(R8),                intent(in)    :: phi0(:,:), beta
    real(R8),                intent(inout) :: out(:,:)
    real(R8),                intent(out)   :: e
    integer(I4),             intent(out)   :: code
    type(morph_ls_grid_t) :: g
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi(:,:,:)
    e = -1.0_R8
    call Morph_LS_Grid_Init ( g, nbn, nblk, H, owner, exch, st, act )
    code = st%code
    if ( .not. st%ok() ) return
    allocate( phi(nbn,nbn,size(g%blk)) )
    call Morph_LS_Grid_Scatter ( g, phi0, phi, st )
    call Morph_LS_Eikonal_Error_Grid ( g, exch, phi, BAND, e, st )
    code = st%code
    if ( .not. st%ok() ) return
    call Morph_LS_Reinit_Grid ( g, exch, phi, iters, beta, st )
    code = st%code
    if ( .not. st%ok() ) return
    call Morph_LS_Grid_Gather ( g, phi, out, st )
    code = st%code
  end subroutine band_run


  !> This rank's slots of the whole grid f (Morph_LS_Grid_Scatter).
  function slots_of ( g, f ) result ( p )
    type(morph_ls_grid_t), intent(in) :: g
    real(R8),              intent(in) :: f(:,:)
    real(R8), allocatable :: p(:,:,:)
    type(morph_status_t)  :: st
    allocate( p(g%nbn, g%nbn, size(g%blk)) )
    call Morph_LS_Grid_Scatter ( g, f, p, st )
  end function slots_of


  subroutine measure_one_rank ()
    type(morph_exchange_serial_t) :: ex
    type(morph_ls_grid_t) :: g
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:), phi(:,:,:)
    real(R8)    :: e, ef, eo, band1
    integer(I4) :: k, nbn, nblk(2), n, b, i, j
    logical     :: ok, okb, ok0

    ! one block covering the whole grid: the additions of the whole grid
    phi0 = field ( 23, 23 )
    call Morph_LS_Grid_Init ( g, 23, [1,1], H, [0], ex, st )
    phi = slots_of ( g, phi0 )
    call Morph_LS_Eikonal_Error_Grid ( g, ex, phi, BAND, e, st )
    call check ( st%ok() .and. same1( e, Morph_LS_Eikonal_Error(phi0, H, BAND) ) .and. e > 0.0_R8, &
                 'E on one block covering a 23 x 23 grid: the bits of the whole grid' )

    ! one node in the band, the corner (1,1): E is its term alone, the
    ! gradient one-sided along both axes. Across the linearly extrapolated
    ! ghost the central difference is the same in exact arithmetic; with
    ! this field the two round differently along x there (an exact
    ! emulation: no fused multiply-add enters them), so a gradient that
    ! ignored the edge would show.
    deallocate( phi0 )
    allocate( phi0(8,8) )
    do j = 1, 8
      do i = 1, 8
        phi0(i,j) = (H*real(i-1, R8) - 0.41_R8*H) * (1.5_R8 + 0.5_R8*real(j-1, R8)/7.0_R8)
      enddo
    enddo
    band1 = abs( phi0(1,1) )
    call Morph_LS_Grid_Init ( g, 8, [1,1], H, [0], ex, st )
    phi = slots_of ( g, phi0 )
    call Morph_LS_Eikonal_Error_Grid ( g, ex, phi, band1, e, st )
    call check ( st%ok() .and. count( abs(phi0) <= band1 ) == 1 .and. e > 0.0_R8 .and. &
                 same1( e, Morph_LS_Eikonal_Error(phi0, H, band1) ), &
                 'E of one node in the band, at a corner of the grid: the bits of the whole grid (one-sided)' )

    ! every block active: the oracle's bits; within 2 (n+1) u E of the whole grid
    ok  = .true.
    okb = .true.
    do k = 1, 3
      nbn  = NBNS(k)
      nblk = [ 40/nbn, 37/nbn ]
      phi0 = field ( nbn*nblk(1), nbn*nblk(2) )
      call Morph_LS_Grid_Init ( g, nbn, nblk, H, [ (0, b = 1, nblk(1)*nblk(2)) ], ex, st )
      phi = slots_of ( g, phi0 )
      call Morph_LS_Eikonal_Error_Grid ( g, ex, phi, BAND, e, st )
      eo = oracle_measure ( nbn, nblk, pattern(1, nbn, nblk), phi0 )
      ef = Morph_LS_Eikonal_Error ( phi0, H, BAND )
      n  = count( abs(phi0) <= BAND )
      if ( .not. ( st%ok() .and. same1( e, eo ) ) ) ok = .false.
      if ( .not. ( n > 0 .and. abs(e - ef) <= 2.0_R8*real(n+1, R8)*U*ef ) ) okb = .false.
    enddo
    call check ( ok, 'E on every block of 3, 5 and 8 nodes: the oracle''s sum in block order, bit for bit' )
    call check ( okb, 'E on every block of 3, 5 and 8 nodes: within 2 (n+1) u E of the whole grid' )

    ! no node in the band; then no active block
    call Morph_LS_Eikonal_Error_Grid ( g, ex, phi, -1.0_R8, e, st )
    ok0 = st%ok() .and. same1( e, 0.0_R8 )
    call Morph_LS_Grid_Init ( g, 5, [3,3], H, [ (0, b = 1, 9) ], ex, st, pattern(2, 5, [3,3]) )
    deallocate( phi )
    allocate( phi(5,5,0) )
    call Morph_LS_Eikonal_Error_Grid ( g, ex, phi, BAND, e, st )
    call check ( ok0 .and. st%ok() .and. same1( e, 0.0_R8 ), 'E with no node in the band, then with no active block: +0' )
    call Morph_LS_Reinit_Grid ( g, ex, phi, 4, BETA, st )
    call check ( st%ok(), 'reinitialisation with no active block: done, nothing to do' )
  end subroutine measure_one_rank


  subroutine reinit_every_block ()
    type(morph_exchange_serial_t) :: ex
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:), full(:,:), out(:,:)
    real(R8)    :: e
    integer(I4) :: k, nbn, nblk(2), code, b
    logical     :: ok4, ok3, ok0

    ok4 = .true.
    ok3 = .true.
    ok0 = .true.
    do k = 1, 3
      nbn  = NBNS(k)
      nblk = [ 40/nbn, 37/nbn ]
      phi0 = field ( nbn*nblk(1), nbn*nblk(2) )

      full = phi0
      call Morph_LS_Reinit_Full ( full, H, 4, BETA, st )
      out = phi0
      out = SENTINEL
      call band_run ( ex, nbn, nblk, [ (0, b = 1, nblk(1)*nblk(2)) ], pattern(1, nbn, nblk), phi0, 4, BETA, &
                      out, e, code )
      if ( .not. ( st%ok() .and. code == MORPH_OK .and. same( out, full ) .and. .not. same( out, phi0 ) ) ) &
        ok4 = .false.

      full = phi0
      call Morph_LS_Reinit_Full ( full, H, 3, 0.0_R8, st )
      out = SENTINEL
      call band_run ( ex, nbn, nblk, [ (0, b = 1, nblk(1)*nblk(2)) ], pattern(1, nbn, nblk), phi0, 3, 0.0_R8, &
                      out, e, code )
      if ( .not. ( st%ok() .and. code == MORPH_OK .and. same( out, full ) .and. .not. same( out, phi0 ) ) ) &
        ok3 = .false.

      out = SENTINEL
      call band_run ( ex, nbn, nblk, [ (0, b = 1, nblk(1)*nblk(2)) ], pattern(1, nbn, nblk), phi0, 0, BETA, &
                      out, e, code )
      if ( .not. ( code == MORPH_OK .and. same( out, phi0 ) ) ) ok0 = .false.
    enddo
    call check ( ok4, 'every block of 3, 5 and 8 nodes, 4 steps with the correction: the bits of the whole grid' )
    call check ( ok3, 'every block of 3, 5 and 8 nodes, 3 steps without it: the bits of the whole grid' )
    call check ( ok0, 'every block of 3, 5 and 8 nodes, 0 steps: phi untouched, bit for bit' )
  end subroutine reinit_every_block


  !> One rank, a band of kind: out and E of band_run against the oracles.
  subroutine one_rank_band ( nbn, nblk, kind, iters, beta, name )
    integer(I4),      intent(in) :: nbn, nblk(2), kind, iters
    real(R8),         intent(in) :: beta
    character(len=*), intent(in) :: name
    type(morph_exchange_serial_t) :: ex
    real(R8), allocatable :: phi0(:,:), ref(:,:), out(:,:)
    logical,  allocatable :: act(:)
    real(R8)    :: e, eo
    integer(I4) :: code, b
    phi0 = field ( nbn*nblk(1), nbn*nblk(2) )
    act  = pattern ( kind, nbn, nblk )
    ref  = masked ( nbn, nblk, act, phi0 )
    eo   = oracle_measure ( nbn, nblk, act, ref )
    call oracle_reinit ( nbn, nblk, act, ref, iters, beta )
    out = phi0
    out = SENTINEL
    call band_run ( ex, nbn, nblk, [ (0, b = 1, size(act)) ], act, phi0, iters, beta, out, e, code )
    call check ( code == MORPH_OK .and. same( out, ref ) .and. same1( e, eo ) .and. &
                 .not. same( out, masked(nbn, nblk, act, phi0) ), name )
  end subroutine one_rank_band


  subroutine reinit_band ()
    call one_rank_band ( 5, [3,3], 3, 4, BETA, 'band, one block of 5 alone: the oracle''s bits, E too' )
    call one_rank_band ( 3, [9,8], 4, 4, BETA, 'band, a checkerboard of blocks of 3: the oracle''s bits, E too' )
    call one_rank_band ( 5, [7,6], 5, 4, BETA, 'band, the frame of blocks of 5: the oracle''s bits, E too' )
    call one_rank_band ( 5, [8,7], 6, 4, BETA, 'band, scattered blocks of 5: the oracle''s bits, E too' )
    call one_rank_band ( 5, [8,7], 6, 3, 0.0_R8, 'band, scattered blocks of 5, no correction: the oracle''s bits' )
    call one_rank_band ( 8, [8,8], 7, 4, BETA, 'band, a ring of blocks of 8: the oracle''s bits, E too' )
  end subroutine reinit_band


  !> On a ring of blocks of 8, from the same phi: after 2 and 4 steps, the
  !> nodes farther than 4 iters (L1, in nodes) from every inactive node have
  !> the bits of the whole grid, and there are such nodes.
  subroutine reach ()
    type(morph_exchange_serial_t) :: ex
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:), full(:,:), out(:,:)
    integer(I4), allocatable :: d(:,:)
    logical,  allocatable :: act(:)
    real(R8)    :: e
    integer(I4) :: nbn, nblk(2), nx, ny, i, j, b, i0, j0, code, it, nfar, ndiff, m
    logical     :: changed, ok

    nbn  = 8
    nblk = [10,10]
    nx   = nbn*nblk(1)
    ny   = nbn*nblk(2)
    phi0 = field ( nx, ny )
    act  = pattern ( 7, nbn, nblk )
    allocate( d(nx,ny) )
    d = nx + ny                                   ! farther than any node
    do b = 1, size(act)
      if ( act(b) ) cycle
      i0 = mod(b-1, nblk(1))*nbn
      j0 = ((b-1)/nblk(1))*nbn
      d(i0+1:i0+nbn, j0+1:j0+nbn) = 0
    enddo
    changed = .true.
    do while ( changed )                          ! the L1 distance, by relaxation
      changed = .false.
      do j = 1, ny
        do i = 1, nx
          m = d(i,j)
          if ( i > 1 )  m = min( m, d(i-1,j) + 1 )
          if ( i < nx ) m = min( m, d(i+1,j) + 1 )
          if ( j > 1 )  m = min( m, d(i,j-1) + 1 )
          if ( j < ny ) m = min( m, d(i,j+1) + 1 )
          if ( m < d(i,j) ) then
            d(i,j) = m
            changed = .true.
          endif
        enddo
      enddo
    enddo
    ok = any( .not. act ) .and. any( act )
    do it = 2, 4, 2
      full = phi0
      call Morph_LS_Reinit_Full ( full, H, it, BETA, st )
      out = phi0
      out = SENTINEL
      call band_run ( ex, nbn, nblk, [ (0, b = 1, size(act)) ], act, phi0, it, BETA, out, e, code )
      nfar  = 0
      ndiff = 0
      do j = 1, ny
        do i = 1, nx
          if ( d(i,j) > 4*it ) then
            nfar = nfar + 1
            if ( .not. same1( out(i,j), full(i,j) ) ) ndiff = ndiff + 1
          endif
        enddo
      enddo
      if ( .not. ( st%ok() .and. code == MORPH_OK .and. nfar > 0 .and. ndiff == 0 ) ) ok = .false.
    enddo
    call check ( ok, 'a ring of blocks of 8, 2 and 4 steps: farther than 4 iters nodes from the inactive blocks, &
                 &the bits of the whole grid' )
  end subroutine reach


  subroutine if_needed_one_rank ()
    type(morph_exchange_serial_t) :: ex
    type(morph_ls_grid_t) :: g
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:), phi(:,:,:), ref(:,:,:), orig(:,:,:)
    real(R8)    :: e, e0
    logical     :: done, above, at, below
    integer(I4) :: b

    phi0 = field ( 40, 35 )
    call Morph_LS_Grid_Init ( g, 5, [8,7], H, [ (0, b = 1, 56) ], ex, st, pattern(6, 5, [8,7]) )
    ref = slots_of ( g, phi0 )
    call Morph_LS_Eikonal_Error_Grid ( g, ex, ref, BAND, e0, st )
    call Morph_LS_Reinit_Grid ( g, ex, ref, 4, BETA, st )

    orig = slots_of ( g, phi0 )
    phi  = orig
    call Morph_LS_Reinit_If_Needed_Grid ( g, ex, phi, BAND, 2.0_R8*e0, 4, BETA, done, e, st )
    below = st%ok() .and. .not. done .and. same1( e, e0 ) .and. same3( phi, orig )
    call Morph_LS_Reinit_If_Needed_Grid ( g, ex, phi, BAND, e0, 4, BETA, done, e, st )
    at = st%ok() .and. .not. done .and. same1( e, e0 ) .and. same3( phi, orig )
    call Morph_LS_Reinit_If_Needed_Grid ( g, ex, phi, BAND, 0.5_R8*e0, 4, BETA, done, e, st )
    above = st%ok() .and. done .and. same1( e, e0 ) .and. same3( phi, ref ) .and. .not. same3( phi, orig )
    call check ( below, 'if needed, E below the threshold: phi untouched, done false' )
    call check ( at, 'if needed, E equal to the threshold: phi untouched, done false' )
    call check ( above, 'if needed, E above the threshold: the bits of Morph_LS_Reinit_Grid, done true' )
  end subroutine if_needed_one_rank


  subroutine refusals_one_rank ()
    type(morph_exchange_serial_t) :: ex
    type(morph_ls_grid_t) :: g, g0
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:), phi(:,:,:), keep(:,:,:), bad(:,:,:)
    real(R8)    :: e
    integer(I4) :: b, k, shp(3,3)
    logical     :: done, okr, oke, okn

    phi0 = field ( 25, 25 )
    call Morph_LS_Grid_Init ( g, 5, [5,5], H, [ (0, b = 1, 25) ], ex, st )
    phi  = slots_of ( g, phi0 )
    keep = phi

    shp = reshape( [4,5,25, 5,4,25, 5,5,24], [3,3] )
    okr = .true.
    oke = .true.
    do k = 1, 3
      if ( allocated(bad) ) deallocate( bad )
      allocate( bad(shp(1,k), shp(2,k), shp(3,k)) )
      bad = 1.0_R8
      call Morph_LS_Reinit_Grid ( g, ex, bad, 4, BETA, st )
      if ( .not. ( st%code == MORPH_ERR_SHAPE .and. all( bad == 1.0_R8 ) ) ) okr = .false.
      if ( k == 3 .and. index(st%message, 'phi is 5 x 5 x 24') == 0 ) okr = .false.
      e = 1.0_R8
      call Morph_LS_Eikonal_Error_Grid ( g, ex, bad, BAND, e, st )
      if ( .not. ( st%code == MORPH_ERR_SHAPE .and. same1( e, 0.0_R8 ) ) ) oke = .false.
    enddo
    call check ( okr, 'reinitialisation of a phi of 4 x 5 x 25, 5 x 4 x 25, 5 x 5 x 24 for 5 x 5 x 25: refused, &
                 &phi untouched' )
    call check ( oke, 'E of those shapes: refused, E = 0' )
    call Morph_LS_Reinit_Grid ( g, ex, bad, -1, BETA, st )
    call check ( st%code == MORPH_ERR_SHAPE .and. all( bad == 1.0_R8 ), &
                 'reinitialisation of a phi of 5 x 5 x 24 with iters = -1: the shape reported' )

    call Morph_LS_Reinit_Grid ( g, ex, phi, -1, BETA, st )
    call check ( st%code == MORPH_ERR_INPUT .and. index(st%message, 'iters = -1') > 0 .and. same3( phi, keep ), &
                 'reinitialisation with iters = -1: refused, phi untouched' )
    call Morph_LS_Reinit_Grid ( g, ex, phi, 4, -0.5_R8, st )
    call check ( st%code == MORPH_ERR_INPUT .and. same3( phi, keep ), &
                 'reinitialisation with beta = -0.5: refused, phi untouched' )
    call Morph_LS_Reinit_Grid ( g, ex, phi, 4, ieee_value(0.0_R8, ieee_quiet_nan), st )
    call check ( st%code == MORPH_ERR_INPUT .and. same3( phi, keep ), &
                 'reinitialisation with beta = NaN: refused, phi untouched' )

    g%nranks = 2
    call Morph_LS_Reinit_Grid ( g, ex, phi, 4, BETA, st )
    call check ( st%code == MORPH_ERR_EXCHANGE .and. index(st%message, 'the exchange is rank 0 of 1') > 0 .and. &
                 same3( phi, keep ), 'reinitialisation with an exchange of another number of ranks: refused, &
                 &phi untouched' )
    e = 1.0_R8
    call Morph_LS_Eikonal_Error_Grid ( g, ex, phi, BAND, e, st )
    call check ( st%code == MORPH_ERR_EXCHANGE .and. same1( e, 0.0_R8 ), &
                 'E with an exchange of another number of ranks: refused, E = 0' )
    call Morph_LS_Reinit_If_Needed_Grid ( g, ex, phi, BAND, 0.0_R8, 4, BETA, done, e, st )
    call check ( st%code == MORPH_ERR_EXCHANGE .and. .not. done .and. same3( phi, keep ), &
                 'if needed with that exchange: refused, phi untouched, done false' )
    g%nranks = 1
    g%me = 1
    call Morph_LS_Reinit_Grid ( g, ex, phi, 4, BETA, st )
    call check ( st%code == MORPH_ERR_EXCHANGE .and. same3( phi, keep ), &
                 'reinitialisation with an exchange of another rank: refused, phi untouched' )

    ! a grid whose Init was refused (2 nodes per block side): refused before
    ! any size of its arrays is taken
    call Morph_LS_Grid_Init ( g0, 2, [3,3], H, [ (0, b = 1, 9) ], ex, st )
    okn = .not. st%ok()
    e = 1.0_R8
    call Morph_LS_Eikonal_Error_Grid ( g0, ex, phi, BAND, e, st )
    okn = okn .and. st%code == MORPH_ERR_SHAPE .and. index(st%message, 'not initialised') > 0 .and. same1( e, 0.0_R8 )
    call Morph_LS_Reinit_Grid ( g0, ex, phi, 4, BETA, st )
    okn = okn .and. st%code == MORPH_ERR_SHAPE .and. index(st%message, 'not initialised') > 0 .and. same3( phi, keep )
    call Morph_LS_Reinit_If_Needed_Grid ( g0, ex, phi, BAND, 0.0_R8, 4, BETA, done, e, st )
    call check ( okn .and. st%code == MORPH_ERR_SHAPE .and. index(st%message, 'not initialised') > 0 .and. &
                 .not. done .and. same3( phi, keep ), &
                 'a grid whose Init was refused: E, the reinitialisation and if needed refuse it, phi untouched' )
  end subroutine refusals_one_rank


  !> One rank: each collective of each kind that the calls make made to fail
  !> in turn (phi untouched, the failure returned), then one past the last
  !> of each kind (the call done).
  subroutine failing_exchange ()
    type(failing_exchange_t) :: ex
    type(morph_ls_grid_t) :: g
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:), phi(:,:,:), keep(:,:,:), ref(:,:,:)
    real(R8)    :: e, e0
    integer(I4) :: b, k, kind
    integer(I4), parameter :: NR(3) = [1, 5, 0]    ! the reinitialisation, 2 steps: 1 allgather, 1 + 2 x 2 swaps
    integer(I4), parameter :: NE(3) = [1, 1, 2]    ! E: 1 allgather, 1 swap, 2 ordered sums
    integer(I4), parameter :: NI(3) = [3, 6, 2]    ! if needed, E above: E, its allgather, the reinitialisation
    logical     :: okf, okp, done, ok0

    phi0 = field ( 30, 25 )
    call Morph_LS_Grid_Init ( g, 5, [6,5], H, [ (0, b = 1, 30) ], ex, st, pattern(6, 5, [6,5]) )
    keep = slots_of ( g, phi0 )
    ref  = keep
    call Morph_LS_Eikonal_Error_Grid ( g, ex, keep, BAND, e0, st )
    call Morph_LS_Reinit_Grid ( g, ex, ref, 2, BETA, st )

    okf = .true.
    okp = .true.
    do kind = 1, 3
      do k = 1, NR(kind) + 1
        ex%fail = 0
        ex%calls = 0
        ex%fail(kind) = k
        phi = keep
        call Morph_LS_Reinit_Grid ( g, ex, phi, 2, BETA, st )
        if ( k <= NR(kind) ) then
          if ( .not. ( st%code == MORPH_ERR_EXCHANGE .and. same3( phi, keep ) ) ) okf = .false.
        else
          if ( .not. ( st%ok() .and. same3( phi, ref ) .and. all( ex%calls == NR ) ) ) okp = .false.
        endif
      enddo
    enddo
    call check ( okf, 'reinitialisation, 2 steps, its allgather or any of its 5 swaps failing: returned, &
                 &phi untouched' )
    call check ( okp, 'reinitialisation, 2 steps: 1 allgather, 5 swaps, no ordered sum (one more fails nothing)' )
    ex%fail = 0
    ex%calls = 0
    phi = keep
    call Morph_LS_Reinit_Grid ( g, ex, phi, 0, BETA, st )
    ok0 = st%ok() .and. same3( phi, keep ) .and. all( ex%calls == [1, 0, 0] )
    ex%calls = 0
    call Morph_LS_Reinit_Grid ( g, ex, phi, 2, 0.0_R8, st )
    call check ( ok0 .and. st%ok() .and. all( ex%calls == [1, 3, 0] ), 'reinitialisation, 0 steps: 1 allgather &
                 &and nothing else; 2 steps without the correction: 1 allgather, 3 swaps' )

    okf = .true.
    okp = .true.
    do kind = 1, 3
      do k = 1, NE(kind) + 1
        ex%fail = 0
        ex%calls = 0
        ex%fail(kind) = k
        e = 1.0_R8
        call Morph_LS_Eikonal_Error_Grid ( g, ex, keep, BAND, e, st )
        if ( k <= NE(kind) ) then
          if ( .not. ( st%code == MORPH_ERR_EXCHANGE .and. same1( e, 0.0_R8 ) ) ) okf = .false.
        else
          if ( .not. ( st%ok() .and. same1( e, e0 ) .and. all( ex%calls == NE ) ) ) okp = .false.
        endif
      enddo
    enddo
    call check ( okf, 'E, its allgather, its swap or either ordered sum failing: returned, E = 0' )
    call check ( okp, 'E: 1 allgather, 1 swap, 2 ordered sums (one more fails nothing)' )

    okf = .true.
    okp = .true.
    do kind = 1, 3
      do k = 1, NI(kind) + 1
        ex%fail = 0
        ex%calls = 0
        ex%fail(kind) = k
        phi = keep
        call Morph_LS_Reinit_If_Needed_Grid ( g, ex, phi, BAND, 0.5_R8*e0, 2, BETA, done, e, st )
        if ( k <= NI(kind) ) then
          if ( .not. ( st%code == MORPH_ERR_EXCHANGE .and. .not. done .and. same3( phi, keep ) ) ) okf = .false.
        else
          if ( .not. ( st%ok() .and. done .and. same3( phi, ref ) .and. all( ex%calls == NI ) ) ) okp = .false.
        endif
      enddo
    enddo
    call check ( okf, 'if needed, E above, any of its collectives failing: returned, phi untouched, done false' )
    call check ( okp, 'if needed, E above: 3 allgathers, 6 swaps, 2 ordered sums (one more fails nothing)' )
  end subroutine failing_exchange


  !> n simulated ranks, owner(b), a band of kind, 4 steps: out and E of
  !> every rank against one rank's.
  subroutine ranks_case ( n, nbn, nblk, owner, kind, name )
    integer(I4),      intent(in) :: n, nbn, nblk(2), owner(:), kind
    character(len=*), intent(in) :: name
    type(morph_sim_world_t), target :: world
    type(morph_exchange_serial_t)   :: ex1
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:), one(:,:), out(:,:)
    logical,  allocatable :: act(:)
    real(R8)    :: e1, es(0:n-1)
    integer(I4) :: codes(0:n-1), code1, b

    phi0 = field ( nbn*nblk(1), nbn*nblk(2) )
    act  = pattern ( kind, nbn, nblk )
    one  = phi0
    one  = SENTINEL
    call band_run ( ex1, nbn, nblk, [ (0, b = 1, size(act)) ], act, phi0, 4, BETA, one, e1, code1 )
    out = one
    out = SENTINEL
    codes = -1
    es = -1.0_R8
    call world%init ( n, st )
    !$omp parallel num_threads(n)
    block
      type(morph_exchange_sim_t) :: ex
      integer(I4) :: t, c
      real(R8) :: e
      t = 0
      !$ t = omp_get_thread_num()
      ex%me = t
      ex%world => world
      call band_run ( ex, nbn, nblk, owner, act, phi0, 4, BETA, out, e, c )
      codes(t) = c
      es(t) = e
    end block
    !$omp end parallel
    call check ( st%ok() .and. code1 == MORPH_OK .and. all( codes == MORPH_OK ) .and. same( out, one ) .and. &
                 all( transfer(es, 0_I8, n) == transfer(e1, 0_I8) ), name )
  end subroutine ranks_case


  subroutine ranks ()
    call ranks_case ( 2, 5, [8,7], owners(2, 56, 1), 1, &
                      '2 ranks, every block of 5: the bits of one rank, E the same on both' )
    call ranks_case ( 2, 5, [8,7], owners(2, 56, 1), 6, &
                      '2 ranks, scattered blocks of 5: the bits of one rank, E the same on both' )
    call ranks_case ( 4, 5, [8,7], owners(4, 56, 3), 6, &
                      '4 ranks, scattered blocks of 5: the bits of one rank, E the same on all' )
    call ranks_case ( 3, 3, [13,12], owners(3, 156, 1), 4, &
                      '3 ranks, a checkerboard of blocks of 3: the bits of one rank, E the same on all' )
    call ranks_case ( 4, 8, [8,8], owners(4, 64, 5), 7, &
                      '4 ranks, a ring of blocks of 8: the bits of one rank, E the same on all' )
    call ranks_case ( 3, 5, [8,7], owners_no_active(3, pattern(6, 5, [8,7])), 6, &
                      '3 ranks, rank 2 with no active block: the bits of one rank, E the same on all' )
  end subroutine ranks


  !> 2 simulated ranks, E above the threshold: both reinitialise, with the
  !> bits and E of one rank.
  subroutine if_needed_ranks ()
    type(morph_sim_world_t), target :: world
    type(morph_exchange_serial_t)   :: ex1
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:), one(:,:), out(:,:)
    logical,  allocatable :: act(:)
    real(R8)    :: e1, es(0:1)
    integer(I4) :: code1, b
    logical     :: dn(0:1)

    phi0 = field ( 40, 35 )
    act  = pattern ( 6, 5, [8,7] )
    one  = phi0
    one  = SENTINEL
    call band_run ( ex1, 5, [8,7], [ (0, b = 1, 56) ], act, phi0, 4, BETA, one, e1, code1 )
    out = one
    out = SENTINEL
    dn = .false.
    es = -1.0_R8
    call world%init ( 2, st )
    !$omp parallel num_threads(2)
    block
      type(morph_exchange_sim_t) :: ex
      type(morph_ls_grid_t) :: g
      type(morph_status_t)  :: s
      real(R8), allocatable :: phi(:,:,:)
      integer(I4) :: t
      real(R8) :: e
      logical :: done
      t = 0
      !$ t = omp_get_thread_num()
      ex%me = t
      ex%world => world
      call Morph_LS_Grid_Init ( g, 5, [8,7], H, owners(2, 56, 1), ex, s, act )
      phi = slots_of ( g, phi0 )
      call Morph_LS_Reinit_If_Needed_Grid ( g, ex, phi, BAND, 0.5_R8*e1, 4, BETA, done, e, s )
      dn(t) = done .and. s%ok()
      es(t) = e
      call Morph_LS_Grid_Gather ( g, phi, out, s )
    end block
    !$omp end parallel
    call check ( st%ok() .and. code1 == MORPH_OK .and. all( dn ) .and. &
                 all( transfer(es, 0_I8, 2) == transfer(e1, 0_I8) ) .and. same( out, one ), &
                 '2 ranks, if needed with E above the threshold: both reinitialise, the bits and E of one rank' )
  end subroutine if_needed_ranks


  !> 2 simulated ranks: rank 1 gives iters = -1; then the ranks give other
  !> iters, other beta, other thresholds.
  subroutine refusals_ranks ()
    type(morph_sim_world_t), target :: world
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:)
    logical,  allocatable :: act(:)
    integer(I4) :: codes(0:1,4)
    logical     :: kept(0:1,4), msgs(0:1,4)

    phi0 = field ( 40, 35 )
    act  = pattern ( 6, 5, [8,7] )
    codes = -1
    kept  = .false.
    msgs  = .false.
    call world%init ( 2, st )
    !$omp parallel num_threads(2)
    block
      type(morph_exchange_sim_t) :: ex
      type(morph_ls_grid_t) :: g
      type(morph_status_t)  :: s
      real(R8), allocatable :: phi(:,:,:), keep(:,:,:)
      integer(I4) :: t
      real(R8) :: e
      logical :: done
      t = 0
      !$ t = omp_get_thread_num()
      ex%me = t
      ex%world => world
      call Morph_LS_Grid_Init ( g, 5, [8,7], H, owners(2, 56, 1), ex, s, act )
      phi  = slots_of ( g, phi0 )
      keep = phi
      call Morph_LS_Reinit_Grid ( g, ex, phi, merge(-1, 4, t == 1), BETA, s )
      codes(t,1) = s%code
      kept(t,1)  = same3( phi, keep )
      if ( t == 1 ) then
        msgs(t,1) = index(s%message, 'iters = -1') > 0
      else
        msgs(t,1) = index(s%message, 'refused on rank 1') > 0
      endif
      call Morph_LS_Reinit_Grid ( g, ex, phi, 3 + t, BETA, s )
      codes(t,2) = s%code
      kept(t,2)  = same3( phi, keep )
      msgs(t,2)  = index(s%message, 'rank 1 was given other arguments') > 0
      call Morph_LS_Reinit_Grid ( g, ex, phi, 4, BETA + 0.25_R8*real(t, R8), s )
      codes(t,3) = s%code
      kept(t,3)  = same3( phi, keep )
      msgs(t,3)  = index(s%message, 'rank 1 was given other arguments') > 0
      call Morph_LS_Reinit_If_Needed_Grid ( g, ex, phi, BAND, 0.01_R8*real(t+1, R8), 4, BETA, done, e, s )
      codes(t,4) = s%code
      kept(t,4)  = same3( phi, keep ) .and. .not. done
      msgs(t,4)  = index(s%message, 'rank 1 was given other arguments') > 0
    end block
    !$omp end parallel
    call check ( st%ok() .and. all( codes(:,1) == MORPH_ERR_INPUT ) .and. all( kept(:,1) ) .and. all( msgs(:,1) ), &
                 '2 ranks, rank 1 gives iters = -1: both refuse (its message there, its code here), phi untouched' )
    call check ( all( codes(:,2) == MORPH_ERR_INPUT ) .and. all( kept(:,2) ) .and. all( msgs(:,2) ), &
                 '2 ranks given 3 and 4 steps: both refuse, phi untouched' )
    call check ( all( codes(:,3) == MORPH_ERR_INPUT ) .and. all( kept(:,3) ) .and. all( msgs(:,3) ), &
                 '2 ranks given other beta: both refuse, phi untouched' )
    call check ( all( codes(:,4) == MORPH_ERR_INPUT ) .and. all( kept(:,4) ) .and. all( msgs(:,4) ), &
                 '2 ranks given other thresholds: both refuse, phi untouched, done false' )
  end subroutine refusals_ranks



  !> 3 simulated ranks: rank 1 gives a phi of one slot too few, rank 2
  !> iters = -1. Each refuses with its own code and message; rank 0 with the
  !> first refusing rank's code (rank 1's, MORPH_ERR_SHAPE), phi untouched.
  subroutine two_refusals_ranks ()
    type(morph_sim_world_t), target :: world
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:)
    logical,  allocatable :: act(:)
    integer(I4) :: codes(0:2)
    logical     :: msgs(0:2), kept(0:2)

    phi0 = field ( 40, 35 )
    act  = pattern ( 6, 5, [8,7] )
    codes = -1
    msgs  = .false.
    kept  = .false.
    call world%init ( 3, st )
    !$omp parallel num_threads(3)
    block
      type(morph_exchange_sim_t) :: ex
      type(morph_ls_grid_t) :: g
      type(morph_status_t)  :: s
      real(R8), allocatable :: phi(:,:,:), keep(:,:,:)
      integer(I4) :: t
      t = 0
      !$ t = omp_get_thread_num()
      ex%me = t
      ex%world => world
      call Morph_LS_Grid_Init ( g, 5, [8,7], H, owners(3, 56, 1), ex, s, act )
      if ( t == 1 ) then
        allocate( phi(5,5,size(g%blk)-1) )
        phi = 1.0_R8
      else
        phi = slots_of ( g, phi0 )
      endif
      keep = phi
      call Morph_LS_Reinit_Grid ( g, ex, phi, merge(-1, 4, t == 2), BETA, s )
      codes(t) = s%code
      kept(t)  = same3( phi, keep )
      select case ( t )
      case ( 0 )
        msgs(t) = index(s%message, 'refused on rank 1') > 0
      case ( 1 )
        msgs(t) = index(s%message, 'phi is 5 x 5 x') > 0
      case default
        msgs(t) = index(s%message, 'iters = -1') > 0
      end select
    end block
    !$omp end parallel
    call check ( st%ok() .and. codes(0) == MORPH_ERR_SHAPE .and. codes(1) == MORPH_ERR_SHAPE .and. &
                 codes(2) == MORPH_ERR_INPUT .and. all( msgs ) .and. all( kept ), &
                 '3 ranks, rank 1 gives a phi too small, rank 2 iters = -1: each its own refusal, rank 0 rank 1''s' )
  end subroutine two_refusals_ranks


  !> 2 simulated ranks, the Init of rank 1 refused (2 nodes per block side):
  !> rank 1 refuses as not initialised, rank 0 with its code; phi untouched.
  subroutine uninitialised_rank ()
    type(morph_sim_world_t), target :: world
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:)
    logical,  allocatable :: act(:)
    integer(I4) :: codes(0:1)
    logical     :: msgs(0:1), kept(0:1), inits(0:1)

    phi0 = field ( 40, 35 )
    act  = pattern ( 6, 5, [8,7] )
    codes = -1
    msgs  = .false.
    kept  = .false.
    inits = .false.
    call world%init ( 2, st )
    !$omp parallel num_threads(2)
    block
      type(morph_exchange_sim_t) :: ex
      type(morph_ls_grid_t) :: g
      type(morph_status_t)  :: s
      real(R8), allocatable :: phi(:,:,:), keep(:,:,:)
      integer(I4) :: t
      t = 0
      !$ t = omp_get_thread_num()
      ex%me = t
      ex%world => world
      call Morph_LS_Grid_Init ( g, merge(2, 5, t == 1), [8,7], H, owners(2, 56, 1), ex, s, act )
      inits(t) = s%ok() .neqv. ( t == 1 )
      if ( t == 1 ) then
        allocate( phi(5,5,3) )
        phi = 1.0_R8
      else
        phi = slots_of ( g, phi0 )
      endif
      keep = phi
      call Morph_LS_Reinit_Grid ( g, ex, phi, 4, BETA, s )
      codes(t) = s%code
      kept(t)  = same3( phi, keep )
      if ( t == 1 ) then
        msgs(t) = index(s%message, 'not initialised') > 0
      else
        msgs(t) = index(s%message, 'refused on rank 1') > 0
      endif
    end block
    !$omp end parallel
    call check ( st%ok() .and. all( inits ) .and. all( codes == MORPH_ERR_SHAPE ) .and. all( msgs ) .and. &
                 all( kept ), '2 ranks, the Init of rank 1 refused: both refuse (not initialised there), phi untouched' )
  end subroutine uninitialised_rank

end program test_morph_ls_reinit_blocks
