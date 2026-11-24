!>@brief Unit test of the level-set field on a grid of blocks split among
!> ranks (Morph_LS_Grid) and of the exchanges it uses (Morph_Exchange,
!> Morph_Exchange_Sim). The reference is the step on the whole grid,
!> Morph_LS_HJ_Step_Full, which test_morph_ls checks against the theory; here
!> every comparison is bit for bit (the 64 bits of each value):
!>   * the serial exchange: copies, the sum in block order, its refusals;
!>   * the simulated ranks, one thread of a parallel region per rank: the
!>     swap, the allgather and the ordered sum give what the serial exchange
!>     gives, over thousands of rounds whose counts change; a team without one
!>     thread per rank, a rank that is not its thread, and counts or owners
!>     that do not match are refused on every rank, none left waiting;
!>   * the grid: its refusals; scatter then gather gives the whole grid back,
!>     a gather writes the nodes of its own blocks only, and the slot of each
!>     block is its place on its rank and 0 on the others;
!>   * the ghosts after one fill: the values of the whole grid (the
!>     neighbours' nodes inside, the extrapolation at the edges), the corners
!>     untouched;
!>   * the step: the bits of the step on the whole grid, for blocks of 3, 5
!>     and 8 nodes, 1 to 20 blocks, 1 to 4 ranks split in six ways, one of
!>     them a rank with no block; V = 0 on one rank's blocks only; V = 0
!>     everywhere leaves phi untouched, a -0 under a V of -0 too (the stages
!>     would make it +0); refused arguments, on one rank or on all, leave phi
!>     untouched on every rank, each refusing rank reporting its own code and
!>     the others the first refusing rank's; another rank's grid is refused.
!> Built without OpenMP, the simulated ranks are skipped, and the output says
!> so. The values printed must be the same at 1 and 8 threads.
program test_morph_ls_blocks
  use iso_fortran_env, only: I4 => int32, I8 => int64, R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  !$ use omp_lib, only: omp_get_thread_num
  use Morph_Types_m
  use Morph_Exchange
  use Morph_Exchange_Sim
  use Morph_LS_HJ,   only: MORPH_LS_NG, Morph_LS_Ghosts_Extrap_All, Morph_LS_HJ_Step_Full
  use Morph_LS_Grid

  implicit none

  integer(I4), parameter :: NG = MORPH_LS_NG
  real(R8),    parameter :: H  = 1.0_R8 / 32.0_R8
  real(R8),    parameter :: DT = 0.4_R8 * H
  real(R8),    parameter :: SENTINEL = -7.0_R8      ! no field value of the test
  integer(I4), parameter :: NSTEP = 4

  integer(I4) :: nfail
  logical     :: omp

  nfail = 0
  omp = .false.
  !$ omp = .true.
  write(*,'(A)') '== Morph level-set grid of blocks and its exchanges =='

  call serial_exchange ()
  call grid_refusals ()
  call scatter_gather ()
  call steps_one_rank ()
  call step_refusals ()
  if ( omp ) then
    call sim_collectives ()
    call sim_stress ()
    call sim_refusals ()
    call steps_ranks ()
    call refused_on_one_rank ()
    call refused_on_two_ranks ()
  else
    write(*,'(A)') '   SKIP the simulated ranks: built without OpenMP'
  endif

  if ( nfail > 0 ) then
    write(*,'(I0,A)') nfail, ' BLOCK-GRID CHECK(S) FAILED'
    error stop 1
  endif
  write(*,'(A)') 'RESULT: all block-grid checks PASS'

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


  !> NSTEP steps on the whole grid.
  function reference ( phi0, vls ) result ( phi )
    real(R8), intent(in)  :: phi0(:,:), vls(:,:)
    real(R8), allocatable :: phi(:,:)
    type(morph_status_t)  :: st
    integer(I4) :: k
    phi = phi0
    do k = 1, NSTEP
      call Morph_LS_HJ_Step_Full ( phi, vls, H, DT, st )
    enddo
  end function reference


  !> The whole grid with its ghosts, extrapolated.
  function with_ghosts ( phi ) result ( wf )
    real(R8), intent(in)  :: phi(:,:)
    real(R8), allocatable :: wf(:,:)
    allocate( wf(1-NG:size(phi,1)+NG, 1-NG:size(phi,2)+NG) )
    wf = 0.0_R8
    wf(1:size(phi,1), 1:size(phi,2)) = phi
    call Morph_LS_Ghosts_Extrap_All ( wf )
  end function with_ghosts


  !> The ghosts of side d of a block, (i1, i2, j1, j2), written here again
  !> from the layout of the kernels.
  function side_box ( nbn, d ) result ( g )
    integer(I4), intent(in) :: nbn, d
    integer(I4) :: g(4)
    g = [ 1, nbn, 1, nbn ]
    select case ( d )
    case ( 1 ); g(1:2) = [ 1-NG, 0 ]
    case ( 2 ); g(1:2) = [ nbn+1, nbn+NG ]
    case ( 3 ); g(3:4) = [ 1-NG, 0 ]
    case ( 4 ); g(3:4) = [ nbn+1, nbn+NG ]
    end select
  end function side_box


  !> After one fill: the four sides of every block of this rank hold the
  !> values of the whole grid wf at the same nodes, the four corners still
  !> hold SENTINEL.
  logical function ghosts_ok ( grid, w, wf ) result ( ok )
    type(morph_ls_grid_t), intent(in) :: grid
    real(R8), intent(in) :: w(1-NG:,1-NG:,:), wf(1-NG:,1-NG:)
    integer(I4) :: s, d, b, i0, j0, g(4), n
    real(R8), allocatable :: corner(:,:)
    n = grid%nbn
    allocate( corner(NG,NG) )
    corner = SENTINEL
    ok = .true.
    do s = 1, size(grid%blk)
      b  = grid%blk(s)
      i0 = mod(b-1, grid%nblk(1)) * n
      j0 = ((b-1)/grid%nblk(1)) * n
      do d = 1, 4
        g = side_box ( n, d )
        if ( .not. same( w(g(1):g(2), g(3):g(4), s), wf(i0+g(1):i0+g(2), j0+g(3):j0+g(4)) ) ) ok = .false.
      enddo
      if ( .not. ( same( w(1-NG:0, 1-NG:0, s), corner )     .and. same( w(n+1:n+NG, 1-NG:0, s), corner ) .and. &
                   same( w(1-NG:0, n+1:n+NG, s), corner )   .and. same( w(n+1:n+NG, n+1:n+NG, s), corner ) ) ) &
        ok = .false.
    enddo
  end function ghosts_ok


  !> One rank's part: its grid, its blocks of phi0 and v0, one fill of the
  !> ghosts (halo: they are those of wf), NSTEP steps, its blocks gathered
  !> into out (shared by the ranks). Stops at the first refusal.
  subroutine advance ( exch, nbn, nblk, owner, phi0, v0, wf, out, halo, st )
    class(morph_exchange_t), intent(inout) :: exch
    integer(I4),             intent(in)    :: nbn, nblk(2), owner(:)
    real(R8),                intent(in)    :: phi0(:,:), v0(:,:), wf(1-NG:,1-NG:)
    real(R8),                intent(inout) :: out(:,:)
    logical,                 intent(out)   :: halo
    type(morph_status_t),    intent(out)   :: st
    type(morph_ls_grid_t) :: grid
    real(R8), allocatable :: phi(:,:,:), vb(:,:,:), w(:,:,:)
    integer(I4) :: k, ns

    halo = .false.
    call Morph_LS_Grid_Init ( grid, nbn, nblk, H, owner, exch, st )
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
    halo = ghosts_ok ( grid, w, wf )
    do k = 1, NSTEP
      call Morph_LS_HJ_Step ( grid, exch, phi, vb, DT, st )
      if ( .not. st%ok() ) return
    enddo
    call Morph_LS_Grid_Gather ( grid, phi, out, st )
  end subroutine advance


  ! ===========================================================================
  subroutine serial_exchange ()
    type(morph_exchange_serial_t) :: ex
    type(morph_status_t) :: st, s1, s2, s3, s4
    real(R8) :: r(4), q(4), rbig(8), qbig(8), tot
    type(morph_status_t) :: s5
    logical :: ok

    q    = 1.0_R8
    qbig = 1.0_R8
    call check ( ex%rank() == 0 .and. ex%nranks() == 1, 'serial exchange: rank 0 of 1' )

    r = 0.0_R8
    call ex%swap ( [3], [1.5_R8, -2.0_R8, 0.25_R8, 9.0_R8], [3], r, st )
    ok = st%ok() .and. all( r == [1.5_R8, -2.0_R8, 0.25_R8, 0.0_R8] )
    r = 0.0_R8
    call ex%allgather ( [4.0_R8, 5.0_R8], r, [2], s1 )
    call check ( ok .and. s1%ok() .and. all( r == [4.0_R8, 5.0_R8, 0.0_R8, 0.0_R8] ), &
                 'serial exchange: swap and allgather copy the counted values' )

    ! the order matters: (1e16 + 1) - 1e16 + 1 = 1 in block order, 2 or 0 in others
    call ex%ordered_sum ( 4, [0,0,0,0], [1.0e16_R8, 1.0_R8, -1.0e16_R8, 1.0_R8], tot, st )
    call check ( st%ok() .and. tot == 1.0_R8, 'serial exchange: the sum in block order (1e16, 1, -1e16, 1 -> 1)' )

    call ex%swap ( [1,1], q, [1,1], r, s1 )
    call ex%swap ( [2], q, [1], r, s2 )
    call ex%swap ( [5], q, [5], rbig, s3 )
    call ex%swap ( [-1], q, [-1], r, s4 )
    call ex%swap ( [5], qbig, [5], r, s5 )
    call check ( s1%code == MORPH_ERR_EXCHANGE .and. s2%code == MORPH_ERR_EXCHANGE .and. &
                 s3%code == MORPH_ERR_EXCHANGE .and. s4%code == MORPH_ERR_EXCHANGE .and. &
                 s5%code == MORPH_ERR_EXCHANGE, &
                 'serial exchange: swap refuses two ranks, counts that differ or are negative, more than '// &
                 'the send or the receive buffer' )
    call ex%allgather ( [1.0_R8], r, [1,1], s1 )
    call ex%allgather ( [1.0_R8], r, [2], s2 )
    call ex%allgather ( [1.0_R8, 2.0_R8, 3.0_R8, 4.0_R8, 5.0_R8], r, [5], s3 )
    call check ( s1%code == MORPH_ERR_EXCHANGE .and. s2%code == MORPH_ERR_EXCHANGE .and. &
                 s3%code == MORPH_ERR_EXCHANGE, &
                 'serial exchange: allgather refuses two ranks, a count not its values, more than the buffer' )
    call ex%ordered_sum ( 2, [0,0,0], [1.0_R8, 1.0_R8], tot, s1 )
    call ex%ordered_sum ( 2, [0,0], [1.0_R8], tot, s2 )
    call ex%ordered_sum ( 2, [0,1], [1.0_R8, 1.0_R8], tot, s3 )
    call ex%ordered_sum ( -1, [integer(I4) ::], [real(R8) ::], tot, s4 )
    call check ( s1%code == MORPH_ERR_EXCHANGE .and. s2%code == MORPH_ERR_EXCHANGE .and. &
                 s3%code == MORPH_ERR_EXCHANGE .and. s4%code == MORPH_ERR_EXCHANGE .and. tot == 0.0_R8, &
                 'serial exchange: ordered sum refuses owners or partials not one per block, a rank 1, -1 blocks' )
  end subroutine serial_exchange


  ! ===========================================================================
  subroutine grid_refusals ()
    type(morph_exchange_serial_t) :: ex
    type(morph_exchange_sim_t)    :: lone, two
    type(morph_sim_world_t), target :: world
    type(morph_ls_grid_t) :: grid
    type(morph_status_t)  :: st, s1, s2, s3, s4, s5
    real(R8) :: nan
    integer(I4) :: z6(6)

    nan = ieee_value(1.0_R8, ieee_quiet_nan)
    z6 = 0
    call Morph_LS_Grid_Init ( grid, 2, [2,3], H, z6, ex, s1 )
    call Morph_LS_Grid_Init ( grid, 8, [0,3], H, [integer(I4) ::], ex, s2 )
    call Morph_LS_Grid_Init ( grid, 8, [2,1], H, z6, ex, s3 )
    call Morph_LS_Grid_Init ( grid, 8, [2,0], H, [integer(I4) ::], ex, s4 )
    call check ( s1%code == MORPH_ERR_SHAPE .and. s2%code == MORPH_ERR_SHAPE .and. s3%code == MORPH_ERR_SHAPE &
                 .and. s4%code == MORPH_ERR_SHAPE, &
                 'grid: refuses blocks of 2 nodes, no block along x or y, owners not one per block' )
    call Morph_LS_Grid_Init ( grid, 8, [2,3], 0.0_R8, z6, ex, s1 )
    call Morph_LS_Grid_Init ( grid, 8, [2,3], nan, z6, ex, s2 )
    call Morph_LS_Grid_Init ( grid, 8, [2,3], H, [0,0,1,0,0,0], ex, s3 )
    call Morph_LS_Grid_Init ( grid, 8, [2,3], H, [0,0,0,0,-1,0], ex, s4 )
    call check ( s1%code == MORPH_ERR_INPUT .and. s2%code == MORPH_ERR_INPUT .and. &
                 s3%code == MORPH_ERR_INPUT .and. s4%code == MORPH_ERR_INPUT, &
                 'grid: refuses h = 0, h = NaN, a block on rank 1 of 1, a block on rank -1' )
    call world%init ( 2, st )
    two%world => world
    two%me = 2
    call Morph_LS_Grid_Init ( grid, 8, [2,3], H, z6, lone, s1 )
    call Morph_LS_Grid_Init ( grid, 8, [2,3], H, z6, two, s2 )
    two%me = -1
    call Morph_LS_Grid_Init ( grid, 8, [2,3], H, z6, two, s3 )
    call world%init ( 0, s5 )
    call check ( st%ok() .and. lone%nranks() == 0 .and. s1%code == MORPH_ERR_EXCHANGE .and. &
                 s2%code == MORPH_ERR_EXCHANGE .and. s3%code == MORPH_ERR_EXCHANGE .and. &
                 s5%code == MORPH_ERR_EXCHANGE .and. .not. allocated(grid%blk), &
                 'grid: refuses an exchange with no world, rank 2 or -1 of 2; a world of 0 ranks refused' )
  end subroutine grid_refusals


  ! ===========================================================================
  subroutine scatter_gather ()
    type(morph_exchange_serial_t) :: ex
    type(morph_exchange_sim_t)    :: r0, r1
    type(morph_sim_world_t), target :: world
    type(morph_ls_grid_t) :: grid, g0, g1
    type(morph_status_t)  :: st, s0, s1, s2, s3, s4, s5, s6, s7, s8
    real(R8), allocatable :: phi0(:,:), v0(:,:), full(:,:), b(:,:,:), b0(:,:,:), b1(:,:,:), wf(:,:), w(:,:,:)
    real(R8), allocatable :: bx(:,:,:), by(:,:,:)
    integer(I4) :: owner(12), i, j, k
    logical :: ok

    ! one rank, 4 x 3 blocks of 3 nodes
    call fields ( 12, 9, phi0, v0 )
    call Morph_LS_Grid_Init ( grid, 3, [4,3], H, [(0, k = 1, 12)], ex, st )
    allocate( b(3,3,12) )
    call Morph_LS_Grid_Scatter ( grid, phi0, b, s0 )
    full = 0.0_R8 * phi0
    call Morph_LS_Grid_Gather ( grid, b, full, s1 )
    call check ( st%ok() .and. s0%ok() .and. s1%ok() .and. same(full, phi0) .and. &
                 all( grid%blk == [(k, k = 1, 12)] ), &
                 'grid: one rank, 4 x 3 blocks of 3: scatter then gather gives the grid back' )

    ! the ghosts of one rank after one fill
    wf = with_ghosts ( phi0 )
    allocate( w(1-NG:3+NG, 1-NG:3+NG, 12) )
    w = SENTINEL
    w(1:3,1:3,:) = b
    call Morph_LS_Halo_Fill ( grid, ex, w, s2 )
    call check ( s2%ok() .and. ghosts_ok ( grid, w, wf ), &
                 'ghosts: one rank, 4 x 3 blocks of 3: the whole grid''s values, corners untouched' )

    ! two ranks built one after the other (the grid talks to no rank): a
    ! gather writes its own blocks only, and the two give the grid back
    owner = [0,1,1,0, 1,0,0,1, 0,0,1,1]
    call world%init ( 2, st )
    r0%world => world
    r1%world => world
    r0%me = 0
    r1%me = 1
    call Morph_LS_Grid_Init ( g0, 3, [4,3], H, [(0, k = 1, 12)], r0, s0 )
    call Morph_LS_Grid_Init ( g0, 3, [4,3], H, owner, r0, s0 )
    call Morph_LS_Grid_Init ( g1, 3, [4,3], H, owner, r1, s1 )
    allocate( b0(3,3,size(g0%blk)), b1(3,3,size(g1%blk)) )
    call Morph_LS_Grid_Scatter ( g0, phi0, b0, s2 )
    call Morph_LS_Grid_Scatter ( g1, phi0, b1, s3 )
    ok = st%ok() .and. s0%ok() .and. s1%ok() .and. s2%ok() .and. s3%ok() .and. &
         all( g0%blk == [1,4,6,7,9,10] ) .and. all( g1%blk == [2,3,5,8,11,12] )
    full = SENTINEL + 0.0_R8 * phi0
    call Morph_LS_Grid_Gather ( g0, b0, full, s0 )
    do j = 1, 9
      do i = 1, 12
        k = (i-1)/3 + 1 + 4*((j-1)/3)
        if ( owner(k) == 0 ) then
          if ( .not. same( full(i:i,j:j), phi0(i:i,j:j) ) ) ok = .false.
        else
          if ( full(i,j) /= SENTINEL ) ok = .false.
        endif
      enddo
    enddo
    call check ( ok .and. s0%ok(), 'grid: rank 0 of 2 gathers the nodes of its blocks only' )
    call check ( all( g0%slot([1,4,6,7,9,10]) == [1,2,3,4,5,6] ) .and. all( g0%slot([2,3,5,8,11,12]) == 0 ) &
                 .and. all( g1%slot([2,3,5,8,11,12]) == [1,2,3,4,5,6] ) .and. &
                 all( g1%slot([1,4,6,7,9,10]) == 0 ), &
                 'grid: the slot of each block on its rank, 0 on the other (a grid made again too)' )
    call Morph_LS_Grid_Gather ( g1, b1, full, s1 )
    call check ( s1%ok() .and. same(full, phi0), 'grid: the two ranks'' gathers give the grid back' )

    ! the strips one rank sends are those the other expects
    ok = all( g0%scount == [0, g1%rcount(0)] ) .and. all( g1%scount == [g0%rcount(1), 0] ) .and. &
         g0%rcount(1) > 0 .and. g1%rcount(0) > 0
    call check ( ok, 'grid: rank 0 sends rank 1 what rank 1 expects, and back' )

    ! refused shapes
    deallocate( b )
    allocate( b(3,3,11) )
    call Morph_LS_Grid_Scatter ( grid, phi0, b, s0 )
    call Morph_LS_Grid_Gather ( grid, b, full, s1 )
    deallocate( b )
    allocate( b(3,3,12) )
    call Morph_LS_Grid_Scatter ( grid, phi0(1:size(phi0,1)-1,:), b, s2 )
    call Morph_LS_Halo_Fill ( grid, ex, w(:,:,1:11), s3 )
    allocate( bx(2,3,12), by(3,2,12) )
    call Morph_LS_Grid_Scatter ( grid, phi0(:,1:size(phi0,2)-1), b, s4 )
    call Morph_LS_Grid_Scatter ( grid, phi0, bx, s5 )
    call Morph_LS_Grid_Scatter ( grid, phi0, by, s6 )
    call Morph_LS_Halo_Fill ( grid, ex, w(-2:5,:,:), s7 )
    call Morph_LS_Halo_Fill ( grid, ex, w(:,-2:5,:), s8 )
    call check ( s0%code == MORPH_ERR_SHAPE .and. s1%code == MORPH_ERR_SHAPE .and. &
                 s2%code == MORPH_ERR_SHAPE .and. s3%code == MORPH_ERR_SHAPE .and. &
                 s4%code == MORPH_ERR_SHAPE .and. s5%code == MORPH_ERR_SHAPE .and. &
                 s6%code == MORPH_ERR_SHAPE .and. s7%code == MORPH_ERR_SHAPE .and. &
                 s8%code == MORPH_ERR_SHAPE, &
                 'grid: scatter, gather and the ghosts refuse arrays of other shapes, axis by axis' )
  end subroutine scatter_gather


  ! ===========================================================================
  !> One rank, the serial exchange: the bits of the whole grid.
  subroutine steps_one_rank ()
    integer(I4), parameter :: NC = 6
    integer(I4) :: nbn(NC), nb1(NC), nb2(NC), c
    character(len=64) :: name

    nbn = [ 8, 8, 8, 3, 5, 3 ]
    nb1 = [ 1, 3, 5, 4, 2, 1 ]
    nb2 = [ 1, 2, 4, 3, 3, 1 ]
    do c = 1, NC
      write(name,'(A,I0,A,I0,A,I0)') 'step, one rank: ', nb1(c), ' x ', nb2(c), ' blocks of ', nbn(c)
      call one_rank_case ( nbn(c), [nb1(c), nb2(c)], trim(name) )
    enddo
    call zero_velocity_one_rank ()
  end subroutine steps_one_rank


  subroutine one_rank_case ( nbn, nblk, name )
    integer(I4),      intent(in) :: nbn, nblk(2)
    character(len=*), intent(in) :: name
    type(morph_exchange_serial_t) :: ex
    type(morph_status_t) :: st
    real(R8), allocatable :: phi0(:,:), v0(:,:), ref(:,:), out(:,:)
    integer(I4) :: k
    logical :: halo

    call fields ( nbn*nblk(1), nbn*nblk(2), phi0, v0 )
    ref = reference ( phi0, v0 )
    out = SENTINEL + 0.0_R8 * phi0
    call advance ( ex, nbn, nblk, [(0, k = 1, nblk(1)*nblk(2))], phi0, v0, with_ghosts(phi0), out, halo, st )
    call check ( st%ok() .and. halo .and. same(out, ref) .and. .not. same(ref, phi0), &
                 name//': the whole grid''s bits' )
  end subroutine one_rank_case


  !> V = 0 everywhere, -0 included, leaves phi untouched, a -0 too.
  subroutine zero_velocity_one_rank ()
    type(morph_exchange_serial_t) :: ex
    type(morph_ls_grid_t) :: grid
    type(morph_status_t)  :: st, s1
    real(R8), allocatable :: phi0(:,:), v0(:,:), b(:,:,:), vb(:,:,:), out(:,:)
    integer(I4) :: k

    call fields ( 16, 24, phi0, v0 )
    phi0(5,7) = -0.0_R8
    v0 = 0.0_R8
    v0(5,7) = -0.0_R8
    call Morph_LS_Grid_Init ( grid, 8, [2,3], H, [(0, k = 1, 6)], ex, st )
    allocate( b(8,8,6), vb(8,8,6) )
    call Morph_LS_Grid_Scatter ( grid, phi0, b, st )
    call Morph_LS_Grid_Scatter ( grid, v0, vb, st )
    call Morph_LS_HJ_Step ( grid, ex, b, vb, DT, s1 )
    out = 0.0_R8 * phi0
    call Morph_LS_Grid_Gather ( grid, b, out, st )
    call check ( st%ok() .and. s1%ok() .and. same(out, phi0), &
                 'step, one rank: V = 0 everywhere (a -0 among them) leaves phi untouched, its -0 too' )
  end subroutine zero_velocity_one_rank


  ! ===========================================================================
  !> Refused on one rank: the codes, phi untouched.
  subroutine step_refusals ()
    type(morph_exchange_serial_t) :: ex
    type(morph_exchange_sim_t)    :: r0
    type(morph_sim_world_t), target :: world
    type(morph_ls_grid_t) :: grid, g2, none
    type(morph_status_t)  :: st, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11
    real(R8), allocatable :: phi0(:,:), v0(:,:), b(:,:,:), vb(:,:,:), b0(:,:,:), bad(:,:,:)
    real(R8), allocatable :: e3(:,:,:), f3(:,:,:), e2(:,:), ew(:,:,:)
    type(morph_status_t)  :: z1, z2, z3, z4
    real(R8) :: nan
    integer(I4) :: k
    logical :: ok

    nan = ieee_value(1.0_R8, ieee_quiet_nan)
    call fields ( 16, 24, phi0, v0 )
    call Morph_LS_Grid_Init ( grid, 8, [2,3], H, [(0, k = 1, 6)], ex, st )
    allocate( b(8,8,6), vb(8,8,6), bad(8,8,5) )
    call Morph_LS_Grid_Scatter ( grid, phi0, b, st )
    call Morph_LS_Grid_Scatter ( grid, v0, vb, st )
    b0 = b
    bad = 1.0_R8
    call Morph_LS_HJ_Step ( grid, ex, bad, vb, DT, s1 )
    call Morph_LS_HJ_Step ( grid, ex, b, vb(:,:,1:5), DT, s2 )
    call Morph_LS_HJ_Step ( grid, ex, b, vb, -DT, s3 )
    call Morph_LS_HJ_Step ( grid, ex, b, vb, nan, s4 )
    call Morph_LS_HJ_Step ( none, ex, b, vb, DT, s5 )
    allocate( e3(0,0,0), f3(0,0,0), e2(0,0) )
    call Morph_LS_Grid_Scatter ( none, e2, e3, z1 )
    call Morph_LS_Grid_Gather ( none, e3, e2, z2 )
    call Morph_LS_HJ_Step ( none, ex, e3, f3, DT, z3 )
    allocate( ew(1-NG:NG, 1-NG:NG, 0) )
    call Morph_LS_Halo_Fill ( none, ex, ew, z4 )
    call check ( z1%code == MORPH_ERR_SHAPE .and. z2%code == MORPH_ERR_SHAPE .and. &
                 z3%code == MORPH_ERR_SHAPE .and. z4%code == MORPH_ERR_SHAPE, &
                 'a grid not made: scatter, gather, the ghosts and the step refuse arrays of size 0' )
    call Morph_LS_HJ_Step ( grid, ex, b(1:7,:,:), vb, DT, s7 )
    call Morph_LS_HJ_Step ( grid, ex, b(:,1:7,:), vb, DT, s8 )
    call Morph_LS_HJ_Step ( grid, ex, b, vb(1:7,:,:), DT, s9 )
    call Morph_LS_HJ_Step ( grid, ex, b, vb(:,1:7,:), DT, s10 )
    ok = s7%code == MORPH_ERR_SHAPE .and. s8%code == MORPH_ERR_SHAPE .and. &
         s9%code == MORPH_ERR_SHAPE .and. s10%code == MORPH_ERR_SHAPE
    do k = 1, 6
      if ( .not. same( b(:,:,k), b0(:,:,k) ) ) ok = .false.
    enddo
    call check ( ok .and. all( bad == 1.0_R8 ) .and. s1%code == MORPH_ERR_SHAPE .and. &
                 s2%code == MORPH_ERR_SHAPE .and. s3%code == MORPH_ERR_INPUT .and. &
                 s4%code == MORPH_ERR_INPUT .and. s5%code == MORPH_ERR_SHAPE, &
                 'step refuses phi or V of other shapes (each axis), dt < 0, dt = NaN, a grid not made; '// &
                 'phi untouched' )

    ! the grid of rank 0 of 2, stepped with the exchange of one rank
    call world%init ( 2, st )
    r0%world => world
    r0%me = 0
    call Morph_LS_Grid_Init ( g2, 8, [2,3], H, [0,0,0,1,1,1], r0, st )
    deallocate( b )
    allocate( b(8,8,3) )
    b = 2.0_R8
    call Morph_LS_HJ_Step ( g2, ex, b, vb(:,:,1:3), DT, s6 )
    call Morph_LS_HJ_Step ( g2, ex, b, 0.0_R8*vb(:,:,1:3), DT, s11 )
    call check ( st%ok() .and. s6%code == MORPH_ERR_EXCHANGE .and. s11%code == MORPH_ERR_EXCHANGE .and. &
                 all( b == 2.0_R8 ), &
                 'step refuses an exchange that is not the grid''s, V = 0 too; phi untouched' )
  end subroutine step_refusals


  ! ===========================================================================
  !> swap, allgather and ordered sum on n simulated ranks against the serial.
  subroutine sim_collectives ()
    integer(I4) :: n
    do n = 1, 4
      call sim_collectives_n ( n )
    enddo
  end subroutine sim_collectives


  subroutine sim_collectives_n ( n )
    integer(I4), intent(in) :: n
    type(morph_sim_world_t), target :: world
    type(morph_exchange_serial_t)   :: ser
    type(morph_status_t) :: st
    real(R8)    :: parts(7), serial_total, totals(0:n-1)
    integer(I4) :: owner(7), b, bad(0:n-1)
    character(len=96) :: name

    ! in block order 2.5; with these owners the loop over the blocks reversed,
    ! each rank's partials summed first, or the ranks' partials one rank after
    ! the other, give 3.5 / 3.0 / 4.5 (2 ranks), 4.5 / 3.0 / 1.0 (3 ranks),
    ! 4.5 / 4.0 / 5.0 (4 ranks): computed in Python, IEEE doubles
    parts = [ 1.0_R8, 2.0_R8**53, 1.0_R8, -(2.0_R8**53), 1.0_R8, 1.0_R8, 0.5_R8 ]
    select case ( n )
    case ( 1 ); owner = [0,0,0,0,0,0,0]
    case ( 2 ); owner = [0,1,1,0,0,0,1]
    case ( 3 ); owner = [0,0,1,2,2,1,0]
    case default; owner = [0,3,1,2,2,1,0]
    end select
    call ser%ordered_sum ( 7, [(0, b = 1, 7)], parts, serial_total, st )
    call world%init ( n, st )
    bad = 0
    !$omp parallel num_threads(n)
    block
      type(morph_exchange_sim_t) :: ex
      type(morph_status_t) :: s1, s2, s3
      real(R8), allocatable :: sb(:), rb(:), mine(:), all_of(:)
      integer(I4), allocatable :: sc(:), rc(:), gc(:)
      integer(I4) :: t, p, q, k, m, nerr
      t = 0
      !$ t = omp_get_thread_num()
      ex%me = t
      ex%world => world
      nerr = 0
      ! swap: t sends q the values 1000 t + 10 q + k, k = 1 .. c(t,q)
      allocate( sc(0:n-1), rc(0:n-1), gc(0:n-1) )
      do q = 0, n-1
        sc(q) = mod(t + 2*q, 3) + merge(1, 0, q == t)
        rc(q) = mod(q + 2*t, 3) + merge(1, 0, q == t)
      enddo
      allocate( sb(sum(sc)), rb(sum(rc)) )
      m = 0
      do q = 0, n-1
        do k = 1, sc(q)
          m = m + 1
          sb(m) = real(1000*t + 10*q + k, R8)
        enddo
      enddo
      call ex%swap ( sc, sb, rc, rb, s1 )
      m = 0
      do p = 0, n-1
        do k = 1, rc(p)
          m = m + 1
          if ( rb(m) /= real(1000*p + 10*t + k, R8) ) nerr = nerr + 1
        enddo
      enddo
      ! allgather: rank p gives mod(p, 3) values 100 (p+1) + k
      do p = 0, n-1
        gc(p) = mod(p, 3)
      enddo
      allocate( mine(gc(t)), all_of(sum(gc)) )
      do k = 1, gc(t)
        mine(k) = real(100*(t+1) + k, R8)
      enddo
      call ex%allgather ( mine, all_of, gc, s2 )
      m = 0
      do p = 0, n-1
        do k = 1, gc(p)
          m = m + 1
          if ( all_of(m) /= real(100*(p+1) + k, R8) ) nerr = nerr + 1
        enddo
      enddo
      ! ordered sum: the partials of t's blocks, in block order
      call ex%ordered_sum ( 7, owner, pack(parts, owner == t), totals(t), s3 )
      if ( .not. ( s1%ok() .and. s2%ok() .and. s3%ok() ) ) nerr = nerr + 1
      bad(t) = nerr
    end block
    !$omp end parallel
    write(name,'(A,I0,A,A)') 'simulated exchange, ', n, trim(merge(' rank: ', ' ranks:', n == 1)), &
                             ' swap and allgather in rank order, the serial sum''s bits'
    call check ( all( bad == 0 ) .and. all( transfer(totals, 0_I8, n) == transfer(serial_total, 0_I8) ) .and. &
                 serial_total == 2.5_R8, trim(name) )
  end subroutine sim_collectives_n


  !> Thousands of rounds of swaps and allgathers whose counts change from one
  !> round to the next, on 4 ranks: a missing barrier shows here.
  subroutine sim_stress ()
    integer(I4), parameter :: N = 4, ROUNDS = 3000
    type(morph_sim_world_t), target :: world
    type(morph_status_t) :: st
    integer(I4) :: bad(0:N-1)

    call world%init ( N, st )
    bad = 0
    !$omp parallel num_threads(N)
    block
      type(morph_exchange_sim_t) :: ex
      type(morph_status_t) :: s1, s2
      real(R8) :: sb(4*N), rb(4*N), mine(3), all_of(3*N), parts(7), tot
      integer(I4) :: sc(0:N-1), rc(0:N-1), gc(0:N-1), owner(7)
      integer(I4) :: t, p, q, k, m, r, nerr
      type(morph_status_t) :: s3
      t = 0
      !$ t = omp_get_thread_num()
      ex%me = t
      ex%world => world
      nerr = 0
      do r = 1, ROUNDS
        do q = 0, N-1
          sc(q) = mod(t + 2*q + r, 4)
          rc(q) = mod(q + 2*t + r, 4)
          gc(q) = mod(q + r, 3)
        enddo
        m = 0
        do q = 0, N-1
          do k = 1, sc(q)
            m = m + 1
            sb(m) = real(10000*r + 100*t + 10*q + k, R8)
          enddo
        enddo
        call ex%swap ( sc, sb, rc, rb, s1 )
        m = 0
        do p = 0, N-1
          do k = 1, rc(p)
            m = m + 1
            if ( rb(m) /= real(10000*r + 100*p + 10*t + k, R8) ) nerr = nerr + 1
          enddo
        enddo
        do k = 1, gc(t)
          mine(k) = real(10000*r + 100*t + k, R8)
        enddo
        call ex%allgather ( mine(1:gc(t)), all_of, gc, s2 )
        m = 0
        do p = 0, N-1
          do k = 1, gc(p)
            m = m + 1
            if ( all_of(m) /= real(10000*r + 100*p + k, R8) ) nerr = nerr + 1
          enddo
        enddo
        ! an ordered sum whose owners turn with the round: 2.5 in block order
        parts = [ 1.0_R8, 2.0_R8**53, 1.0_R8, -(2.0_R8**53), 1.0_R8, 1.0_R8, 0.5_R8 ]
        do k = 1, 7
          owner(k) = mod(k*k + r, N)
        enddo
        call ex%ordered_sum ( 7, owner, pack(parts, owner == t), tot, s3 )
        if ( tot /= 2.5_R8 ) nerr = nerr + 1
        if ( .not. ( s1%ok() .and. s2%ok() .and. s3%ok() ) ) nerr = nerr + 1
      enddo
      bad(t) = nerr
    end block
    !$omp end parallel
    call check ( all( bad == 0 ), &
                 'simulated exchange, 4 ranks: 3000 rounds of a swap, an allgather and an ordered sum, every value' )
  end subroutine sim_stress


  !> Refusals that every rank takes: none is left waiting at a barrier. Each
  !> world is first used once by all its ranks, so that the mailbox of a
  !> missing rank holds a valid post from before and only the team check
  !> refuses; each refusal below is the only check of its call that fails.
  subroutine sim_refusals ()
    integer(I4), parameter :: NK = 20
    type(morph_sim_world_t), target :: world, three
    type(morph_exchange_sim_t) :: out, lone
    type(morph_status_t) :: st, s1, s2, s3
    real(R8) :: r(4), q(4), tot
    integer(I4) :: codes(0:1,NK), setup(0:2), k

    call world%init ( 2, st )
    call three%init ( 3, s1 )
    setup = -1
    do k = 2, 3
      !$omp parallel num_threads(k)
      block
        type(morph_exchange_sim_t) :: ex
        type(morph_status_t) :: a
        real(R8) :: sb(1), rb(1)
        integer(I4) :: t, c(0:2), p
        t = 0
        !$ t = omp_get_thread_num()
        ex%me = t
        if ( k == 2 ) ex%world => world
        if ( k == 3 ) ex%world => three
        do p = 0, k-1
          c(p) = merge(1, 0, p == t)
        enddo
        sb = real(t, R8)
        call ex%swap ( c(0:k-1), sb, c(0:k-1), rb, a )
        if ( a%ok() .and. rb(1) == real(t, R8) .and. k == 3 ) setup(t) = 0
        if ( .not. ( a%ok() .and. rb(1) == real(t, R8) ) ) setup(t) = 1
      end block
      !$omp end parallel
    enddo
    call check ( st%ok() .and. s1%ok() .and. all( setup == 0 ), &
                 'simulated exchange: each world used once by all its ranks' )

    ! outside the parallel region: a team of one thread for two ranks
    out%world => world
    r = 0.0_R8
    q = 0.0_R8
    call out%swap ( [0,0], q, [0,0], r, s1 )
    call out%allgather ( q(1:0), r, [0,0], s2 )
    call out%ordered_sum ( 0, [integer(I4) ::], [real(R8) ::], tot, s3 )
    call check ( s1%code == MORPH_ERR_EXCHANGE .and. s2%code == MORPH_ERR_EXCHANGE .and. &
                 s3%code == MORPH_ERR_EXCHANGE, &
                 'simulated exchange: called outside the ranks'' region, refused' )
    call lone%swap ( [0], q, [0], r, s1 )
    call lone%allgather ( q(1:0), r, [0], s2 )
    call lone%ordered_sum ( 0, [integer(I4) ::], [real(R8) ::], tot, s3 )
    call check ( s1%code == MORPH_ERR_EXCHANGE .and. s2%code == MORPH_ERR_EXCHANGE .and. &
                 s3%code == MORPH_ERR_EXCHANGE, &
                 'simulated exchange: an instance without a world, refused' )

    codes = -1
    !$omp parallel num_threads(2)
    block
      type(morph_exchange_sim_t) :: ex
      type(morph_status_t) :: c(NK)
      real(R8) :: buf(4), got(4), big(8), tot1
      integer(I4) :: t, u, i
      buf = 0.0_R8
      big = 0.0_R8
      t = 0
      !$ t = omp_get_thread_num()
      u = 1 - t
      ! 1: a team of 2 threads for a world of 3 ranks
      ex%world => three
      ex%me = t
      call ex%swap ( [0,0,0], buf, [0,0,0], got, c(1) )
      ex%world => world
      ! 2-4: two threads that swap their ranks
      ex%me = u
      call ex%swap ( [0,0], buf, [0,0], got, c(2) )
      call ex%allgather ( buf(1:1), got, [1,1], c(3) )
      call ex%ordered_sum ( 2, [0,1], buf(1:1), tot1, c(4) )
      ex%me = t
      ! 5: rank 0 sends rank 1 two values, rank 1 expects one
      call ex%swap ( [0, 2*u], buf, [t, 0], got, c(5) )
      ! 6: a negative count that both ranks agree on
      call ex%swap ( [0, -u], buf, [-t, 0], got, c(6) )
      ! 7: rank 0 sends more than its buffer holds
      if ( t == 0 ) then
        call ex%swap ( [0,5], buf, [0,0], big, c(7) )
      else
        call ex%swap ( [0,0], buf, [5,0], big, c(7) )
      endif
      ! 8: rank 1 receives more than its buffer holds
      if ( t == 0 ) then
        call ex%swap ( [0,5], big, [0,0], got, c(8) )
      else
        call ex%swap ( [0,0], big, [5,0], got, c(8) )
      endif
      ! 9: counts for three ranks, on both ranks
      call ex%swap ( [0,0,0], buf, [0,0,0], got, c(9) )
      ! 10: allgather counts that differ between the ranks
      call ex%allgather ( buf(1:1+t), got, [1, 1+t], c(10) )
      ! 11: rank 1's count is not its values
      call ex%allgather ( buf(1:1), got, [1,2], c(11) )
      ! 12: a gather larger than rank 1's buffer
      if ( t == 0 ) then
        call ex%allgather ( buf(1:3), big, [3,3], c(12) )
      else
        call ex%allgather ( buf(1:3), got, [3,3], c(12) )
      endif
      ! 13: a block on rank 2 of 2, on both ranks
      if ( t == 0 ) then
        call ex%ordered_sum ( 2, [0,2], buf(1:1), tot1, c(13) )
      else
        call ex%ordered_sum ( 2, [0,2], buf(1:0), tot1, c(13) )
      endif
      ! 14: rank 1 gives two partials for one block
      call ex%ordered_sum ( 2, [0,1], buf(1:1+t), tot1, c(14) )
      ! 15: owners of other lengths, each consistent with its own partials
      if ( t == 0 ) then
        call ex%ordered_sum ( 2, [0,1], buf(1:1), tot1, c(15) )
      else
        call ex%ordered_sum ( 3, [0,1,1], buf(1:2), tot1, c(15) )
      endif
      ! 16: rank 1 has owners for 2 blocks of 3
      call ex%ordered_sum ( 2+t, [0,1], buf(1:1), tot1, c(16) )
      ! 17: owners that differ: [0,1] on rank 0, [1,0] on rank 1
      call ex%ordered_sum ( 2, [t, u], buf(1:1), tot1, c(17) )
      ! 18: allgather counts for three ranks, on both ranks
      call ex%allgather ( buf(1:1), got, [1,1,0], c(18) )
      ! 20: owners of other lengths, the longer on rank 0
      if ( t == 0 ) then
        call ex%ordered_sum ( 3, [0,1,1], buf(1:1), tot1, c(20) )
      else
        call ex%ordered_sum ( 2, [0,1], buf(1:1), tot1, c(20) )
      endif
      ! 19: a block on rank -1, on both ranks
      if ( t == 0 ) then
        call ex%ordered_sum ( 2, [0,-1], buf(1:1), tot1, c(19) )
      else
        call ex%ordered_sum ( 2, [0,-1], buf(1:0), tot1, c(19) )
      endif
      do i = 1, NK
        codes(t,i) = c(i)%code
      enddo
    end block
    !$omp end parallel
    call check ( all( codes(:,1) == MORPH_ERR_EXCHANGE ), &
                 'simulated exchange: a team of 2 threads for 3 ranks, refused on both' )
    call check ( all( codes(:,2:4) == MORPH_ERR_EXCHANGE ), &
                 'simulated exchange: ranks that are not their threads, refused on both (swap, allgather, sum)' )
    call check ( all( codes(:,5:9) == MORPH_ERR_EXCHANGE ), &
                 'simulated swap: counts the ranks disagree on, a negative one, more than a buffer holds, '// &
                 'counts for three ranks of two: refused on both' )
    call check ( all( codes(:,10:12) == MORPH_ERR_EXCHANGE ) .and. all( codes(:,18) == MORPH_ERR_EXCHANGE ), &
                 'simulated allgather: counts that differ, a count not the rank''s values, more than the '// &
                 'buffer, counts for three ranks of two: refused on both' )
    call check ( all( codes(:,13:17) == MORPH_ERR_EXCHANGE ) .and. all( codes(:,19:20) == MORPH_ERR_EXCHANGE ), &
                 'simulated ordered sum: a block on rank 2 or -1 of 2, partials or owners not one per block, '// &
                 'owners of other lengths or values: refused on both' )
  end subroutine sim_refusals


  ! ===========================================================================
  !> n simulated ranks, the blocks split in several ways: the bits of the
  !> whole grid, and the ghosts of each rank after one fill.
  subroutine steps_ranks ()
    integer(I4) :: own20(20), own12(12), k, ib, jb

    do k = 1, 20
      own20(k) = 0
    enddo
    call ranks_case ( 1, 8, [5,4], own20, -1, 'step, 1 simulated rank: 5 x 4 blocks of 8' )
    do k = 1, 20
      ib = mod(k-1, 5) + 1
      jb = (k-1)/5 + 1
      own20(k) = mod(ib + jb, 2)
    enddo
    call ranks_case ( 2, 8, [5,4], own20, -1, 'step, 2 ranks, a checkerboard of 5 x 4 blocks of 8' )
    do k = 1, 20
      own20(k) = (k-1)/10
    enddo
    call ranks_case ( 2, 8, [5,4], own20, 0, 'step, 2 ranks, halves, V = 0 on rank 0''s blocks' )
    do k = 1, 20
      own20(k) = mod(k-1, 3)
    enddo
    call ranks_case ( 3, 8, [5,4], own20, -1, 'step, 3 ranks, blocks dealt in turn' )
    do k = 1, 20
      own20(k) = 2*mod(k-1, 2)
    enddo
    call ranks_case ( 3, 8, [5,4], own20, -1, 'step, 3 ranks, rank 1 without blocks' )
    do k = 1, 20
      own20(k) = (k-1)/5
    enddo
    call ranks_case ( 4, 8, [5,4], own20, -1, 'step, 4 ranks, one row of blocks each' )
    do k = 1, 12
      own12(k) = mod(5*k, 4)
    enddo
    call ranks_case ( 4, 3, [4,3], own12, -1, 'step, 4 ranks, 4 x 3 blocks of 3, scattered' )
    call zero_velocity_ranks ()
  end subroutine steps_ranks


  !> vzero >= 0: V = 0 on the blocks of that rank.
  subroutine ranks_case ( n, nbn, nblk, owner, vzero, name )
    integer(I4),      intent(in) :: n, nbn, nblk(2), owner(:), vzero
    character(len=*), intent(in) :: name
    type(morph_sim_world_t), target :: world
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:), v0(:,:), ref(:,:), out(:,:), wf(:,:)
    integer(I4) :: codes(0:n-1), b, i0, j0
    logical :: halos(0:n-1), moved

    call fields ( nbn*nblk(1), nbn*nblk(2), phi0, v0 )
    if ( vzero >= 0 ) then
      do b = 1, size(owner)
        if ( owner(b) /= vzero ) cycle
        i0 = mod(b-1, nblk(1)) * nbn
        j0 = ((b-1)/nblk(1)) * nbn
        v0(i0+1:i0+nbn, j0+1:j0+nbn) = 0.0_R8
      enddo
    endif
    ref = reference ( phi0, v0 )
    wf  = with_ghosts ( phi0 )
    out = SENTINEL + 0.0_R8 * phi0
    call world%init ( n, st )
    codes = -1
    halos = .false.
    !$omp parallel num_threads(n)
    block
      type(morph_exchange_sim_t) :: ex
      type(morph_status_t) :: s
      integer(I4) :: t
      logical :: hk
      t = 0
      !$ t = omp_get_thread_num()
      ex%me = t
      ex%world => world
      call advance ( ex, nbn, nblk, owner, phi0, v0, wf, out, hk, s )
      codes(t) = s%code
      halos(t) = hk
    end block
    !$omp end parallel
    moved = .not. same(ref, phi0)
    call check ( st%ok() .and. all( codes == MORPH_OK ) .and. all( halos ) .and. same(out, ref) .and. moved, &
                 name//': the whole grid''s bits' )
  end subroutine ranks_case


  !> V = 0 everywhere on 2 ranks: phi untouched on both, a -0 too.
  subroutine zero_velocity_ranks ()
    type(morph_sim_world_t), target :: world
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:), v0(:,:), out(:,:)
    integer(I4) :: codes(0:1)

    call fields ( 16, 24, phi0, v0 )
    phi0(12,20) = -0.0_R8
    v0 = 0.0_R8
    v0(12,20) = -0.0_R8
    out = SENTINEL + 0.0_R8 * phi0
    call world%init ( 2, st )
    codes = -1
    !$omp parallel num_threads(2)
    block
      type(morph_exchange_sim_t) :: ex
      type(morph_ls_grid_t) :: grid
      type(morph_status_t) :: s
      real(R8), allocatable :: b(:,:,:), vb(:,:,:)
      integer(I4) :: t
      t = 0
      !$ t = omp_get_thread_num()
      ex%me = t
      ex%world => world
      call Morph_LS_Grid_Init ( grid, 8, [2,3], H, [0,1,0,1,0,1], ex, s )
      allocate( b(8,8,3), vb(8,8,3) )
      call Morph_LS_Grid_Scatter ( grid, phi0, b, s )
      call Morph_LS_Grid_Scatter ( grid, v0, vb, s )
      call Morph_LS_HJ_Step ( grid, ex, b, vb, DT, s )
      codes(t) = s%code
      call Morph_LS_Grid_Gather ( grid, b, out, s )
    end block
    !$omp end parallel
    call check ( st%ok() .and. all( codes == MORPH_OK ) .and. same(out, phi0), &
                 'step, 2 ranks: V = 0 everywhere leaves phi untouched, its -0 too' )
  end subroutine zero_velocity_ranks


  !> A refusal on one rank is taken by both, phi untouched on both.
  subroutine refused_on_one_rank ()
    type(morph_sim_world_t), target :: world
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:), v0(:,:)
    integer(I4) :: codes(0:1,3)
    logical :: kept(0:1), named(0:1)

    call fields ( 16, 24, phi0, v0 )
    call world%init ( 2, st )
    codes = -1
    kept  = .false.
    named = .false.
    !$omp parallel num_threads(2)
    block
      type(morph_exchange_sim_t) :: ex, oth
      type(morph_ls_grid_t) :: grid, ogrid
      type(morph_status_t) :: s, a, b, c
      real(R8), allocatable :: phi(:,:,:), vb(:,:,:), phi_in(:,:,:)
      real(R8) :: dt_t
      integer(I4) :: t
      t = 0
      !$ t = omp_get_thread_num()
      ex%me = t
      ex%world => world
      call Morph_LS_Grid_Init ( grid, 8, [2,3], H, [0,0,0,1,1,1], ex, s )
      allocate( phi(8,8,3), vb(8,8,3) )
      call Morph_LS_Grid_Scatter ( grid, phi0, phi, s )
      call Morph_LS_Grid_Scatter ( grid, v0, vb, s )
      phi_in = phi
      ! rank 1 passes V of another shape
      if ( t == 1 ) then
        call Morph_LS_HJ_Step ( grid, ex, phi, vb(:,:,1:2), DT, a )
      else
        call Morph_LS_HJ_Step ( grid, ex, phi, vb, DT, a )
      endif
      ! rank 0 passes dt < 0
      dt_t = DT
      if ( t == 0 ) dt_t = -DT
      call Morph_LS_HJ_Step ( grid, ex, phi, vb, dt_t, b )
      ! each rank steps with the other rank's grid
      oth%me = 1 - t
      oth%world => world
      call Morph_LS_Grid_Init ( ogrid, 8, [2,3], H, [0,0,0,1,1,1], oth, s )
      call Morph_LS_HJ_Step ( ogrid, ex, phi, vb, DT, c )
      codes(t,:) = [ a%code, b%code, c%code ]
      kept(t) = all( transfer(phi, 0_I8, size(phi)) == transfer(phi_in, 0_I8, size(phi_in)) )
      if ( t == 0 ) named(t) = index(a%message, 'refused on rank 1') > 0 .and. index(b%message, 'dt =') > 0
      if ( t == 1 ) named(t) = index(a%message, 'V 8 x 8 x 2') > 0 .and. index(b%message, 'refused on rank 0') > 0
    end block
    !$omp end parallel
    call check ( all( codes(:,1) == MORPH_ERR_SHAPE ) .and. all( codes(:,2) == MORPH_ERR_INPUT ) .and. &
                 all( kept ) .and. all( named ), &
                 'step, 2 ranks: a refusal on one rank is taken by both, phi untouched on both' )
    call check ( all( codes(:,3) == MORPH_ERR_EXCHANGE ), &
                 'step, 2 ranks: each with the other rank''s grid, refused on both' )
  end subroutine refused_on_one_rank


  !> Three ranks, ranks 1 and 2 refusing with different codes: each refusing
  !> rank reports its own, rank 0 the first in rank order.
  subroutine refused_on_two_ranks ()
    type(morph_sim_world_t), target :: world
    type(morph_status_t)  :: st
    real(R8), allocatable :: phi0(:,:), v0(:,:)
    integer(I4) :: codes(0:2)
    logical :: named(0:2)

    call fields ( 24, 16, phi0, v0 )
    call world%init ( 3, st )
    codes = -1
    named = .false.
    !$omp parallel num_threads(3)
    block
      type(morph_exchange_sim_t) :: ex
      type(morph_ls_grid_t) :: grid
      type(morph_status_t) :: s, a
      real(R8), allocatable :: phi(:,:,:), vb(:,:,:)
      real(R8) :: dt_t
      integer(I4) :: t, nv
      t = 0
      !$ t = omp_get_thread_num()
      ex%me = t
      ex%world => world
      call Morph_LS_Grid_Init ( grid, 8, [3,2], H, [0,1,2,0,1,2], ex, s )
      allocate( phi(8,8,2), vb(8,8,2) )
      call Morph_LS_Grid_Scatter ( grid, phi0, phi, s )
      call Morph_LS_Grid_Scatter ( grid, v0, vb, s )
      dt_t = DT
      if ( t == 1 ) dt_t = -DT
      nv = 2
      if ( t == 2 ) nv = 1
      call Morph_LS_HJ_Step ( grid, ex, phi, vb(:,:,1:nv), dt_t, a )
      codes(t) = a%code
      if ( t == 0 ) named(t) = index(a%message, 'refused on rank 1') > 0
      if ( t == 1 ) named(t) = index(a%message, 'dt =') > 0
      if ( t == 2 ) named(t) = index(a%message, 'V 8 x 8 x 1') > 0
    end block
    !$omp end parallel
    call check ( st%ok() .and. all( codes == [MORPH_ERR_INPUT, MORPH_ERR_INPUT, MORPH_ERR_SHAPE] ) .and. &
                 all( named ), &
                 'step, 3 ranks, two refusing: each its own code, the third the first refusing rank''s' )
  end subroutine refused_on_two_ranks

end program test_morph_ls_blocks
