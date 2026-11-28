!>@brief The safety reinitialisation of a level-set field toward a signed
!> distance, and the measure that decides when it is needed, on a whole 2-D
!> grid of nodes x0 + (i,j)*h and on the grid of blocks of Morph_LS_Grid.
!>
!> THE MEASURE. E = mean of | |grad phi| - 1 | over the nodes with
!> |phi| <= band, the gradient by central differences (one-sided on the edge
!> nodes). The sum runs in a fixed node order, so E does not depend on the
!> number of threads. Morph_LS_Reinit_If_Needed_Full reinitialises only when
!> E exceeds a threshold: a reinitialisation moves the zero level by a small
!> amount every time it runs, so it must not run at every step.
!>
!> THE REINITIALISATION. iters pseudo-time steps of length 0.45 h of
!>   phi_tau = S(phi0) (1 - |grad phi|),  S = phi0 / sqrt(phi0^2 + |grad phi0|^2 h^2)
!> (Sussman, Smereka and Osher 1994), |grad phi| by Godunov's scheme with the
!> WENO-Z derivatives of Morph_LS_HJ (the branch S >= 0), with two corrections
!> at the nodes whose four neighbours include a sign change of phi0 (the cut
!> nodes; the edge nodes are never cut):
!>   * the subcell fix (Russo and Smereka 2000): there the update is
!>     -(phi - D)/h, D = phi0/|grad phi0| with |grad phi0| from the largest of
!>     the central and the two one-sided differences of each axis, so that the
!>     zero level stays where phi0 put it;
!>   * a correction of the residual displacement of the zero level after the
!>     step (Hartmann, Meinke and Schroeder 2008), weight beta: for each
!>     neighbour across a sign change, the zero of the linear interpolant of
!>     the new values at the position of the old zero, averaged and scaled by
!>     |grad phi|.
!> The edge nodes take the plain update. The order of operations is fixed (the
!> parentheses), as in Morph_LS_HJ: a build without fused multiply-adds gives
!> the bits of any implementation of the same formulas in the same order.
!>
!> THE BLOCKS. The computations at the nodes are kernels on one block with
!> its ghost layers, w(1-NG:nx+NG, 1-NG:ny+NG), and the flags edge(4) of its
!> sides W, E, S, N (as in Morph_LS_Grid): a side is an edge when it faces the
!> edge of the grid or an inactive block, and there the gradients are
!> one-sided and the nodes of the side are never cut, as on the edges of the
!> whole grid; across the other sides the ghosts hold the neighbour's nodes.
!> The routines on the whole grid call them on one block with four edges, the
!> routines on the grid of blocks on each active block, the ghosts filled by
!> Morph_LS_Halo_Fill without V (linear toward an inactive block, as at the
!> edge of the grid). So, with every block active, Morph_LS_Reinit_Grid gives
!> the bits of Morph_LS_Reinit_Full for any blocks, any split among ranks and
!> any number of threads. On a band, a node's value after iters steps depends
!> only on the nodes within L1 distance 4 iters (WENO-Z reads 3 nodes along
!> each axis, the correction one more), so the nodes farther than that from
!> every inactive node get the bits of the whole grid started from the same
!> phi. The measure on blocks sums each block in node order, then the blocks
!> in block order (ordered_sum of Morph_Exchange): the same bits on every rank
!> and for any number of ranks and threads; on a single block, those of the
!> whole grid; otherwise the order of the additions differs from the whole
!> grid's, which changes the last bits.
!>
!> Never stops, no I/O: inconsistent arguments return MORPH_ERR_SHAPE,
!> MORPH_ERR_INPUT or MORPH_ERR_EXCHANGE in the status, phi untouched.
module Morph_LS_Reinit
  use iso_fortran_env, only: I4 => int32, I8 => int64, R8 => real64
  use Morph_Types_m
  use Morph_Exchange,  only: morph_exchange_t
  use Morph_LS_HJ,     only: MORPH_LS_NG, Morph_LS_Ghosts_Extrap_All, Morph_LS_Grad_Godunov
  use Morph_LS_Grid,   only: morph_ls_grid_t, Morph_LS_Halo_Fill

  implicit none
  private

  integer(I4), parameter :: NG = MORPH_LS_NG
  logical,     parameter :: ALL_EDGES(4) = .true.

  public :: Morph_LS_Eikonal_Error, Morph_LS_Reinit_Full, Morph_LS_Reinit_If_Needed_Full
  public :: Morph_LS_Eikonal_Terms, Morph_LS_Reinit_Setup, Morph_LS_Reinit_Advance, Morph_LS_Reinit_Correct
  public :: Morph_LS_Eikonal_Error_Grid, Morph_LS_Reinit_Grid, Morph_LS_Reinit_If_Needed_Grid

contains

  !> The gradient at the nodes of a block f(1-NG:nx+NG, 1-NG:ny+NG):
  !> (f(i+1) - f(i-1))/(2h) where both neighbours are there (across a side
  !> that is not an edge, a ghost), (f(2) - f(1))/h and (f(n) - f(n-1))/h on
  !> the nodes of an edge side. At least 2 nodes along each axis.
  subroutine block_gradient ( f, h, edge, gx, gy )
    real(R8),    intent(in)  :: f(1-NG:, 1-NG:)
    real(R8),    intent(in)  :: h
    logical,     intent(in)  :: edge(4)
    real(R8),    intent(out) :: gx(:,:), gy(:,:)
    integer(I4) :: nx, ny, i, j

    nx = size(f,1) - 2*NG
    ny = size(f,2) - 2*NG
    !$omp parallel do private(i,j)
    do j = 1, ny
      do i = 1, nx
        if ( i == 1 .and. edge(1) ) then
          gx(i,j) = (f(2,j) - f(1,j)) / h
        else if ( i == nx .and. edge(2) ) then
          gx(i,j) = (f(nx,j) - f(nx-1,j)) / h
        else
          gx(i,j) = (f(i+1,j) - f(i-1,j)) / (2.0_R8*h)
        endif
        if ( j == 1 .and. edge(3) ) then
          gy(i,j) = (f(i,2) - f(i,1)) / h
        else if ( j == ny .and. edge(4) ) then
          gy(i,j) = (f(i,ny) - f(i,ny-1)) / h
        else
          gy(i,j) = (f(i,j+1) - f(i,j-1)) / (2.0_R8*h)
        endif
      enddo
    enddo
    !$omp end parallel do
  end subroutine block_gradient


  !> The measure's part from one block w(1-NG:nx+NG, 1-NG:ny+NG): the sum s of
  !> | |grad phi| - 1 | over its nodes with |phi| <= band, in node order (i
  !> fastest), and their number n. The ghosts across the sides that are not
  !> edges must be filled.
  subroutine Morph_LS_Eikonal_Terms ( w, h, band, edge, s, n )
    real(R8),    intent(in)  :: w(1-NG:, 1-NG:)
    real(R8),    intent(in)  :: h, band
    logical,     intent(in)  :: edge(4)
    real(R8),    intent(out) :: s
    integer(I4), intent(out) :: n
    real(R8), allocatable :: gx(:,:), gy(:,:)
    integer(I4) :: nx, ny, i, j

    nx = size(w,1) - 2*NG
    ny = size(w,2) - 2*NG
    allocate( gx(nx,ny), gy(nx,ny) )
    call block_gradient ( w, h, edge, gx, gy )
    s = 0.0_R8
    n = 0
    do j = 1, ny
      do i = 1, nx
        if ( abs(w(i,j)) <= band ) then
          s = s + abs( sqrt(gx(i,j)*gx(i,j) + gy(i,j)*gy(i,j)) - 1.0_R8 )
          n = n + 1
        endif
      enddo
    enddo
  end subroutine Morph_LS_Eikonal_Terms


  !> E = mean | |grad phi| - 1 | over the nodes with |phi| <= band; 0 when no
  !> node is in the band. Summed in node order (i fastest).
  real(R8) function Morph_LS_Eikonal_Error ( phi, h, band ) result ( e )
    real(R8), intent(in) :: phi(:,:)
    real(R8), intent(in) :: h, band
    real(R8), allocatable :: w(:,:)
    real(R8)    :: s
    integer(I4) :: n

    e = 0.0_R8
    if ( size(phi,1) < 2 .or. size(phi,2) < 2 ) return
    allocate( w(1-NG:size(phi,1)+NG, 1-NG:size(phi,2)+NG) )
    w = 0.0_R8
    w(1:size(phi,1), 1:size(phi,2)) = phi
    call Morph_LS_Eikonal_Terms ( w, h, band, ALL_EDGES, s, n )
    if ( n > 0 ) e = s / real(n,R8)
  end function Morph_LS_Eikonal_Error


  !> From phi0 on one block with its ghosts, p0w(1-NG:nx+NG, 1-NG:ny+NG): S at
  !> the nodes, its branch pos = S >= 0, the cut nodes and the distance dd of
  !> the subcell fix. The nodes of an edge side are never cut (dd = 0 there).
  subroutine Morph_LS_Reinit_Setup ( p0w, h, edge, s, pos, cut, dd )
    real(R8), intent(in)  :: p0w(1-NG:, 1-NG:)
    real(R8), intent(in)  :: h
    logical,  intent(in)  :: edge(4)
    real(R8), intent(out) :: s(:,:), dd(:,:)
    logical,  intent(out) :: pos(:,:), cut(:,:)
    real(R8), allocatable :: gx(:,:), gy(:,:)
    integer(I4) :: nx, ny, i, j
    real(R8)    :: gdx, gdy, f, b, a

    nx = size(p0w,1) - 2*NG
    ny = size(p0w,2) - 2*NG
    allocate( gx(nx,ny), gy(nx,ny) )

    ! The sign of phi0, smoothed by the grid: S
    call block_gradient ( p0w, h, edge, gx, gy )
    !$omp parallel do private(j)
    do j = 1, ny
      s(:,j) = p0w(1:nx,j) / sqrt( (p0w(1:nx,j)*p0w(1:nx,j) + ((gx(:,j)*gx(:,j) + gy(:,j)*gy(:,j))*h)*h) &
               + 1.0e-30_R8 )
    enddo
    !$omp end parallel do
    pos = s >= 0.0_R8

    ! The cut nodes and the distance D of the subcell fix (never on an edge)
    cut = .false.
    dd  = 0.0_R8
    !$omp parallel do private(i,j,f,b,a,gdx,gdy)
    do j = 1, ny
      do i = 1, nx
        if ( (i == 1 .and. edge(1)) .or. (i == nx .and. edge(2)) .or. &
             (j == 1 .and. edge(3)) .or. (j == ny .and. edge(4)) ) cycle
        a = p0w(i,j)
        cut(i,j) = a*p0w(i-1,j) <= 0.0_R8 .or. a*p0w(i+1,j) <= 0.0_R8 .or. &
                   a*p0w(i,j-1) <= 0.0_R8 .or. a*p0w(i,j+1) <= 0.0_R8
        f = p0w(i+1,j)
        b = p0w(i-1,j)
        gdx = max( max( abs(f - b)/(2.0_R8*h), abs(f - a)/h ), abs(a - b)/h )
        f = p0w(i,j+1)
        b = p0w(i,j-1)
        gdy = max( max( abs(f - b)/(2.0_R8*h), abs(f - a)/h ), abs(a - b)/h )
        dd(i,j) = a / sqrt( (gdx*gdx + gdy*gdy) + 1.0e-30_R8 )
      enddo
    enddo
    !$omp end parallel do
  end subroutine Morph_LS_Reinit_Setup


  !> One pseudo-time step of one block without the correction: from phi with
  !> NG layers of ghosts, w(1-NG:nx+NG, 1-NG:ny+NG), st = phi + 0.45 h R, R
  !> the plain update, or the subcell fix at the cut nodes.
  subroutine Morph_LS_Reinit_Advance ( w, pos, s, cut, dd, h, st )
    real(R8), intent(in)  :: w(1-NG:, 1-NG:)
    logical,  intent(in)  :: pos(:,:), cut(:,:)
    real(R8), intent(in)  :: s(:,:), dd(:,:)
    real(R8), intent(in)  :: h
    real(R8), intent(out) :: st(:,:)
    real(R8), allocatable :: g(:,:)
    integer(I4) :: nx, ny, i, j
    real(R8)    :: r

    nx = size(w,1) - 2*NG
    ny = size(w,2) - 2*NG
    allocate( g(nx,ny) )
    call Morph_LS_Grad_Godunov ( w, pos, h, g )
    !$omp parallel do private(i,j,r)
    do j = 1, ny
      do i = 1, nx
        if ( cut(i,j) ) then
          r = -(w(i,j) - dd(i,j)) / h
        else
          r = -s(i,j) * (g(i,j) - 1.0_R8)
        endif
        st(i,j) = w(i,j) + (0.45_R8*h)*r
      enddo
    enddo
    !$omp end parallel do
  end subroutine Morph_LS_Reinit_Advance


  !> The correction of the zero level at the cut nodes of one block, from the
  !> step st and phi0 with their ghosts, stw and p0w(1-NG:nx+NG, 1-NG:ny+NG):
  !> corr = st elsewhere.
  subroutine Morph_LS_Reinit_Correct ( stw, p0w, cut, h, beta, edge, corr )
    real(R8), intent(in)  :: stw(1-NG:, 1-NG:), p0w(1-NG:, 1-NG:)
    logical,  intent(in)  :: cut(:,:)
    real(R8), intent(in)  :: h, beta
    logical,  intent(in)  :: edge(4)
    real(R8), intent(out) :: corr(:,:)
    real(R8), allocatable :: gx(:,:), gy(:,:)
    integer(I4) :: nx, ny, i, j
    real(R8)    :: a, fs, nc, pj, sj, gn

    nx = size(stw,1) - 2*NG
    ny = size(stw,2) - 2*NG
    allocate( gx(nx,ny), gy(nx,ny) )
    call block_gradient ( stw, h, edge, gx, gy )
    !$omp parallel do private(i,j,a,fs,nc,pj,sj,gn)
    do j = 1, ny
      do i = 1, nx
        corr(i,j) = stw(i,j)
        if ( .not. cut(i,j) ) cycle
        a  = p0w(i,j)
        fs = 0.0_R8
        nc = 0.0_R8
        pj = p0w(i-1,j); sj = stw(i-1,j)
        call add_crossing ( a, pj, stw(i,j), sj, fs, nc )
        pj = p0w(i+1,j); sj = stw(i+1,j)
        call add_crossing ( a, pj, stw(i,j), sj, fs, nc )
        pj = p0w(i,j-1); sj = stw(i,j-1)
        call add_crossing ( a, pj, stw(i,j), sj, fs, nc )
        pj = p0w(i,j+1); sj = stw(i,j+1)
        call add_crossing ( a, pj, stw(i,j), sj, fs, nc )
        if ( nc > 0.0_R8 ) then
          gn = sqrt( (gx(i,j)**2 + gy(i,j)**2) + 1.0e-30_R8 )
          corr(i,j) = stw(i,j) - ((beta*gn)*fs) / max(nc, 1.0_R8)
        endif
      enddo
    enddo
    !$omp end parallel do
  end subroutine Morph_LS_Reinit_Correct


  !> iters steps of the reinitialisation of the header, on phi(nx,ny) in place.
  !> Refused, phi untouched: fewer than 3 nodes along an axis
  !> (MORPH_ERR_SHAPE); h <= 0, iters < 0 or beta < 0 (MORPH_ERR_INPUT).
  subroutine Morph_LS_Reinit_Full ( phi, h, iters, beta, status )
    real(R8),             intent(inout) :: phi(:,:)
    real(R8),             intent(in)    :: h
    integer(I4),          intent(in)    :: iters
    real(R8),             intent(in)    :: beta
    type(morph_status_t), intent(out)   :: status
    real(R8), allocatable :: p0w(:,:), w(:,:), stw(:,:), s(:,:), dd(:,:), st(:,:), corr(:,:)
    logical,  allocatable :: cut(:,:), pos(:,:)
    integer(I4) :: nx, ny, it
    character(len=MORPH_MSGLEN) :: msg

    call status%clear()
    nx = size(phi,1)
    ny = size(phi,2)
    if ( nx < 3 .or. ny < 3 ) then
      write(msg,'(A,I0,A,I0,A)') 'level set reinitialisation: phi is ', nx, ' x ', ny, ' (at least 3 x 3)'
      call status%fail ( MORPH_ERR_SHAPE, msg )
      return
    endif
    if ( .not. (h > 0.0_R8) .or. iters < 0 .or. .not. (beta >= 0.0_R8) ) then
      write(msg,'(A,ES10.3,A,I0,A,ES10.3,A)') 'level set reinitialisation: h = ', h, ', iters = ', iters, &
            ', beta = ', beta, ' (h > 0, iters >= 0, beta >= 0)'
      call status%fail ( MORPH_ERR_INPUT, msg )
      return
    endif

    allocate( s(nx,ny), dd(nx,ny), st(nx,ny), corr(nx,ny), cut(nx,ny), pos(nx,ny) )
    allocate( p0w(1-NG:nx+NG, 1-NG:ny+NG), w(1-NG:nx+NG, 1-NG:ny+NG), stw(1-NG:nx+NG, 1-NG:ny+NG) )
    p0w = 0.0_R8
    w   = 0.0_R8
    stw = 0.0_R8
    p0w(1:nx,1:ny) = phi
    call Morph_LS_Reinit_Setup ( p0w, h, ALL_EDGES, s, pos, cut, dd )

    do it = 1, iters
      ! The plain update everywhere, the subcell fix at the cut nodes
      w(1:nx,1:ny) = phi
      call Morph_LS_Ghosts_Extrap_All ( w )
      call Morph_LS_Reinit_Advance ( w, pos, s, cut, dd, h, st )
      ! The correction of the zero level at the cut nodes, from the values of
      ! st before any correction
      if ( beta > 0.0_R8 ) then
        stw(1:nx,1:ny) = st
        call Morph_LS_Reinit_Correct ( stw, p0w, cut, h, beta, ALL_EDGES, corr )
        phi = corr
      else
        phi = st
      endif
    enddo
  end subroutine Morph_LS_Reinit_Full


  !> One neighbour of a cut node: across a sign change of phi0 (a*pj < 0) adds
  !> the value at the old zero of the line through the new values, and counts
  !> it. (Adding zero elsewhere would change nothing: fs starts at +0, and a
  !> sum that starts at +0 is never -0.)
  pure subroutine add_crossing ( a, pj, sa, sj, fs, nc )
    real(R8), intent(in)    :: a, pj, sa, sj
    real(R8), intent(inout) :: fs, nc
    if ( a*pj < 0.0_R8 ) then
      fs = fs + (pj*sa - a*sj) / (pj - a)
      nc = nc + 1.0_R8
    endif
  end subroutine add_crossing


  !> E over |phi| <= band; when E > threshold, Morph_LS_Reinit_Full(iters,
  !> beta) and done = .true.; otherwise phi untouched and done = .false.
  subroutine Morph_LS_Reinit_If_Needed_Full ( phi, h, band, threshold, iters, beta, done, e, status )
    real(R8),             intent(inout) :: phi(:,:)
    real(R8),             intent(in)    :: h, band, threshold
    integer(I4),          intent(in)    :: iters
    real(R8),             intent(in)    :: beta
    logical,              intent(out)   :: done
    real(R8),             intent(out)   :: e
    type(morph_status_t), intent(out)   :: status

    call status%clear()
    done = .false.
    e = Morph_LS_Eikonal_Error ( phi, h, band )
    if ( .not. (e > threshold) ) return
    call Morph_LS_Reinit_Full ( phi, h, iters, beta, status )
    done = status%ok()
  end subroutine Morph_LS_Reinit_If_Needed_Full


  !> The verdicts of every rank on its own arguments, shared (an allgather):
  !> code and msg are this rank's; extra(:) is compared between the ranks bit
  !> by bit, and a rank whose extra differs from rank 0's makes every rank
  !> refuse with MORPH_ERR_INPUT. Refused on every rank when one rank refuses (its
  !> own message there, the first refusing rank's code elsewhere). ok on
  !> return: go on.
  subroutine share_verdicts ( exch, what, code, msg, extra, status, ok )
    class(morph_exchange_t),    intent(inout) :: exch
    character(len=*),           intent(in)    :: what
    integer(I4),                intent(in)    :: code
    character(len=*),           intent(in)    :: msg
    real(R8),                   intent(in)    :: extra(:)
    type(morph_status_t),       intent(out)   :: status
    logical,                    intent(out)   :: ok
    real(R8), allocatable :: got(:)
    integer(I4), allocatable :: counts(:)
    integer(I4) :: nr, m, p
    character(len=MORPH_MSGLEN) :: txt

    ok = .false.
    nr = exch%nranks()
    m  = 1 + size(extra)
    allocate( counts(0:nr-1), got(m*nr) )
    counts = m
    call exch%allgather ( [ real(code,R8), extra ], got, counts, status )
    if ( .not. status%ok() ) return
    if ( code /= MORPH_OK ) then
      call status%fail ( code, msg )
      return
    endif
    do p = 0, nr-1
      if ( got(m*p+1) /= 0.0_R8 ) then
        write(txt,'(A,A,I0,A,I0,A)') what, ': refused on rank ', p, ' (status ', nint(got(m*p+1)), ')'
        call status%fail ( nint(got(m*p+1), I4), txt )
        return
      endif
    enddo
    do p = 1, nr-1                ! compared as bits: the same NaN is the same argument
      if ( any( transfer(got(m*p+2:m*p+m), 0_I8, m-1) /= transfer(got(2:m), 0_I8, m-1) ) ) then
        write(txt,'(A,A,I0,A)') what, ': rank ', p, ' was given other arguments than rank 0'
        call status%fail ( MORPH_ERR_INPUT, txt )
        return
      endif
    enddo
    ok = .true.
  end subroutine share_verdicts


  !> This rank's verdict on the grid, phi(nbn, nbn, slot) and the exchange.
  !> A grid not initialised (its Init refused or never called) is refused
  !> before any size of its arrays is taken.
  subroutine check_grid_args ( grid, exch, phi, what, code, msg )
    type(morph_ls_grid_t),   intent(in)    :: grid
    class(morph_exchange_t), intent(inout) :: exch
    real(R8),                intent(in)    :: phi(:,:,:)
    character(len=*),        intent(in)    :: what
    integer(I4),             intent(out)   :: code
    character(len=*),        intent(out)   :: msg
    integer(I4) :: n, ns

    code = MORPH_OK
    msg  = ''
    if ( .not. allocated(grid%blk) ) then
      code = MORPH_ERR_SHAPE
      msg  = what // ': the grid is not initialised'
      return
    endif
    n  = grid%nbn
    ns = size(grid%blk)
    if ( size(phi,1) /= n .or. size(phi,2) /= n .or. size(phi,3) /= ns ) then
      code = MORPH_ERR_SHAPE
      write(msg,'(A,A,I0,A,I0,A,I0,A,I0,A,I0,A)') what, ': phi is ', size(phi,1), ' x ', size(phi,2), ' x ', &
            size(phi,3), ' (the grid: ', n, ' nodes per block side, ', ns, ' blocks here)'
    else if ( exch%rank() /= grid%me .or. exch%nranks() /= grid%nranks ) then
      code = MORPH_ERR_EXCHANGE
      write(msg,'(A,A,I0,A,I0,A,I0,A,I0,A)') what, ': the exchange is rank ', exch%rank(), ' of ', &
            exch%nranks(), ', the grid rank ', grid%me, ' of ', grid%nranks
    endif
  end subroutine check_grid_args


  !> E on the grid of blocks, this rank's blocks phi(nbn, nbn, slot): each
  !> block's part (Morph_LS_Eikonal_Terms, the ghosts from Morph_LS_Halo_Fill,
  !> the sides that face the edge of the grid or an inactive block as edges),
  !> then the sums and the counts of the active blocks in block order
  !> (ordered_sum); e = 0 when no node is in the band. A collective; the same
  !> e on every rank. Refused on every rank, e = 0, when one rank refuses:
  !> a grid not initialised or shapes that are not the grid's
  !> (MORPH_ERR_SHAPE), an exchange that is not the grid's
  !> (MORPH_ERR_EXCHANGE); an exchange that fails returns its status.
  subroutine Morph_LS_Eikonal_Error_Grid ( grid, exch, phi, band, e, status )
    type(morph_ls_grid_t),   intent(in)    :: grid
    class(morph_exchange_t), intent(inout) :: exch
    real(R8),                intent(in)    :: phi(:,:,:)
    real(R8),                intent(in)    :: band
    real(R8),                intent(out)   :: e
    type(morph_status_t),    intent(out)   :: status
    real(R8), allocatable :: w(:,:,:), ps(:), pn(:)
    integer(I4), allocatable :: act(:), pc(:)
    integer(I4) :: n, ns, code, s, b
    real(R8)    :: stot, ntot
    logical     :: ok
    character(len=MORPH_MSGLEN) :: msg

    call status%clear()
    e = 0.0_R8
    call check_grid_args ( grid, exch, phi, 'level-set eikonal error', code, msg )
    call share_verdicts ( exch, 'level-set eikonal error', code, msg, [ real(R8) :: ], status, ok )
    if ( .not. ok ) return

    n  = grid%nbn
    ns = size(grid%blk)
    allocate( w(1-NG:n+NG, 1-NG:n+NG, ns), ps(ns), pc(ns), pn(ns) )
    w = 0.0_R8
    w(1:n,1:n,:) = phi
    call Morph_LS_Halo_Fill ( grid, exch, w, status )
    if ( .not. status%ok() ) return
    !$omp parallel do private(s)
    do s = 1, ns
      call Morph_LS_Eikonal_Terms ( w(:,:,s), grid%h, band, grid%nbr(:,s) == 0, ps(s), pc(s) )
    enddo
    !$omp end parallel do
    pn  = real(pc, R8)
    act = pack( [ (b, b = 1, size(grid%owner)) ], grid%active )
    call exch%ordered_sum ( size(act), grid%owner(act), ps, stot, status )
    if ( .not. status%ok() ) return
    call exch%ordered_sum ( size(act), grid%owner(act), pn, ntot, status )
    if ( .not. status%ok() ) return
    if ( ntot > 0.0_R8 ) e = stot / ntot
  end subroutine Morph_LS_Eikonal_Error_Grid


  !> iters steps of the reinitialisation on this rank's blocks phi(nbn, nbn,
  !> slot), the kernels of Morph_LS_Reinit_Full on each block (THE BLOCKS):
  !> the ghosts of phi0 once, then at every step those of phi and, with
  !> beta > 0, those of the step. A collective. phi is assigned at the end
  !> only: when one rank refuses or an exchange fails, phi is untouched on
  !> every rank. Refused on every rank (its own message on the refusing rank,
  !> the first refusing rank's code elsewhere): a grid not initialised or
  !> shapes that are not the grid's (MORPH_ERR_SHAPE), an exchange that is
  !> not the grid's (MORPH_ERR_EXCHANGE), iters < 0 or beta < 0, or iters
  !> and beta that are not the same on every rank (MORPH_ERR_INPUT: the ranks
  !> would make different numbers of exchanges).
  subroutine Morph_LS_Reinit_Grid ( grid, exch, phi, iters, beta, status )
    type(morph_ls_grid_t),   intent(in)    :: grid
    class(morph_exchange_t), intent(inout) :: exch
    real(R8),                intent(inout) :: phi(:,:,:)
    integer(I4),             intent(in)    :: iters
    real(R8),                intent(in)    :: beta
    type(morph_status_t),    intent(out)   :: status
    real(R8), allocatable :: u(:,:,:), p0w(:,:,:), w(:,:,:), stw(:,:,:), s(:,:,:), dd(:,:,:), st(:,:,:)
    logical,  allocatable :: cut(:,:,:), pos(:,:,:)
    integer(I4) :: n, ns, code, sl, it
    logical     :: ok
    character(len=MORPH_MSGLEN) :: msg

    call status%clear()
    call check_grid_args ( grid, exch, phi, 'level-set reinitialisation', code, msg )
    if ( code == MORPH_OK .and. ( iters < 0 .or. .not. (beta >= 0.0_R8) ) ) then
      code = MORPH_ERR_INPUT
      write(msg,'(A,I0,A,ES10.3,A)') 'level-set reinitialisation: iters = ', iters, ', beta = ', beta, &
            ' (iters >= 0, beta >= 0)'
    endif
    call share_verdicts ( exch, 'level-set reinitialisation', code, msg, [ real(iters,R8), beta ], status, ok )
    if ( .not. ok ) return
    if ( iters == 0 ) return

    n  = grid%nbn
    ns = size(grid%blk)
    allocate( p0w(1-NG:n+NG, 1-NG:n+NG, ns), w(1-NG:n+NG, 1-NG:n+NG, ns), stw(1-NG:n+NG, 1-NG:n+NG, ns) )
    allocate( s(n,n,ns), dd(n,n,ns), st(n,n,ns), cut(n,n,ns), pos(n,n,ns) )
    p0w = 0.0_R8
    w   = 0.0_R8
    stw = 0.0_R8
    u   = phi
    p0w(1:n,1:n,:) = phi
    call Morph_LS_Halo_Fill ( grid, exch, p0w, status )
    if ( .not. status%ok() ) return
    !$omp parallel do private(sl)
    do sl = 1, ns
      call Morph_LS_Reinit_Setup ( p0w(:,:,sl), grid%h, grid%nbr(:,sl) == 0, s(:,:,sl), pos(:,:,sl), &
                                   cut(:,:,sl), dd(:,:,sl) )
    enddo
    !$omp end parallel do

    do it = 1, iters
      w(1:n,1:n,:) = u
      call Morph_LS_Halo_Fill ( grid, exch, w, status )
      if ( .not. status%ok() ) return
      !$omp parallel do private(sl)
      do sl = 1, ns
        call Morph_LS_Reinit_Advance ( w(:,:,sl), pos(:,:,sl), s(:,:,sl), cut(:,:,sl), dd(:,:,sl), grid%h, &
                                       st(:,:,sl) )
      enddo
      !$omp end parallel do
      if ( beta > 0.0_R8 ) then
        stw(1:n,1:n,:) = st
        call Morph_LS_Halo_Fill ( grid, exch, stw, status )
        if ( .not. status%ok() ) return
        !$omp parallel do private(sl)
        do sl = 1, ns
          call Morph_LS_Reinit_Correct ( stw(:,:,sl), p0w(:,:,sl), cut(:,:,sl), grid%h, beta, &
                                         grid%nbr(:,sl) == 0, u(:,:,sl) )
        enddo
        !$omp end parallel do
      else
        u = st
      endif
    enddo
    phi = u
  end subroutine Morph_LS_Reinit_Grid


  !> E on the grid of blocks (Morph_LS_Eikonal_Error_Grid); when E >
  !> threshold, Morph_LS_Reinit_Grid(iters, beta) and done = .true.;
  !> otherwise phi untouched and done = .false. A collective: the decision is
  !> shared, and a threshold that is not the same on every rank makes every
  !> rank refuse (MORPH_ERR_INPUT), phi untouched.
  subroutine Morph_LS_Reinit_If_Needed_Grid ( grid, exch, phi, band, threshold, iters, beta, done, e, status )
    type(morph_ls_grid_t),   intent(in)    :: grid
    class(morph_exchange_t), intent(inout) :: exch
    real(R8),                intent(inout) :: phi(:,:,:)
    real(R8),                intent(in)    :: band, threshold
    integer(I4),             intent(in)    :: iters
    real(R8),                intent(in)    :: beta
    logical,                 intent(out)   :: done
    real(R8),                intent(out)   :: e
    type(morph_status_t),    intent(out)   :: status
    logical :: ok

    call status%clear()
    done = .false.
    call Morph_LS_Eikonal_Error_Grid ( grid, exch, phi, band, e, status )
    if ( .not. status%ok() ) return
    call share_verdicts ( exch, 'level-set reinitialisation if needed', MORPH_OK, '', [ threshold ], status, ok )
    if ( .not. ok ) return
    if ( .not. (e > threshold) ) return
    call Morph_LS_Reinit_Grid ( grid, exch, phi, iters, beta, status )
    done = status%ok()
  end subroutine Morph_LS_Reinit_If_Needed_Grid

end module Morph_LS_Reinit
