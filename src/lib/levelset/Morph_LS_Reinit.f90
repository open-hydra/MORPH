!>@brief The safety reinitialisation of a level-set field toward a signed
!> distance, and the measure that decides when it is needed, on a whole 2-D
!> grid of nodes x0 + (i,j)*h.
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
!> Never stops, no I/O: inconsistent arguments return MORPH_ERR_SHAPE or
!> MORPH_ERR_INPUT in the status, phi untouched.
module Morph_LS_Reinit
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m
  use Morph_LS_HJ,     only: MORPH_LS_NG, Morph_LS_Ghosts_Extrap_All, Morph_LS_Grad_Godunov

  implicit none
  private

  public :: Morph_LS_Eikonal_Error, Morph_LS_Reinit_Full, Morph_LS_Reinit_If_Needed_Full

contains

  !> The gradient along one axis at every node: (f(i+1) - f(i-1))/(2h)
  !> inside, (f(2) - f(1))/h and (f(n) - f(n-1))/h on the edges.
  subroutine gradient_axis ( f, h, axis, g )
    real(R8),    intent(in)  :: f(:,:)
    real(R8),    intent(in)  :: h
    integer(I4), intent(in)  :: axis
    real(R8),    intent(out) :: g(:,:)
    integer(I4) :: nx, ny, i, j

    nx = size(f,1)
    ny = size(f,2)
    if ( axis == 1 ) then
      !$omp parallel do private(i,j)
      do j = 1, ny
        g(1,j)  = (f(2,j) - f(1,j)) / h
        do i = 2, nx-1
          g(i,j) = (f(i+1,j) - f(i-1,j)) / (2.0_R8*h)
        enddo
        g(nx,j) = (f(nx,j) - f(nx-1,j)) / h
      enddo
      !$omp end parallel do
    else
      !$omp parallel do private(i,j)
      do j = 1, ny
        if ( j == 1 ) then
          g(:,j) = (f(:,2) - f(:,1)) / h
        else if ( j == ny ) then
          g(:,j) = (f(:,ny) - f(:,ny-1)) / h
        else
          g(:,j) = (f(:,j+1) - f(:,j-1)) / (2.0_R8*h)
        endif
      enddo
      !$omp end parallel do
    endif
  end subroutine gradient_axis


  !> E = mean | |grad phi| - 1 | over the nodes with |phi| <= band; 0 when no
  !> node is in the band. Summed in node order (i fastest).
  real(R8) function Morph_LS_Eikonal_Error ( phi, h, band ) result ( e )
    real(R8), intent(in) :: phi(:,:)
    real(R8), intent(in) :: h, band
    real(R8), allocatable :: gx(:,:), gy(:,:)
    real(R8)    :: s
    integer(I4) :: i, j, n

    e = 0.0_R8
    if ( size(phi,1) < 2 .or. size(phi,2) < 2 ) return
    allocate( gx(size(phi,1),size(phi,2)), gy(size(phi,1),size(phi,2)) )
    call gradient_axis ( phi, h, 1, gx )
    call gradient_axis ( phi, h, 2, gy )
    s = 0.0_R8
    n = 0
    do j = 1, size(phi,2)
      do i = 1, size(phi,1)
        if ( abs(phi(i,j)) <= band ) then
          s = s + abs( sqrt(gx(i,j)*gx(i,j) + gy(i,j)*gy(i,j)) - 1.0_R8 )
          n = n + 1
        endif
      enddo
    enddo
    if ( n > 0 ) e = s / real(n,R8)
  end function Morph_LS_Eikonal_Error


  !> iters steps of the reinitialisation of the header, on phi(nx,ny) in place.
  !> Refused, phi untouched: fewer than 3 nodes along an axis
  !> (MORPH_ERR_SHAPE); h <= 0, iters < 0 or beta < 0 (MORPH_ERR_INPUT).
  subroutine Morph_LS_Reinit_Full ( phi, h, iters, beta, status )
    real(R8),             intent(inout) :: phi(:,:)
    real(R8),             intent(in)    :: h
    integer(I4),          intent(in)    :: iters
    real(R8),             intent(in)    :: beta
    type(morph_status_t), intent(out)   :: status
    real(R8), allocatable :: p0(:,:), gx(:,:), gy(:,:), s(:,:), dd(:,:), w(:,:), g(:,:), st(:,:), corr(:,:)
    logical,  allocatable :: cut(:,:), pos(:,:)
    integer(I4) :: nx, ny, i, j, it
    real(R8)    :: gdx, gdy, f, b, a, fs, nc, pj, sj, gn, r
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

    allocate( p0(nx,ny), gx(nx,ny), gy(nx,ny), s(nx,ny), dd(nx,ny), g(nx,ny), st(nx,ny), corr(nx,ny) )
    allocate( cut(nx,ny), pos(nx,ny) )
    allocate( w(1-MORPH_LS_NG:nx+MORPH_LS_NG, 1-MORPH_LS_NG:ny+MORPH_LS_NG) )
    w  = 0.0_R8
    p0 = phi

    ! The sign of phi0, smoothed by the grid: S
    call gradient_axis ( p0, h, 1, gx )
    call gradient_axis ( p0, h, 2, gy )
    !$omp parallel do private(j)
    do j = 1, ny
      s(:,j) = p0(:,j) / sqrt( (p0(:,j)*p0(:,j) + ((gx(:,j)*gx(:,j) + gy(:,j)*gy(:,j))*h)*h) + 1.0e-30_R8 )
    enddo
    !$omp end parallel do
    pos = s >= 0.0_R8

    ! The cut nodes and the distance D of the subcell fix (inside only: the
    ! edge nodes are never cut)
    cut = .false.
    dd  = 0.0_R8
    !$omp parallel do private(i,j,f,b,a,gdx,gdy)
    do j = 2, ny-1
      do i = 2, nx-1
        a = p0(i,j)
        cut(i,j) = a*p0(i-1,j) <= 0.0_R8 .or. a*p0(i+1,j) <= 0.0_R8 .or. &
                   a*p0(i,j-1) <= 0.0_R8 .or. a*p0(i,j+1) <= 0.0_R8
        f = p0(i+1,j)
        b = p0(i-1,j)
        gdx = max( max( abs(f - b)/(2.0_R8*h), abs(f - a)/h ), abs(a - b)/h )
        f = p0(i,j+1)
        b = p0(i,j-1)
        gdy = max( max( abs(f - b)/(2.0_R8*h), abs(f - a)/h ), abs(a - b)/h )
        dd(i,j) = a / sqrt( (gdx*gdx + gdy*gdy) + 1.0e-30_R8 )
      enddo
    enddo
    !$omp end parallel do

    do it = 1, iters
      ! The plain update everywhere, the subcell fix at the cut nodes
      w(1:nx,1:ny) = phi
      call Morph_LS_Ghosts_Extrap_All ( w )
      call Morph_LS_Grad_Godunov ( w, pos, h, g )
      !$omp parallel do private(i,j,r)
      do j = 1, ny
        do i = 1, nx
          if ( cut(i,j) ) then
            r = -(phi(i,j) - dd(i,j)) / h
          else
            r = -s(i,j) * (g(i,j) - 1.0_R8)
          endif
          st(i,j) = phi(i,j) + (0.45_R8*h)*r
        enddo
      enddo
      !$omp end parallel do

      ! The correction of the zero level at the cut nodes, from the values of
      ! st before any correction
      if ( beta > 0.0_R8 ) then
        call gradient_axis ( st, h, 1, gx )
        call gradient_axis ( st, h, 2, gy )
        !$omp parallel do private(i,j,a,fs,nc,pj,sj,gn)
        do j = 1, ny
          do i = 1, nx
            corr(i,j) = st(i,j)
            if ( .not. cut(i,j) ) cycle
            a  = p0(i,j)
            fs = 0.0_R8
            nc = 0.0_R8
            pj = p0(i-1,j); sj = st(i-1,j)
            call add_crossing ( a, pj, st(i,j), sj, fs, nc )
            pj = p0(i+1,j); sj = st(i+1,j)
            call add_crossing ( a, pj, st(i,j), sj, fs, nc )
            pj = p0(i,j-1); sj = st(i,j-1)
            call add_crossing ( a, pj, st(i,j), sj, fs, nc )
            pj = p0(i,j+1); sj = st(i,j+1)
            call add_crossing ( a, pj, st(i,j), sj, fs, nc )
            if ( nc > 0.0_R8 ) then
              gn = sqrt( (gx(i,j)**2 + gy(i,j)**2) + 1.0e-30_R8 )
              corr(i,j) = st(i,j) - ((beta*gn)*fs) / max(nc, 1.0_R8)
            endif
          enddo
        enddo
        !$omp end parallel do
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

end module Morph_LS_Reinit
