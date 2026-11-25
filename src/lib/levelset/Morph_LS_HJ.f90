!>@brief The level-set equation  phi_t + V |grad phi| = 0  on a 2-D grid of
!> nodes x0 + (i,j)*h: fifth-order WENO-Z one-sided derivatives, Godunov's
!> upwind Hamiltonian, and the three-stage strong-stability-preserving
!> Runge-Kutta step (Shu and Osher). phi < 0 is the solid, phi > 0 the fluid,
!> and a regressing wall has V = -(regression rate) <= 0.
!>
!> THE ARRAYS. The kernels work on one array of interior nodes 1..nx, 1..ny
!> with MORPH_LS_NG = 3 ghost layers on each side of each axis (the WENO5
!> stencil): u(1-NG:nx+NG, 1-NG:ny+NG). They do not care who fills the ghosts.
!> Over a whole grid the ghosts are the linear extrapolation of the last two
!> nodes of each line (Morph_LS_Ghosts_Extrap), and Morph_LS_HJ_Step_Full
!> makes one step. A grid split into blocks fills them from the neighbouring
!> blocks between the stages, so the three stages are exposed one by one:
!> fill the ghosts, Morph_LS_Rate, Morph_LS_RK3_Combine.
!>
!> THE SCHEME, at an interior node and along one axis, with the backward
!> differences D(m) = (u(m) - u(m-1))/h:
!>   phi_x^-  = WENO-Z( D(i-2), D(i-1), D(i),   D(i+1), D(i+2) )
!>   phi_x^+  = WENO-Z( D(i+3), D(i+2), D(i+1), D(i),   D(i-1) )
!> WENO-Z (Borges et al. 2008) combines the three third-order candidates with
!> the weights c_k (1 + (tau5/(beta_k + eps))^2), c = (.1, .6, .3),
!> tau5 = |beta_1 - beta_3| and eps = 1e-6 max(v_k^2) + 1e-40: the epsilon
!> scales with the data, so the weights do not depend on the units of phi.
!> Godunov's |grad phi| takes, per axis, max(max(phi^-,0)^2, min(phi^+,0)^2)
!> where the front moves along +n (pos = .true.) and
!> max(min(phi^-,0)^2, max(phi^+,0)^2) where it moves along -n; the caller
!> says which, per node. The step is
!>   u1 = u + dt L(u),  u2 = 3/4 u + 1/4 (u1 + dt L(u1)),
!>   u_new = u/3 + 2/3 (u2 + dt L(u2)),         L(w) = -V |grad w|.
!>
!> THE ORDER OF OPERATIONS is part of the contract: every sum and product is
!> written in the order the parentheses show, so that a build without fused
!> multiply-adds (-ffp-contract=off) gives the same bits as any other
!> implementation that evaluates the same formulas in the same order. A build
!> that fuses them differs in the last bits.
!>
!> Every node is computed from the previous stage alone: no reduction, so the
!> result is the same at any number of OpenMP threads. A zero velocity
!> everywhere leaves phi untouched, bit for bit: Morph_LS_HJ_Step_Full returns
!> before computing anything (the stage arithmetic alone would move the last
!> bits of phi).
!>
!> Never stops, no I/O: inconsistent arguments return MORPH_ERR_SHAPE or
!> MORPH_ERR_INPUT in the status.
module Morph_LS_HJ
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m

  implicit none
  private

  !> Ghost layers on each side of each axis: the WENO5 stencil reaches 3 nodes.
  integer(I4), parameter, public :: MORPH_LS_NG = 3

  real(R8), parameter :: C1312 = 13.0_R8 / 12.0_R8
  real(R8), parameter :: C23   = 2.0_R8 / 3.0_R8

  public :: Morph_LS_WenoZ
  public :: Morph_LS_Ghosts_Extrap, Morph_LS_Ghosts_Extrap_All
  public :: Morph_LS_Derivs, Morph_LS_Godunov, Morph_LS_Grad_Godunov
  public :: Morph_LS_Rate, Morph_LS_RK3_Combine
  public :: Morph_LS_HJ_Step_Full

contains

  !> WENO-Z value of a derivative from five differences, v3 the central one;
  !> see the module header for the weights.
  elemental real(R8) function Morph_LS_WenoZ ( v1, v2, v3, v4, v5 ) result ( d )
    real(R8), intent(in) :: v1, v2, v3, v4, v5
    real(R8) :: p1, p2, p3, s1, s2, b1, b2, b3, e, t, a1, a2, a3

    p1 = ((2.0_R8*v1 - 7.0_R8*v2) + 11.0_R8*v3) / 6.0_R8
    p2 = ((-v2 + 5.0_R8*v3) + 2.0_R8*v4) / 6.0_R8
    p3 = ((2.0_R8*v3 + 5.0_R8*v4) - v5) / 6.0_R8

    s1 = (v1 - 2.0_R8*v2) + v3
    s2 = (v1 - 4.0_R8*v2) + 3.0_R8*v3
    b1 = (C1312*s1)*s1 + (0.25_R8*s2)*s2
    s1 = (v2 - 2.0_R8*v3) + v4
    s2 = v2 - v4
    b2 = (C1312*s1)*s1 + (0.25_R8*s2)*s2
    s1 = (v3 - 2.0_R8*v4) + v5
    s2 = (3.0_R8*v3 - 4.0_R8*v4) + v5
    b3 = (C1312*s1)*s1 + (0.25_R8*s2)*s2

    e = 1.0e-6_R8 * max( max( max( max( v1*v1, v2*v2 ), v3*v3 ), v4*v4 ), v5*v5 ) + 1.0e-40_R8
    t = abs(b1 - b3)
    a1 = 0.1_R8 * (1.0_R8 + (t/(b1 + e))**2)
    a2 = 0.6_R8 * (1.0_R8 + (t/(b2 + e))**2)
    a3 = 0.3_R8 * (1.0_R8 + (t/(b3 + e))**2)
    d = ((a1*p1 + a2*p2) + a3*p3) / ((a1 + a2) + a3)
  end function Morph_LS_WenoZ


  !> The ghosts of one side of one axis (side 1: low, 2: high) by linear
  !> extrapolation of the last two interior nodes of each line: on the low
  !> side u(1-k) = u(1) - k (u(2) - u(1)), k = 1..NG. Lines across the
  !> interior only; the corner ghosts are left as they are (no stencil reads
  !> them). Needs at least two interior nodes along the axis.
  subroutine Morph_LS_Ghosts_Extrap ( u, axis, side )
    real(R8),    intent(inout) :: u(1-MORPH_LS_NG:, 1-MORPH_LS_NG:)
    integer(I4), intent(in)    :: axis, side
    integer(I4) :: nx, ny, i, j, k
    real(R8)    :: d

    nx = size(u,1) - 2*MORPH_LS_NG
    ny = size(u,2) - 2*MORPH_LS_NG
    if ( axis == 1 ) then
      !$omp parallel do private(j,k,d)
      do j = 1, ny
        if ( side == 1 ) then
          d = u(2,j) - u(1,j)
          do k = 1, MORPH_LS_NG
            u(1-k,j) = u(1,j) - real(k,R8)*d
          end do
        else
          d = u(nx,j) - u(nx-1,j)
          do k = 1, MORPH_LS_NG
            u(nx+k,j) = u(nx,j) + real(k,R8)*d
          end do
        endif
      enddo
      !$omp end parallel do
    else
      !$omp parallel do private(i,k,d)
      do i = 1, nx
        if ( side == 1 ) then
          d = u(i,2) - u(i,1)
          do k = 1, MORPH_LS_NG
            u(i,1-k) = u(i,1) - real(k,R8)*d
          end do
        else
          d = u(i,ny) - u(i,ny-1)
          do k = 1, MORPH_LS_NG
            u(i,ny+k) = u(i,ny) + real(k,R8)*d
          end do
        endif
      enddo
      !$omp end parallel do
    endif
  end subroutine Morph_LS_Ghosts_Extrap


  !> All four sides by linear extrapolation: the ghosts of a whole grid.
  subroutine Morph_LS_Ghosts_Extrap_All ( u )
    real(R8), intent(inout) :: u(1-MORPH_LS_NG:, 1-MORPH_LS_NG:)
    integer(I4) :: axis, side
    do axis = 1, 2
      do side = 1, 2
        call Morph_LS_Ghosts_Extrap ( u, axis, side )
      enddo
    enddo
  end subroutine Morph_LS_Ghosts_Extrap_All


  !> The one-sided WENO-Z derivatives phi^- (dm) and phi^+ (dp) along one axis
  !> at every interior node; the ghosts must be filled.
  subroutine Morph_LS_Derivs ( u, axis, h, dm, dp )
    real(R8),    intent(in)  :: u(1-MORPH_LS_NG:, 1-MORPH_LS_NG:)
    integer(I4), intent(in)  :: axis
    real(R8),    intent(in)  :: h
    real(R8),    intent(out) :: dm(:,:), dp(:,:)
    integer(I4) :: nx, ny, i, j

    nx = size(u,1) - 2*MORPH_LS_NG
    ny = size(u,2) - 2*MORPH_LS_NG
    if ( axis == 1 ) then
      !$omp parallel do private(j)
      do j = 1, ny
        call line_derivs ( nx, u(:,j), h, dm(:,j), dp(:,j) )
      enddo
      !$omp end parallel do
    else
      !$omp parallel do private(i)
      do i = 1, nx
        call line_derivs ( ny, u(i,:), h, dm(i,:), dp(i,:) )
      enddo
      !$omp end parallel do
    endif
  end subroutine Morph_LS_Derivs


  !> One line: ul(1-NG:n+NG) -> the one-sided derivatives at its n interior nodes.
  pure subroutine line_derivs ( n, ul, h, dml, dpl )
    integer(I4), intent(in)  :: n
    real(R8),    intent(in)  :: ul(1-MORPH_LS_NG:n+MORPH_LS_NG)
    real(R8),    intent(in)  :: h
    real(R8),    intent(out) :: dml(n), dpl(n)
    real(R8)    :: d(-1:n+3)
    integer(I4) :: i, m

    do m = -1, n+3
      d(m) = (ul(m) - ul(m-1)) / h
    enddo
    do i = 1, n
      dml(i) = Morph_LS_WenoZ ( d(i-2), d(i-1), d(i),   d(i+1), d(i+2) )
      dpl(i) = Morph_LS_WenoZ ( d(i+3), d(i+2), d(i+1), d(i),   d(i-1) )
    enddo
  end subroutine line_derivs


  !> Godunov's |grad phi| from the one-sided derivatives of the two axes;
  !> pos = .true. where the front moves along +grad phi (see the header).
  elemental real(R8) function Morph_LS_Godunov ( dmx, dpx, dmy, dpy, pos ) result ( g )
    real(R8), intent(in) :: dmx, dpx, dmy, dpy
    logical,  intent(in) :: pos
    g = sqrt( axis_term( dmx, dpx, pos ) + axis_term( dmy, dpy, pos ) )
  end function Morph_LS_Godunov


  elemental real(R8) function axis_term ( dm, dp, pos ) result ( s )
    real(R8), intent(in) :: dm, dp
    logical,  intent(in) :: pos
    if ( pos ) then
      s = max( max(dm, 0.0_R8)**2, min(dp, 0.0_R8)**2 )
    else
      s = max( min(dm, 0.0_R8)**2, max(dp, 0.0_R8)**2 )
    endif
  end function axis_term


  !> Godunov's |grad u| at the interior nodes of u, whose ghosts must be filled.
  subroutine Morph_LS_Grad_Godunov ( u, pos, h, g )
    real(R8), intent(in)  :: u(1-MORPH_LS_NG:, 1-MORPH_LS_NG:)
    logical,  intent(in)  :: pos(:,:)
    real(R8), intent(in)  :: h
    real(R8), intent(out) :: g(:,:)
    real(R8), allocatable :: dmx(:,:), dpx(:,:), dmy(:,:), dpy(:,:)
    integer(I4) :: nx, ny, j

    nx = size(u,1) - 2*MORPH_LS_NG
    ny = size(u,2) - 2*MORPH_LS_NG
    allocate( dmx(nx,ny), dpx(nx,ny), dmy(nx,ny), dpy(nx,ny) )
    call Morph_LS_Derivs ( u, 1, h, dmx, dpx )
    call Morph_LS_Derivs ( u, 2, h, dmy, dpy )
    !$omp parallel do private(j)
    do j = 1, ny
      g(:,j) = Morph_LS_Godunov ( dmx(:,j), dpx(:,j), dmy(:,j), dpy(:,j), pos(:,j) )
    enddo
    !$omp end parallel do
  end subroutine Morph_LS_Grad_Godunov


  !> L(u) = -V |grad u| at the interior nodes (ghosts filled); the Godunov
  !> branch is V > 0 at each node.
  subroutine Morph_LS_Rate ( u, vls, h, rate )
    real(R8), intent(in)  :: u(1-MORPH_LS_NG:, 1-MORPH_LS_NG:)
    real(R8), intent(in)  :: vls(:,:)
    real(R8), intent(in)  :: h
    real(R8), intent(out) :: rate(:,:)
    integer(I4) :: j

    call Morph_LS_Grad_Godunov ( u, vls > 0.0_R8, h, rate )
    !$omp parallel do private(j)
    do j = 1, size(rate,2)
      rate(:,j) = -vls(:,j) * rate(:,j)
    enddo
    !$omp end parallel do
  end subroutine Morph_LS_Rate


  !> One stage of the step, in place on us (the stage's state on entry, the
  !> next stage's on exit); u0 is the state at the start of the step:
  !>   stage 1: us + dt L,   stage 2: 3/4 u0 + 1/4 (us + dt L),
  !>   stage 3: u0/3 + 2/3 (us + dt L).
  subroutine Morph_LS_RK3_Combine ( stage, u0, rate, dt, us )
    integer(I4), intent(in)    :: stage
    real(R8),    intent(in)    :: u0(:,:), rate(:,:)
    real(R8),    intent(in)    :: dt
    real(R8),    intent(inout) :: us(:,:)
    integer(I4) :: j

    select case ( stage )
    case ( 1 )
      !$omp parallel do private(j)
      do j = 1, size(us,2)
        us(:,j) = us(:,j) + dt*rate(:,j)
      enddo
      !$omp end parallel do
    case ( 2 )
      !$omp parallel do private(j)
      do j = 1, size(us,2)
        us(:,j) = 0.75_R8*u0(:,j) + 0.25_R8*(us(:,j) + dt*rate(:,j))
      enddo
      !$omp end parallel do
    case default
      !$omp parallel do private(j)
      do j = 1, size(us,2)
        us(:,j) = u0(:,j)/3.0_R8 + C23*(us(:,j) + dt*rate(:,j))
      enddo
      !$omp end parallel do
    end select
  end subroutine Morph_LS_RK3_Combine


  !> One step of length dt over a whole grid phi(nx,ny) (nodes x0 + (i,j)*h),
  !> the ghosts extrapolated at every stage; vls(nx,ny) is V at the nodes.
  !> Refused, phi untouched: shapes that differ or fewer than 2 nodes along an
  !> axis (MORPH_ERR_SHAPE); h <= 0 or dt < 0 (MORPH_ERR_INPUT).
  subroutine Morph_LS_HJ_Step_Full ( phi, vls, h, dt, status )
    real(R8),             intent(inout) :: phi(:,:)
    real(R8),             intent(in)    :: vls(:,:)
    real(R8),             intent(in)    :: h, dt
    type(morph_status_t), intent(out)   :: status
    real(R8), allocatable :: w(:,:), u0(:,:), us(:,:), rate(:,:)
    integer(I4) :: nx, ny, stage
    character(len=MORPH_MSGLEN) :: msg

    call status%clear()
    nx = size(phi,1)
    ny = size(phi,2)
    if ( size(vls,1) /= nx .or. size(vls,2) /= ny .or. nx < 2 .or. ny < 2 ) then
      write(msg,'(A,I0,A,I0,A,I0,A,I0,A)') 'level set step: phi is ', nx, ' x ', ny, ', V ', &
            size(vls,1), ' x ', size(vls,2), ' (same shape, at least 2 x 2)'
      call status%fail ( MORPH_ERR_SHAPE, msg )
      return
    endif
    if ( .not. (h > 0.0_R8) .or. .not. (dt >= 0.0_R8) ) then
      write(msg,'(A,ES10.3,A,ES10.3,A)') 'level set step: h = ', h, ', dt = ', dt, ' (h > 0, dt >= 0)'
      call status%fail ( MORPH_ERR_INPUT, msg )
      return
    endif
    if ( all(vls == 0.0_R8) ) return

    allocate( w(1-MORPH_LS_NG:nx+MORPH_LS_NG, 1-MORPH_LS_NG:ny+MORPH_LS_NG) )
    allocate( rate(nx,ny) )
    w  = 0.0_R8
    u0 = phi
    us = phi
    do stage = 1, 3
      w(1:nx,1:ny) = us
      call Morph_LS_Ghosts_Extrap_All ( w )
      call Morph_LS_Rate ( w, vls, h, rate )
      call Morph_LS_RK3_Combine ( stage, u0, rate, dt, us )
    enddo
    phi = us
  end subroutine Morph_LS_HJ_Step_Full

end module Morph_LS_HJ
