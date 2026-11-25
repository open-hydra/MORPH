!>@brief Unit test of the level-set kernels on whole 2-D grids (Morph_LS_HJ,
!> Morph_LS_Reinit). Grid spacing h = 1/32 with power-of-two origins and
!> slopes, so that the checks marked exact are exact by construction:
!>   * WENO-Z: on cubic data, both one-sided derivatives exact up to rounding
!>     (each third-order candidate is exact for a cubic); next to a jump, the
!>     value its weights give in closed form; on data scaled by a power of two,
!>     the value scaled by it, bit for bit (the epsilon scales with the data);
!>   * Godunov: on a linear field with slopes 2 and -1/2, |grad phi| =
!>     sqrt(4.25) at every node, edges included, in both branches, exactly (the
!>     extrapolated ghosts of a linear field are exact); at the bottom of a
!>     valley |x - x_v|, 0 for a front moving along +n and 1 for one moving
!>     along -n (the two branches);
!>   * the step: on a plane front of unit slope with V = -1, phi + dt up to
!>     the last stage's rounding; on a cubic profile, phi0(x + dt) up to
!>     rounding, which only a third-order step gives; a zero velocity
!>     everywhere leaves phi untouched, bit for bit;
!>   * a circle of radius 1 (32 h) centred on a node away from the grid's
!>     centre, solid inside and fluid inside, V = -1 over a stroke of h: the
!>     error at the nodes within 2h of the front is the same for sub-steps of
!>     h/4, h/10 and h/100 (the spatial error dominates) and stays small;
!>   * the eikonal error: exactly 0 on a linear field of unit slope, exactly
!>     1/32 on (x^2 - 1)/2; the reinitialisation runs only strictly above its
!>     threshold;
!>   * the reinitialisation: the signed distance of a plane is a fixed point,
!>     bit for bit; the -0 of an edge node keeps its sign; the circle's
!>     distance is left alone (below the threshold); on a field with the
!>     circle's zero level and a gradient between 0.78 and 1.23, four
!>     iterations bring the eikonal error below the threshold, the zero
!>     crossings along the grid edges hardly move, and two values recorded
!>     when the kernels were written come back within the spread of the builds;
!>   * refused arguments: the status codes, phi untouched; the eikonal error
!>     of a grid with a single node along an axis is 0.
!> Next to each bound that is not exact: its mechanism and the measurement.
!> A NaN fails every check. The values printed must be the same at 1 and 8
!> threads (run_morph_guards.sh units, ctest unit_*_threads_identical).
program test_morph_ls
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_class, ieee_negative_zero, operator(==)
  use Morph_Types_m
  use Morph_LS_HJ
  use Morph_LS_Reinit

  implicit none

  real(R8), parameter :: H   = 1.0_R8 / 32.0_R8
  real(R8), parameter :: EPS = epsilon(1.0_R8)
  real(R8), parameter :: CX  = 0.125_R8, CY = -0.09375_R8, RAD = 1.0_R8   ! the circle
  integer(I4), parameter :: NXC = 117, NYC = 113                          ! x in [-1.75, 1.875], y in [-1.75, 1.75]

  integer(I4) :: nfail
  nfail = 0
  write(*,'(A)') '== Morph level-set kernels =='

  call weno_cubic ()
  call weno_jump ()
  call weno_flat ()
  call weno_scale ()
  call godunov_linear ()
  call godunov_valley ()
  call plane_front ()
  call cubic_step ()
  call zero_velocity ()
  call circle ()
  call eikonal_linear ()
  call eikonal_quadratic ()
  call reinit_plane ()
  call reinit_signed_zero ()
  call reinit_checks ()
  call refusals ()

  if ( nfail > 0 ) then
    write(*,'(I0,A)') nfail, ' LEVEL-SET CHECK(S) FAILED'
    error stop 1
  endif
  write(*,'(A)') 'RESULT: all level-set checks PASS'

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


  !> Node coordinate i (1-based) of an axis with origin x0
  pure real(R8) function xn ( x0, i )
    real(R8),    intent(in) :: x0
    integer(I4), intent(in) :: i
    xn = x0 + real(i-1,R8)*H
  end function xn


  !> The cubic profile of cubic_step: p(s) = s^3/8 + s/2 - 1/16, p' >= 1/2;
  !> exact on these dyadic nodes (every partial result fits in 53 bits)
  pure real(R8) function pcub ( s )
    real(R8), intent(in) :: s
    pcub = (((0.125_R8*s)*s)*s + 0.5_R8*s) - 0.0625_R8
  end function pcub


  ! ---------------------------------------------------------------------------
  !> A cubic: both one-sided derivatives exact up to rounding at the nodes
  !> whose stencils reach no ghost. The differences of these dyadic values are
  !> exact, so the rounding is one per candidate (the division by 6) and about
  !> five in the weighted mean, u = eps/2 of the derivative each: below 4 eps
  !> of |derivative| a priori; bound 8 eps of the largest.
  subroutine weno_cubic ()
    integer(I4), parameter :: NX = 40, NY = 3
    real(R8), parameter :: X0 = -0.625_R8
    real(R8), parameter :: C3 = 0.25_R8, C2 = -0.5_R8, C1 = 0.125_R8, C0 = 0.0625_R8
    real(R8) :: u(1-MORPH_LS_NG:NX+MORPH_LS_NG, 1-MORPH_LS_NG:NY+MORPH_LS_NG), dm(NX,NY), dp(NX,NY)
    real(R8) :: x, ex(NX), worst, tol
    integer(I4) :: i, j
    logical :: ok

    u = 0.0_R8
    do i = 1, NX
      x = xn(X0, i)
      u(i,1:NY) = ((C3*x + C2)*x + C1)*x + C0
      ex(i)     = (3.0_R8*C3*x + 2.0_R8*C2)*x + C1
    enddo
    call Morph_LS_Ghosts_Extrap_All ( u )
    call Morph_LS_Derivs ( u, 1, H, dm, dp )
    tol   = 8.0_R8*EPS*maxval(abs(ex(4:NX-3)))         ! the stencils that reach no ghost
    worst = 0.0_R8
    ok    = .true.
    do j = 1, NY
      do i = 4, NX-3
        worst = max( worst, abs(dm(i,j) - ex(i)), abs(dp(i,j) - ex(i)) )
        ok = ok .and. abs(dm(i,j) - ex(i)) <= tol .and. abs(dp(i,j) - ex(i)) <= tol
      enddo
    enddo
    write(*,'(A,ES10.3,A,F6.2,A)') '   WENO-Z on a cubic: max |derivative - exact| ', worst, ' (', &
          worst/(EPS*maxval(abs(ex(4:NX-3)))), ' eps of the largest derivative)'
    call check ( ok, 'WENO-Z on cubic data: exact up to rounding (8 eps of the largest derivative)' )
  end subroutine weno_cubic


  ! ---------------------------------------------------------------------------
  !> Differences 0, 0, 0, 0, 1: a jump inside the third stencil only. The
  !> smoothness indicators are 0, 0 and 4/3, so tau5 = 4/3, eps = 1e-6 (1e-6
  !> times the largest difference squared), and the weights are
  !> a1 + a2 = 0.7 (1 + (tau5/eps)^2) and a3 = 0.3 (1 + (tau5/(4/3 + eps))^2):
  !> the value is the third candidate, -1/6, times a3/(a1 + a2 + a3), -8.04e-14.
  !> The closed form below evaluates it in another order: about 15 roundings of
  !> u = eps/2 relative each, the squares counting twice, below 10 eps a
  !> priori; bound 32 eps relative.
  subroutine weno_jump ()
    real(R8) :: d, t, e, a12, a3, ex

    d   = Morph_LS_WenoZ ( 0.0_R8, 0.0_R8, 0.0_R8, 0.0_R8, 1.0_R8 )
    t   = 4.0_R8 / 3.0_R8
    e   = 1.0e-6_R8
    a12 = 0.7_R8 * (1.0_R8 + (t/e)**2)
    a3  = 0.3_R8 * (1.0_R8 + (t/(t + e))**2)
    ex  = -(a3/6.0_R8) / (a12 + a3)
    write(*,'(A,ES12.5,A,ES12.5,A,F6.2,A)') '   WENO-Z next to a jump (differences 0, 0, 0, 0, 1): ', d, &
          ', closed form ', ex, ' (', abs(d - ex)/(EPS*abs(ex)), ' eps relative)'
    call check ( abs(d - ex) <= 32.0_R8*EPS*abs(ex), 'WENO-Z next to a jump: the value of its weights (32 eps relative)' )
  end subroutine weno_jump


  ! ---------------------------------------------------------------------------
  !> Flat data, all five differences zero: the derivative is 0, not 0/0 (the
  !> 1e-40 keeps eps positive). Read on the derivatives themselves: Godunov's
  !> clipping max(min(d,0)^2, ...) would hide a NaN here, since gfortran's
  !> min and max return the number when the other argument is a NaN.
  subroutine weno_flat ()
    integer(I4), parameter :: NX = 8, NY = 5
    real(R8) :: u(1-MORPH_LS_NG:NX+MORPH_LS_NG, 1-MORPH_LS_NG:NY+MORPH_LS_NG), dm(NX,NY), dp(NX,NY)

    u = 0.375_R8
    call Morph_LS_Derivs ( u, 1, H, dm, dp )
    call check ( Morph_LS_WenoZ ( 0.0_R8, 0.0_R8, 0.0_R8, 0.0_R8, 0.0_R8 ) == 0.0_R8 .and. &
                 all( dm == 0.0_R8 ) .and. all( dp == 0.0_R8 ), 'WENO-Z on flat data: 0 (not 0/0)' )
  end subroutine weno_flat


  ! ---------------------------------------------------------------------------
  !> Data scaled by a power of two c: WENO-Z(c v) = c WENO-Z(v), bit for bit.
  !> Every operation of the scheme is homogeneous once eps scales with the data
  !> (1e-6 times the largest difference squared; the 1e-40 is far below half an
  !> ulp of it here), and a power-of-two factor commutes with rounding, fused
  !> multiply-adds included. A fixed eps would change the weights with c.
  subroutine weno_scale ()
    real(R8), parameter :: V(5) = [0.3_R8, -0.7_R8, 1.1_R8, 2.9_R8, -0.4_R8]   ! rough: the weights matter
    real(R8), parameter :: C(2) = [2.0_R8**20, 2.0_R8**(-20)]
    real(R8)    :: d, dc
    integer(I4) :: k
    logical     :: ok

    d  = Morph_LS_WenoZ ( V(1), V(2), V(3), V(4), V(5) )
    ok = ieee_is_finite(d)
    do k = 1, 2
      dc = Morph_LS_WenoZ ( C(k)*V(1), C(k)*V(2), C(k)*V(3), C(k)*V(4), C(k)*V(5) )
      ok = ok .and. dc == C(k)*d
    enddo
    call check ( ok, 'WENO-Z on data scaled by 2^20 and by 2^-20: the value scaled, bit for bit' )
  end subroutine weno_scale


  ! ---------------------------------------------------------------------------
  subroutine godunov_linear ()
    integer(I4), parameter :: NX = 20, NY = 17
    real(R8), parameter :: X0 = -0.25_R8, Y0 = -0.375_R8
    real(R8) :: u(1-MORPH_LS_NG:NX+MORPH_LS_NG, 1-MORPH_LS_NG:NY+MORPH_LS_NG), g(NX,NY)
    logical  :: pos(NX,NY)
    integer(I4) :: i, j

    u = 0.0_R8
    do j = 1, NY
      do i = 1, NX
        u(i,j) = (2.0_R8*xn(X0,i) - 0.5_R8*xn(Y0,j)) + 0.125_R8
      enddo
    enddo
    call Morph_LS_Ghosts_Extrap_All ( u )
    pos = .true.
    call Morph_LS_Grad_Godunov ( u, pos, H, g )
    call check ( all( g == sqrt(4.25_R8) ), 'Godunov, front along +n, linear field: |grad phi| exact at every node' )
    pos = .false.
    call Morph_LS_Grad_Godunov ( u, pos, H, g )
    call check ( all( g == sqrt(4.25_R8) ), 'Godunov, front along -n, linear field: |grad phi| exact at every node' )
  end subroutine godunov_linear


  ! ---------------------------------------------------------------------------
  !> phi = |x - x_v| with the bottom of the valley on a node: there phi^- is
  !> -1 and phi^+ is +1, each up to the weight of the WENO-Z stencils that
  !> cross the bottom (2.1e-13, from the weights as for the jump). A front
  !> moving along +n sees |grad phi| = 0 there, exactly (both one-sided terms
  !> are cut to zero); one moving along -n sees 1 - 2.1e-13; bound 1e-12.
  subroutine godunov_valley ()
    integer(I4), parameter :: NX = 21, NY = 7, IV = 11
    real(R8), parameter :: X0 = -0.3125_R8
    real(R8) :: u(1-MORPH_LS_NG:NX+MORPH_LS_NG, 1-MORPH_LS_NG:NY+MORPH_LS_NG), g(NX,NY)
    logical  :: pos(NX,NY)
    integer(I4) :: i

    u = 0.0_R8
    do i = 1, NX
      u(i,1:NY) = abs( xn(X0,i) - xn(X0,IV) )
    enddo
    call Morph_LS_Ghosts_Extrap_All ( u )
    pos = .true.
    call Morph_LS_Grad_Godunov ( u, pos, H, g )
    write(*,'(A,ES11.4)') '   Godunov at the bottom of a valley, front along +n: ', g(IV,4)
    call check ( all( g(IV,:) == 0.0_R8 ), 'Godunov at the bottom of a valley, front along +n: |grad phi| = 0 exactly' )
    pos = .false.
    call Morph_LS_Grad_Godunov ( u, pos, H, g )
    write(*,'(A,ES11.4)') '   Godunov at the bottom of a valley, front along -n: 1 + ', g(IV,4) - 1.0_R8
    call check ( all( abs(g(IV,:) - 1.0_R8) <= 1.0e-12_R8 ), &
                 'Godunov at the bottom of a valley, front along -n: |grad phi| = 1 (1e-12)' )
  end subroutine godunov_valley


  ! ---------------------------------------------------------------------------
  !> A plane front of unit slope, V = -1: the rate is 1 exactly at every node,
  !> so one step gives phi + dt up to the rounding of the last two stages
  !> (u/3 and 2/3 are not exact): about five roundings of u = eps/2 of |phi|,
  !> below 3 eps a priori; bound 4 eps of max |phi|.
  subroutine plane_front ()
    integer(I4), parameter :: NX = 24, NY = 9
    real(R8), parameter :: X0 = -0.5_R8
    real(R8) :: phi(NX,NY), ex(NX,NY), v(NX,NY), dt, tol
    type(morph_status_t) :: st
    integer(I4) :: i

    do i = 1, NX
      phi(i,:) = xn(X0, i) - 0.25_R8
    enddo
    v   = -1.0_R8
    dt  = 0.25_R8*H
    ex  = phi + dt
    tol = 4.0_R8*EPS*maxval(abs(ex))
    call Morph_LS_HJ_Step_Full ( phi, v, H, dt, st )
    write(*,'(A,ES10.3,A,F6.2,A)') '   plane front, one step: max |phi - (phi0 + dt)| ', maxval(abs(phi - ex)), &
          ' (', maxval(abs(phi - ex))/(EPS*maxval(abs(ex))), ' eps of max |phi|)'
    call check ( st%ok() .and. all( abs(phi - ex) <= tol ), &
                 'plane front, V = -1: one step gives phi0 + dt up to rounding (4 eps of max |phi|)' )
  end subroutine plane_front


  ! ---------------------------------------------------------------------------
  !> One step on the cubic profile p (pcub) along one axis, the other flat,
  !> V = -1: phi_t = |grad phi| = p' > 0, so phi(s, t) = p(s + t). On a cubic,
  !> WENO-Z and Godunov return the exact derivative D, so the step applies
  !> 1 + dt D + (dt D)^2/2 + (dt D)^3/6 (the three-stage third-order step on a
  !> linear problem), the Taylor series of p(s + dt), complete for a cubic:
  !> p(s + dt) up to rounding at the nodes the extrapolated ghosts do not reach
  !> within the step (9 from each edge: three stages of three). A second-order
  !> step would miss (dt^3/6) p''' = 5.96e-8. Rounding: about five roundings of
  !> u = eps/2 of |phi| in the last two stages, below 3 eps a priori; bound 8 eps
  !> of max |phi|.
  subroutine cubic_step ()
    integer(I4), parameter :: N = 40, M = 5
    real(R8), parameter :: S0 = -0.625_R8
    real(R8) :: phx(N,M), phy(M,N), vx(N,M), vy(M,N), ex(N), dt, tol, worst
    type(morph_status_t) :: st1, st2
    integer(I4) :: i
    logical :: ok

    dt = 0.25_R8*H
    do i = 1, N
      phx(i,:) = pcub( xn(S0,i) )
      phy(:,i) = pcub( xn(S0,i) )
      ex(i)    = pcub( xn(S0,i) + dt )
    enddo
    vx = -1.0_R8
    vy = -1.0_R8
    call Morph_LS_HJ_Step_Full ( phx, vx, H, dt, st1 )
    call Morph_LS_HJ_Step_Full ( phy, vy, H, dt, st2 )
    tol   = 8.0_R8*EPS*maxval(abs(ex(10:N-9)))
    worst = 0.0_R8
    ok    = st1%ok() .and. st2%ok()
    do i = 10, N-9
      worst = max( worst, maxval(abs(phx(i,:) - ex(i))), maxval(abs(phy(:,i) - ex(i))) )
      ok = ok .and. all( abs(phx(i,:) - ex(i)) <= tol ) .and. all( abs(phy(:,i) - ex(i)) <= tol )
    enddo
    write(*,'(A,ES10.3,A,F6.2,A)') '   cubic profile, one step: max |phi - p(s + dt)| ', worst, ' (', &
          worst/(EPS*maxval(abs(ex(10:N-9)))), ' eps of max |phi|)'
    call check ( ok, 'cubic profile along x and along y, V = -1: one step gives p(s + dt) up to rounding '// &
                 '(8 eps of max |phi|): the step is third order' )
  end subroutine cubic_step


  ! ---------------------------------------------------------------------------
  !> A zero velocity everywhere leaves phi untouched, bit for bit: the step
  !> returns before any arithmetic. The stage arithmetic alone would not
  !> always: built without fused multiply-adds it moves the last bit of 28 of
  !> the circle's 13221 values (with them it gives every value back: the
  !> rounding errors of x/3 and of fl(2/3) x stay within half an ulp of x);
  !> and in any build a velocity of -0 (a zero regression rate, negated)
  !> makes the rate +0, which turns the -0 of a node on the front into +0.
  subroutine zero_velocity ()
    real(R8) :: phi(NXC,NYC), phi0(NXC,NYC), v(NXC,NYC)
    type(morph_status_t) :: st1, st2

    call circle_field ( 1.0_R8, phi )
    phi0 = phi
    v = 0.0_R8
    call Morph_LS_HJ_Step_Full ( phi, v, H, 0.25_R8*H, st1 )
    phi(5,7) = sign( 0.0_R8, -1.0_R8 )
    phi0 = phi
    v = sign( 0.0_R8, -1.0_R8 )
    call Morph_LS_HJ_Step_Full ( phi, v, H, 0.25_R8*H, st2 )
    call check ( st1%ok() .and. st2%ok() .and. all( phi == phi0 ) .and. ieee_class(phi(5,7)) == ieee_negative_zero, &
                 'zero velocity everywhere (+0, then -0): phi untouched, bit for bit, a node''s -0 included' )
  end subroutine zero_velocity


  ! ---------------------------------------------------------------------------
  !> phi = sg (r - RAD) of the circle on the NXC x NYC grid
  subroutine circle_field ( sg, phi )
    real(R8), intent(in)  :: sg
    real(R8), intent(out) :: phi(NXC,NYC)
    integer(I4) :: i, j
    real(R8)    :: dx, dy
    do j = 1, NYC
      do i = 1, NXC
        dx = xn(-1.75_R8, i) - CX
        dy = xn(-1.75_R8, j) - CY
        phi(i,j) = sg * (sqrt(dx*dx + dy*dy) - RAD)
      enddo
    enddo
  end subroutine circle_field


  !> A stroke of h with V = -1 in nsteps sub-steps: the error against the exact
  !> phi0 + t at the nodes within 2h of the exact front, in units of h; fine
  !> is false when a node is not finite or a step failed.
  subroutine circle_run ( sg, nsteps, emax, emean, csum, fine )
    real(R8),    intent(in)  :: sg
    integer(I4), intent(in)  :: nsteps
    real(R8),    intent(out) :: emax, emean, csum
    logical,     intent(out) :: fine
    real(R8) :: phi(NXC,NYC), ex(NXC,NYC), v(NXC,NYC), dt, t
    type(morph_status_t) :: st
    integer(I4) :: s, i, j, n

    call circle_field ( sg, phi )
    v  = -1.0_R8
    dt = H / real(nsteps,R8)
    fine = .true.
    do s = 1, nsteps
      call Morph_LS_HJ_Step_Full ( phi, v, H, dt, st )
      fine = fine .and. st%ok()
    enddo
    fine = fine .and. all( ieee_is_finite(phi) )
    t = real(nsteps,R8)*dt
    call circle_field ( sg, ex )
    ex = ex + t
    emax  = 0.0_R8
    emean = 0.0_R8
    csum  = 0.0_R8
    n = 0
    do j = 1, NYC
      do i = 1, NXC
        csum = csum + phi(i,j)
        if ( abs(ex(i,j)) <= 2.0_R8*H ) then
          emax  = max( emax, abs(phi(i,j) - ex(i,j)) )
          emean = emean + abs(phi(i,j) - ex(i,j))
          n = n + 1
        endif
      enddo
    enddo
    emax  = emax / H
    emean = emean / real(n,R8) / H
  end subroutine circle_run


  !> The bounds: the spread is the time error, which the measurement puts at
  !> 5.5e-6 of the error (bound 1e-3); the error is WENO5's spatial error,
  !> measured 1.75e-8 h (bound 1e-6 h).
  subroutine circle ()
    integer(I4), parameter :: NS(3) = [4, 10, 100]
    real(R8) :: emax(3), emean(3), csum(3), sg, spread
    integer(I4) :: o, k
    logical  :: fine(3)
    character(len=12) :: side

    do o = 1, 2
      if ( o == 1 ) then
        sg = 1.0_R8;  side = 'solid inside'
      else
        sg = -1.0_R8; side = 'fluid inside'
      endif
      do k = 1, 3
        call circle_run ( sg, NS(k), emax(k), emean(k), csum(k), fine(k) )
        write(*,'(A,A,A,I0,A,ES11.4,A,ES11.4,A,ES24.16)') '   circle, ', trim(side), ', sub-step h/', &
              NS(k), ': error max ', emax(k), ' h, mean ', emean(k), ' h; sum of phi ', csum(k)
      enddo
      spread = (maxval(emax) - minval(emax)) / maxval(emax)
      write(*,'(A,ES10.3)') '   relative spread of the max error over the three sub-steps: ', spread
      call check ( all(fine), 'circle, '//trim(side)//': every step taken, every node finite' )
      call check ( spread <= 1.0e-3_R8, 'circle, '//trim(side)//': the error does not depend on the sub-step (1e-3)' )
      call check ( maxval(emax) <= 1.0e-6_R8, 'circle, '//trim(side)//': error within 2h of the front below 1e-6 h' )
    enddo
  end subroutine circle


  ! ---------------------------------------------------------------------------
  subroutine eikonal_linear ()
    integer(I4), parameter :: NX = 30, NY = 11
    real(R8) :: phi(NX,NY)
    integer(I4) :: i

    do i = 1, NX
      phi(i,:) = xn(-0.5_R8, i) - 0.25_R8
    enddo
    call check ( Morph_LS_Eikonal_Error ( phi, H, 2.0_R8*H ) == 0.0_R8, &
                 'eikonal error of a linear field of unit slope: exactly 0' )
    call check ( Morph_LS_Eikonal_Error ( phi, H, -1.0_R8 ) == 0.0_R8, 'eikonal error with no node in the band: 0' )
  end subroutine eikonal_linear


  ! ---------------------------------------------------------------------------
  !> phi = (x^2 - 1)/2, flat in y, on x in [-2, 2]: central differences are
  !> exact on a quadratic, so |grad phi| = |x| at the inner nodes; the band
  !> |phi| <= 2h = 1/16 holds the nodes |x| = 30/32 .. 33/32 of each side,
  !> where ||x| - 1| = 2/32, 1/32, 0, 1/32: E = 1/32 exactly. A band whose
  !> edge is the |phi| of the outermost of those nodes, 124/2048, still holds
  !> it (the band includes its edge). The trigger runs only when E is strictly
  !> above the threshold.
  subroutine eikonal_quadratic ()
    integer(I4), parameter :: NX = 129, NY = 3
    real(R8) :: phi(NX,NY), phi0(NX,NY), e
    type(morph_status_t) :: st
    logical  :: done
    integer(I4) :: i

    do i = 1, NX
      phi(i,:) = 0.5_R8*(xn(-2.0_R8,i)**2 - 1.0_R8)
    enddo
    e = Morph_LS_Eikonal_Error ( phi, H, 2.0_R8*H )
    call check ( e == 1.0_R8/32.0_R8, 'eikonal error of (x^2 - 1)/2 over |phi| <= 2h: exactly 1/32' )
    e = Morph_LS_Eikonal_Error ( phi, H, 124.0_R8/2048.0_R8 )
    call check ( e == 1.0_R8/32.0_R8, 'eikonal error, band edge on a node (|phi| = 124/2048): the node counted, 1/32' )
    phi0 = phi
    call Morph_LS_Reinit_If_Needed_Full ( phi, H, 2.0_R8*H, 1.0_R8/32.0_R8, 4, 0.5_R8, done, e, st )
    call check ( st%ok() .and. .not. done .and. all( phi == phi0 ), &
                 'reinitialisation trigger, threshold equal to the error: not run, phi untouched' )
    call Morph_LS_Reinit_If_Needed_Full ( phi, H, 2.0_R8*H, 0.03_R8, 4, 0.5_R8, done, e, st )
    call check ( st%ok() .and. done, 'reinitialisation trigger, threshold 0.03 below the error: run' )
  end subroutine eikonal_quadratic


  ! ---------------------------------------------------------------------------
  !> The signed distance of a plane, x - x0 and y0 - y with x0 and y0 a
  !> quarter of h off a node, is a fixed point of the reinitialisation, bit
  !> for bit: Godunov's |grad phi| is 1 exactly (a linear field), the subcell
  !> distance is phi itself (its three differences are all 1), and the
  !> correction of the zero level vanishes (pj sa - a sj = 0 when nothing
  !> moved).
  subroutine reinit_plane ()
    integer(I4), parameter :: NX = 24, NY = 19
    real(R8) :: px(NX,NY), py(NX,NY), p0x(NX,NY), p0y(NX,NY)
    type(morph_status_t) :: st1, st2
    integer(I4) :: i, j

    do j = 1, NY
      do i = 1, NX
        px(i,j) = xn(-0.375_R8,i) - (0.25_R8 + 0.25_R8*H)
        py(i,j) = (0.125_R8 - 0.25_R8*H) - xn(-0.3125_R8,j)
      enddo
    enddo
    p0x = px
    p0y = py
    call Morph_LS_Reinit_Full ( px, H, 4, 0.5_R8, st1 )
    call Morph_LS_Reinit_Full ( py, H, 4, 0.5_R8, st2 )
    call check ( st1%ok() .and. st2%ok() .and. all( px == p0x ) .and. all( py == p0y ), &
                 'reinitialisation of the signed distance of a plane, along x and along y: untouched, bit for bit' )
  end subroutine reinit_plane


  ! ---------------------------------------------------------------------------
  !> A -0 on an edge node stays -0. An edge node is never cut, its S is -0,
  !> and its update -S (|grad phi| - 1) is a zero whose sign follows that of
  !> |grad phi| - 1. Here the node sits at the bottom of a valley in y,
  !> phi = |y - y0| + (x - x0)/2 with x0 on the edge: the branch of S >= 0 (a
  !> front moving along +n) cuts the valley's terms to zero and sees
  !> |grad phi| = 1/2, so the update is -0 and the node keeps its -0; the
  !> other branch would see sqrt(1/4 + 1) > 1 and turn it into +0.
  subroutine reinit_signed_zero ()
    integer(I4), parameter :: NX = 12, NY = 13, J0 = 7
    real(R8) :: phi(NX,NY)
    type(morph_status_t) :: st
    integer(I4) :: i, j

    do j = 1, NY
      do i = 1, NX
        phi(i,j) = abs( xn(0.0_R8,j) - xn(0.0_R8,J0) ) + 0.5_R8*xn(0.0_R8,i)
      enddo
    enddo
    phi(1,J0) = sign( 0.0_R8, -1.0_R8 )
    call Morph_LS_Reinit_Full ( phi, H, 4, 0.5_R8, st )
    call check ( st%ok() .and. ieee_class(phi(1,J0)) == ieee_negative_zero, &
                 'reinitialisation: the -0 of an edge node at the bottom of a valley stays -0 (branch of S >= 0)' )
  end subroutine reinit_signed_zero


  ! ---------------------------------------------------------------------------
  !> The recorded values: the eikonal error after the reinitialisation of the
  !> distorted field and the sum of phi, as computed when the kernels were
  !> written (gfortran 13.3, aarch64, release flags). Other builds round
  !> differently: over five gfortran builds (-O0 to -O3, fused multiply-adds
  !> on and off) they moved by 3.1e-16 and 0, relative; a compiler that sums
  !> in another order could move them by 8.9e-14 and 2.1e-12 at worst
  !> ((n - 1) eps/2 sum|x| / |sum x|, n = 800 and 13221 terms). The bound,
  !> 1e-10, is 47 times the larger. Of the 67 deliberate changes of the
  !> kernels tried, 38 move them, the least by 5.8e-8 (the edge gradient
  !> zeroed).
  subroutine reinit_checks ()
    real(R8), parameter :: E1_REC = 2.2633461533705723e-2_R8, SUM_REC = 4.9679098459959087e+03_R8
    real(R8), parameter :: RECORD_TOL = 1.0e-10_R8
    real(R8) :: phi(NXC,NYC), phi0(NXC,NYC), sdf(NXC,NYC), e0, e1, move, sp
    type(morph_status_t) :: st
    logical  :: done
    integer(I4) :: i, j

    ! Below the threshold: nothing happens
    call circle_field ( 1.0_R8, phi )
    phi0 = phi
    call Morph_LS_Reinit_If_Needed_Full ( phi, H, 2.0_R8*H, 0.03_R8, 4, 0.5_R8, done, e0, st )
    write(*,'(A,ES11.4)') '   reinitialisation: eikonal error of the circle''s distance ', e0
    call check ( st%ok() .and. .not. done .and. all( phi == phi0 ), &
                 'reinitialisation below the threshold: not run, phi untouched bit for bit' )

    ! The circle's zero level with a gradient between 0.78 and 1.23
    call circle_field ( 1.0_R8, sdf )
    do j = 1, NYC
      do i = 1, NXC
        phi(i,j) = sdf(i,j) * (1.0_R8 + 0.125_R8*xn(-1.75_R8, i))
      enddo
    enddo
    phi0 = phi
    call Morph_LS_Reinit_If_Needed_Full ( phi, H, 2.0_R8*H, 0.03_R8, 4, 0.5_R8, done, e0, st )
    e1   = Morph_LS_Eikonal_Error ( phi, H, 2.0_R8*H )
    move = crossing_move ( phi0, phi )
    sp   = sum(phi)
    write(*,'(A,ES11.4,A,ES24.16,A,ES11.4,A)') '   reinitialisation of a distorted field: eikonal error ', e0, &
          ' -> ', e1, '; zero crossings moved by ', move/H, ' h at most'
    write(*,'(A,ES24.16,A,ES10.3,A,ES10.3,A)') '   sum of phi after it ', sp, '; against the record: ', &
          abs(e1 - E1_REC)/E1_REC, ' (eikonal error), ', abs(sp - SUM_REC)/SUM_REC, ' (sum), relative'
    call check ( st%ok() .and. done .and. all( ieee_is_finite(phi) ), 'reinitialisation above the threshold: run' )
    call check ( e1 < 0.03_R8, 'reinitialisation: the eikonal error back below the threshold that started it' )
    call check ( move <= 5.0e-3_R8*H, 'reinitialisation: the zero crossings along the grid edges move by 0.005 h at most' )
    call check ( abs(e1 - E1_REC) <= RECORD_TOL*E1_REC .and. abs(sp - SUM_REC) <= RECORD_TOL*SUM_REC, &
                 'reinitialisation: the recorded eikonal error and sum of phi come back (1e-10 relative)' )
  end subroutine reinit_checks


  !> The largest move of a zero crossing along a grid edge (linear
  !> interpolation of the two end values) whose sign changes in both fields
  real(R8) function crossing_move ( p, q ) result ( m )
    real(R8), intent(in) :: p(:,:), q(:,:)
    integer(I4) :: i, j
    real(R8)    :: tp, tq
    m = 0.0_R8
    do j = 1, size(p,2)
      do i = 1, size(p,1)-1
        if ( p(i,j)*p(i+1,j) < 0.0_R8 .and. q(i,j)*q(i+1,j) < 0.0_R8 ) then
          tp = p(i,j) / (p(i,j) - p(i+1,j))
          tq = q(i,j) / (q(i,j) - q(i+1,j))
          m = max( m, abs(tp - tq)*H )
        endif
      enddo
    enddo
    do j = 1, size(p,2)-1
      do i = 1, size(p,1)
        if ( p(i,j)*p(i,j+1) < 0.0_R8 .and. q(i,j)*q(i,j+1) < 0.0_R8 ) then
          tp = p(i,j) / (p(i,j) - p(i,j+1))
          tq = q(i,j) / (q(i,j) - q(i,j+1))
          m = max( m, abs(tp - tq)*H )
        endif
      enddo
    enddo
  end function crossing_move


  ! ---------------------------------------------------------------------------
  subroutine refusals ()
    real(R8) :: phi(5,4), phi0(5,4), v(4,5), v2(5,4), one(1,4)
    type(morph_status_t) :: st

    call some_field ( phi )
    phi0 = phi
    v = -1.0_R8
    call Morph_LS_HJ_Step_Full ( phi, v, H, H, st )
    call check ( st%code == MORPH_ERR_SHAPE .and. all(phi == phi0), 'step: V of another shape refused, phi untouched' )
    v2 = -1.0_R8
    call Morph_LS_HJ_Step_Full ( phi, v2, 0.0_R8, H, st )
    call check ( st%code == MORPH_ERR_INPUT .and. all(phi == phi0), 'step: h = 0 refused, phi untouched' )
    call Morph_LS_HJ_Step_Full ( phi, v2, H, -H, st )
    call check ( st%code == MORPH_ERR_INPUT .and. all(phi == phi0), 'step: dt < 0 refused, phi untouched' )
    one = 0.5_R8
    call Morph_LS_HJ_Step_Full ( one, one, H, H, st )
    call check ( st%code == MORPH_ERR_SHAPE, 'step: a single node along an axis refused' )
    call Morph_LS_Reinit_Full ( one, H, 4, 0.5_R8, st )
    call check ( st%code == MORPH_ERR_SHAPE, 'reinitialisation: fewer than 3 nodes along an axis refused' )
    call Morph_LS_Reinit_Full ( phi, H, -1, 0.5_R8, st )
    call check ( st%code == MORPH_ERR_INPUT .and. all(phi == phi0), 'reinitialisation: iters < 0 refused, phi untouched' )
    call Morph_LS_Reinit_Full ( phi, H, 0, 0.5_R8, st )
    call check ( st%ok() .and. all(phi == phi0), 'reinitialisation: zero iterations leave phi untouched' )
    call Morph_LS_Reinit_Full ( phi, 0.0_R8, 4, 0.5_R8, st )
    call check ( st%code == MORPH_ERR_INPUT .and. all(phi == phi0), 'reinitialisation: h = 0 refused, phi untouched' )
    call Morph_LS_Reinit_Full ( phi, H, 4, -0.5_R8, st )
    call check ( st%code == MORPH_ERR_INPUT .and. all(phi == phi0), 'reinitialisation: beta < 0 refused, phi untouched' )
    call check ( Morph_LS_Eikonal_Error ( one, H, 1.0_R8 ) == 0.0_R8, &
                 'eikonal error on a single node along an axis: 0 (no gradient to take)' )
  end subroutine refusals


  !> Some values of both signs (any would do: the refused calls must leave
  !> them as they are)
  subroutine some_field ( a )
    real(R8), intent(out) :: a(:,:)
    integer(I4) :: i, j
    do j = 1, size(a,2)
      do i = 1, size(a,1)
        a(i,j) = real(mod(7*i + 13*j, 11) - 5, R8) * 0.125_R8
      enddo
    enddo
  end subroutine some_field

end program test_morph_ls
