!>@brief The smoothed interior of a receding block: the nodes between the
!> receding face and the pinned face follow an elliptic smoothing of the
!> boundary displacement, so that a deformation of the boundary spreads
!> gradually into the block instead of being confined to each grid line.
!>
!> The displacement field is the equilibrium of lineal springs on the grid
!> edges (Batina 1990), stiffness 1/|e| measured on node_old: at every solved
!> node p,  sum over its edges e = (p,q) of (delta_q - delta_p)/|e| = 0.  That
!> is the discrete Laplacian with inverse-distance weights. Along one straight
!> grid line it makes delta linear in the arc length, which is exactly the line
!> law (node m moves by disp * s_m/L, plan 08 section 2.3): on a Cartesian block
!> with a uniform face displacement the two agree to solver tolerance, and a
!> graded line keeps its grading. Index-space weights, or 1/|e|^2, would not:
!> they compress the small cells of a graded line by a larger fraction and lose
!> the property that no cell can invert before the whole line is consumed.
!>
!> THE BOUNDARY HIERARCHY, lowest dimension first, so that each entity is fixed
!> before the next is solved and a block face two blocks share gets the same
!> data in both (the nodeshift law builds the first three; this module solves
!> the last two):
!>   the receding face     delta = disp(a,b)       prescribed
!>   the pinned face       delta = 0
!>   the four line edges   delta = disp * s/L      the line law's closed form, the
!>                         same evaluation in every block that holds the edge
!>   the four side faces   the in-face springs only, Dirichlet on the four
!>                         edges: the nodes slide within their face
!>   the interior          the six springs, Dirichlet on the six faces
!>
!> THE SOLVER. Each solved set is iterated with two-colour (red-black)
!> successive over-relaxation, from the field the caller hands in (the line
!> blend on a first call, the previous solution later), until the largest update
!> of a sweep is at or below tol times the largest prescribed displacement. A
!> node of one colour reads only nodes of the other and the stopping test is a
!> max-norm, so the result does not depend on the number of threads. It does
!> depend on the sweep order: two 3-D blocks that share a side face WITH
!> interior face nodes agree only to tolerance, not bitwise. A 2-D mesh (two
!> node planes in k) has no such face: its seam faces hold edge nodes only,
!> and the edges are the closed form; its two planes are solved as two faces
!> with the same data, the same weights and the same colouring, and come out
!> bit for bit the same.
!>
!> Never stops, no I/O: a set that does not converge within maxit sweeps, or a
!> coincident pair of nodes (a spring of zero length), returns MORPH_ERR_SMOOTH.
module Morph_Smooth
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m

  implicit none
  private

  public :: morph_smooth_ctl_t, morph_smooth_info_t, Morph_Smooth_Lineal
  public :: MORPH_INTERIOR_LINE, MORPH_INTERIOR_SMOOTH

  integer(I4), parameter :: MORPH_INTERIOR_LINE   = 1   !> the line law: node m moves by disp * s_m/L
  integer(I4), parameter :: MORPH_INTERIOR_SMOOTH = 2   !> this module

  !> The solver's controls, with the defaults the nodeshift law carries.
  type :: morph_smooth_ctl_t
    real(R8)    :: tol   = 1.0e-12_R8  !> stop when a sweep's largest update <= tol * max |disp|
    integer(I4) :: maxit = 100000      !> sweeps per solved set before MORPH_ERR_SMOOTH
    real(R8)    :: omega = 0.0_R8      !> over-relaxation; 0 = 2/(1 + sin(pi/N)), N the extent of the set
  end type morph_smooth_ctl_t

  !> What the last call did, for the caller's log.
  type :: morph_smooth_info_t
    integer(I4) :: sweeps_face     = 0        !> summed over the four side faces
    integer(I4) :: sweeps_interior = 0
    real(R8)    :: update_face     = 0.0_R8   !> the last sweep's largest update, relative to max |disp|
    real(R8)    :: update_interior = 0.0_R8
    real(R8)    :: omega_face      = 0.0_R8   !> the relaxation used (the largest over the faces)
    real(R8)    :: omega_interior  = 0.0_R8
  end type morph_smooth_info_t

contains

  !> Solve the four side faces and the interior of the block whose face `side`
  !> (1 = low index, 2 = high) of direction `dir` recedes. On entry dsp holds
  !> the displacement of every node: the prescribed values on the receding and
  !> the pinned face and on the four line edges, which are not touched, and a
  !> start for every other node. On exit the side faces and the interior hold
  !> the spring equilibrium. `scale` is the largest prescribed displacement
  !> (the unit of the stopping test); with scale = 0 nothing is solved.
  subroutine Morph_Smooth_Lineal ( node_old, dim, dir, side, scale, dsp, ctl, info, status )
    type(morph_vec3_t),        intent(in)    :: node_old(0:,0:,0:)
    integer(I4),               intent(in)    :: dim(3), dir, side
    real(R8),                  intent(in)    :: scale
    type(morph_vec3_t),        intent(inout) :: dsp(0:,0:,0:)
    type(morph_smooth_ctl_t),  intent(in)    :: ctl
    type(morph_smooth_info_t), intent(out)   :: info
    type(morph_status_t),      intent(out)   :: status
    ! Local
    integer(I4) :: t1, t2, f, lo(3), hi(3), sweeps
    logical     :: stencil(3)
    real(R8)    :: upd, omega

    call status%clear()
    info = morph_smooth_info_t()

    if ( dir < 1 .or. dir > 3 .or. side < 1 .or. side > 2 ) then
      call status%fail ( MORPH_ERR_SHAPE, 'smoothing: dir must be 1..3 and side 1 or 2' )
      return
    endif
    if ( any( shape(dsp) /= dim + 1 ) .or. any( shape(node_old) /= dim + 1 ) ) then
      call status%fail ( MORPH_ERR_SHAPE, 'smoothing: the arrays do not match the block' )
      return
    endif
    if ( scale <= 0.0_R8 ) return

    select case ( dir )
      case (1) ; t1 = 2 ; t2 = 3
      case (2) ; t1 = 1 ; t2 = 3
      case default ; t1 = 1 ; t2 = 2
    end select

    ! ---- the four side faces: a 2-D set each, the in-face springs only
    do f = 1, 4
      lo = 1 ; hi = dim - 1                   ! interior along the line and along the face
      stencil = .true.
      select case ( f )
        case (1) ; lo(t1) = 0       ; hi(t1) = 0       ; stencil(t1) = .false.
        case (2) ; lo(t1) = dim(t1) ; hi(t1) = dim(t1) ; stencil(t1) = .false.
        case (3) ; lo(t2) = 0       ; hi(t2) = 0       ; stencil(t2) = .false.
        case (4) ; lo(t2) = dim(t2) ; hi(t2) = dim(t2) ; stencil(t2) = .false.
      end select
      call sor_box ( node_old, dsp, lo, hi, stencil, ctl, scale, sweeps, upd, omega, status )
      if ( .not. status%ok() ) return
      info%sweeps_face = info%sweeps_face + sweeps
      info%update_face = max( info%update_face, upd )
      info%omega_face  = max( info%omega_face, omega )
    enddo

    ! ---- the interior: the six springs
    lo = 1 ; hi = dim - 1
    stencil = .true.
    call sor_box ( node_old, dsp, lo, hi, stencil, ctl, scale, sweeps, upd, omega, status )
    if ( .not. status%ok() ) return
    info%sweeps_interior = sweeps
    info%update_interior = upd
    info%omega_interior  = omega

  end subroutine Morph_Smooth_Lineal


  ! ---------------------------------------------------------------------------

  !> Red-black SOR on the box lo..hi of nodes, every node of which has both its
  !> neighbours in each direction of `stencil` (the box excludes the Dirichlet
  !> layer in those directions). An empty box is a no-op.
  subroutine sor_box ( node_old, dsp, lo, hi, stencil, ctl, scale, sweeps, upd, omega, status )
    type(morph_vec3_t),       intent(in)    :: node_old(0:,0:,0:)
    type(morph_vec3_t),       intent(inout) :: dsp(0:,0:,0:)
    integer(I4),              intent(in)    :: lo(3), hi(3)
    logical,                  intent(in)    :: stencil(3)
    type(morph_smooth_ctl_t), intent(in)    :: ctl
    real(R8),                 intent(in)    :: scale
    integer(I4),              intent(out)   :: sweeps
    real(R8),                 intent(out)   :: upd
    real(R8),                 intent(out)   :: omega
    type(morph_status_t),     intent(inout) :: status
    ! Local
    real(R8), allocatable :: w(:,:,:,:)      ! (6, i, j, k): -i, +i, -j, +j, -k, +k
    integer(I4) :: i, j, k, colour, nseg, d, bad(3), par(3)
    real(R8)    :: acc(3), sw, nw, ratio
    logical     :: zero
    character(len=MORPH_MSGLEN) :: msg

    sweeps = 0 ; upd = 0.0_R8 ; omega = 0.0_R8
    if ( any( hi < lo ) .or. .not. any( stencil ) ) return

    ! the relaxation: Young's optimum for the five-point Laplacian on the
    ! longest extent of the set (N intervals, N-1 unknowns), a good guess here
    nseg = 0
    do d = 1, 3
      if ( stencil(d) ) nseg = max( nseg, hi(d) - lo(d) + 2 )
    enddo
    if ( ctl%omega > 0.0_R8 ) then
      omega = ctl%omega
    else
      omega = 2.0_R8 / ( 1.0_R8 + sin( acos(-1.0_R8) / real( nseg, R8 ) ) )
    endif

    ! the spring stiffnesses, once: 1/|e| on node_old
    allocate( w(6, lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) )
    zero = .false.
    bad  = 0
    !$omp parallel do collapse(3) private(i,j,k) reduction(.or.:zero)
    do k = lo(3), hi(3)
    do j = lo(2), hi(2)
    do i = lo(1), hi(1)
      w(:,i,j,k) = 0.0_R8
      if ( stencil(1) ) then
        w(1,i,j,k) = stiffness( node_old(i,j,k)%c, node_old(i-1,j,k)%c )
        w(2,i,j,k) = stiffness( node_old(i,j,k)%c, node_old(i+1,j,k)%c )
      endif
      if ( stencil(2) ) then
        w(3,i,j,k) = stiffness( node_old(i,j,k)%c, node_old(i,j-1,k)%c )
        w(4,i,j,k) = stiffness( node_old(i,j,k)%c, node_old(i,j+1,k)%c )
      endif
      if ( stencil(3) ) then
        w(5,i,j,k) = stiffness( node_old(i,j,k)%c, node_old(i,j,k-1)%c )
        w(6,i,j,k) = stiffness( node_old(i,j,k)%c, node_old(i,j,k+1)%c )
      endif
      if ( any( w(:,i,j,k) < 0.0_R8 ) ) zero = .true.
    enddo; enddo; enddo
    !$omp end parallel do
    if ( zero ) then
      outer: do k = lo(3), hi(3)
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        if ( any( w(:,i,j,k) < 0.0_R8 ) ) then
          bad = [ i, j, k ]
          exit outer
        endif
      enddo; enddo; enddo outer
      call status%fail ( MORPH_ERR_SMOOTH, 'smoothing: two neighbouring nodes coincide (a spring of zero length)', &
                         bad(1), bad(2), bad(3) )
      return
    endif

    ! the sweeps. A node's colour is the parity of its indices in the stencil
    ! directions only, so that the two node planes of a 2-D mesh, solved as
    ! two faces, sweep their colours in the same order and come out bit for
    ! bit the same (with k in the parity they would not).
    par = 0
    do d = 1, 3
      if ( stencil(d) ) par(d) = 1
    enddo
    do while ( sweeps < ctl%maxit )
      sweeps = sweeps + 1
      upd = 0.0_R8
      do colour = 0, 1
        !$omp parallel do collapse(3) private(i,j,k,acc,sw,nw) reduction(max:upd)
        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          if ( mod( par(1) * i + par(2) * j + par(3) * k, 2 ) == colour ) then
            acc = 0.0_R8 ; sw = 0.0_R8
            if ( stencil(1) ) then
              acc = acc + w(1,i,j,k) * dsp(i-1,j,k)%c + w(2,i,j,k) * dsp(i+1,j,k)%c
              sw  = sw  + w(1,i,j,k) + w(2,i,j,k)
            endif
            if ( stencil(2) ) then
              acc = acc + w(3,i,j,k) * dsp(i,j-1,k)%c + w(4,i,j,k) * dsp(i,j+1,k)%c
              sw  = sw  + w(3,i,j,k) + w(4,i,j,k)
            endif
            if ( stencil(3) ) then
              acc = acc + w(5,i,j,k) * dsp(i,j,k-1)%c + w(6,i,j,k) * dsp(i,j,k+1)%c
              sw  = sw  + w(5,i,j,k) + w(6,i,j,k)
            endif
            acc = omega * ( acc / sw - dsp(i,j,k)%c )     ! the over-relaxed Gauss-Seidel update
            nw  = maxval( abs( acc ) )
            upd = max( upd, nw )
            dsp(i,j,k)%c = dsp(i,j,k)%c + acc
          endif
        enddo; enddo; enddo
        !$omp end parallel do
      enddo
      if ( upd <= ctl%tol * scale ) exit
    enddo

    ratio = upd / scale
    upd   = ratio
    if ( ratio > ctl%tol ) then
      write(msg,'(A,I0,A,ES9.2,A,ES9.2,A)') 'smoothing: not converged in ', sweeps, &
            ' sweeps: the last update is ', ratio, ' of the displacement, tolerance ', ctl%tol, &
            ' (raise smooth-maxit or the tolerance)'
      call status%fail ( MORPH_ERR_SMOOTH, trim(msg), lo(1), lo(2), lo(3) )
    endif

  contains

    !> 1/|q - p|; -1 flags a zero length (the caller reports it)
    pure function stiffness ( p, q ) result ( k )
      real(R8), intent(in) :: p(3), q(3)
      real(R8) :: k, e(3), l2
      e  = q - p
      l2 = dot_product( e, e )
      if ( l2 > 0.0_R8 ) then
        k = 1.0_R8 / sqrt( l2 )
      else
        k = -1.0_R8
      endif
    end function stiffness

  end subroutine sor_box

end module Morph_Smooth
