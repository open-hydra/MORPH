!>@brief Unit test of the nodeshift law's smoothed interior (Morph_Smooth,
!> interior = MORPH_INTERIOR_SMOOTH), always beside the line law on a copy of
!> the same mesh.
!>
!> A Cartesian 8 x 6 x 4 block, a graded one, a 2-D one (8 x 6 x 1):
!>   * a uniform displacement of the i-high face: the smoothed interior equals
!>     the line law to the solver's tolerance (on a straight line the linear
!>     blend IS the spring equilibrium), on the Cartesian and on the graded
!>     block (the 1/|e| weights keep the grading); the receding face is the
!>     line law's bit for bit, the pinned face is unchanged, the GCL holds;
!>   * a displacement that varies across the i-low face (a bump): the field of
!>     a step satisfies the spring equilibrium at every solved node, lies
!>     between the prescribed values (the maximum principle), moves only the
!>     coordinate the face moves in, and differs from the line law; the second
!>     step starts from the first and needs no more sweeps;
!>   * 2-D: the two node planes in k hold the same values, bit for bit;
!>   * a zero displacement leaves the mesh and the swept volumes bit-identical;
!>   * a displacement equal to the line length returns MORPH_ERR_NEGVOL with the
!>     line law's message; one sweep allowed returns MORPH_ERR_SMOOTH and leaves
!>     the nodes where they were; two coincident nodes return MORPH_ERR_SMOOTH.
!> The law never stops the program. The sweep counts and the checksums printed
!> must be the same at 1 and 8 threads (run_morph_guards.sh units, ctest
!> unit_*_threads_identical).
program test_morph_smooth
  use iso_fortran_env, only: R8 => real64
  use Morph_API

  implicit none

  integer,  parameter :: NI = 8, NJ = 6, NK = 4
  real(R8), parameter :: LX = 0.08_R8, LY = 0.03_R8, LZ = 0.02_R8
  real(R8), parameter :: TOL_GCL = 1.0e-14_R8
  real(R8), parameter :: DELTA   = 1.0e-3_R8       ! the face displacement per step [m]
  real(R8), parameter :: PI      = 3.14159265358979323846_R8

  type(morph_vec3_t), allocatable :: node(:,:,:), node_old(:,:,:), node0(:,:,:), nodel(:,:,:), nodel_old(:,:,:)
  type(morph_geom_t)   :: geom, geoml
  type(morph_status_t) :: st, stl
  type(morph_motion_nodeshift_t) :: law, lawl
  integer  :: dim(3), nfail, s, i, j, k, sweeps1, sweeps2
  real(R8) :: r_abs, r_rel, worst, dmax, t, res, dlo, dhi, expct
  logical  :: ok

  nfail = 0
  t = 0.0_R8
  write(*,'(A)') '== Morph nodeshift law, smoothed interior =='

  ! ---- 1. Cartesian block, i-high face recedes uniformly along -x: the line law
  call new_block ( NI, NJ, NK, graded = .false. )
  call set_laws ( side = 2, profile = 0 )
  call run ( 20 )
  call check ( ok, 'uniform, Cartesian: every update ok, both laws' )
  call check ( worst <= TOL_GCL, 'uniform, Cartesian: GCL residual' )
  write(*,'(A,ES10.2)') '   max relative GCL residual: ', worst
  call check ( same_face( NI ), 'uniform, Cartesian: the receding face is the line law''s, bit for bit' )
  call check ( same_face_as_start( 0 ), 'uniform, Cartesian: the pinned face i = 0 is unchanged' )
  dmax = max_diff()
  call check ( dmax <= 1.0e-13_R8, 'uniform, Cartesian: the smoothed interior is the line law within 1e-13 m' )
  write(*,'(A,ES10.2,A,I0,A,I0)') '   max |smooth - line| ', dmax, ' m; last step''s sweeps: faces ', &
        law%smooth_info%sweeps_face, ', interior ', law%smooth_info%sweeps_interior

  ! ---- 2. graded block (ratio 0.1 along the lines), the same
  call new_block ( NI, NJ, NK, graded = .true. )
  call set_laws ( side = 2, profile = 0 )
  call run ( 20 )
  call check ( ok, 'uniform, graded: every update ok, both laws' )
  call check ( worst <= TOL_GCL, 'uniform, graded: GCL residual' )
  call check ( same_face( NI ), 'uniform, graded: the receding face is the line law''s, bit for bit' )
  dmax = max_diff()
  call check ( dmax <= 1.0e-13_R8, 'uniform, graded: the smoothed interior is the line law within 1e-13 m (grading kept)' )
  write(*,'(A,ES10.2)') '   max |smooth - line| ', dmax

  ! ---- 3. a bump across the i-low face, along +x
  call new_block ( NI, NJ, NK, graded = .false. )
  call set_laws ( side = 1, profile = 1 )
  call Morph_Update ( node, node_old, dim, t, 1.0_R8, law, geom, st )
  sweeps1 = law%smooth_info%sweeps_interior
  call Morph_Update ( node, node_old, dim, t, 1.0_R8, law, geom, st )
  sweeps2 = law%smooth_info%sweeps_interior
  call check ( sweeps2 <= sweeps1, 'bump: the second step starts from the first and needs no more sweeps' )
  write(*,'(A,I0,A,I0)') '   interior sweeps: first step ', sweeps1, ', second step ', sweeps2
  call Morph_Update ( nodel, nodel_old, dim, t, 1.0_R8, lawl, geoml, stl )
  call Morph_Update ( nodel, nodel_old, dim, t, 1.0_R8, lawl, geoml, stl )
  call run ( 8 )
  call check ( ok, 'bump: every update ok, both laws' )
  call check ( worst <= TOL_GCL, 'bump: GCL residual' )
  write(*,'(A,ES10.2)') '   max relative GCL residual: ', worst
  ok = .true.
  do k = 0, NK ; do j = 0, NJ
    expct = node0(0,j,k)%c(1)
    do s = 1, 10
      expct = expct + law%disp(j,k)%c(1)
    enddo
    if ( node(0,j,k)%c(1) /= expct ) ok = .false.
  enddo ; enddo
  call check ( ok, 'bump: the receding face moved by its displacement, 10 times, bit for bit' )
  call check ( same_face_as_start( NI ), 'bump: the pinned face i = NI is unchanged' )
  ok = .true.
  do k = 0, NK ; do j = 0, NJ ; do i = 0, NI
    if ( any( node(i,j,k)%c(2:3) /= node0(i,j,k)%c(2:3) ) ) ok = .false.
  enddo ; enddo ; enddo
  call check ( ok, 'bump: y and z did not move, bit for bit (the face moves in x only)' )
  res = equilibrium()
  call check ( res <= 10.0_R8 * law%smooth%tol, 'bump: the last step''s field is the spring equilibrium (residual <= 10 tol)' )
  write(*,'(A,ES10.2,A,ES10.2)') '   max equilibrium residual, relative to the displacement: ', res, ', tolerance ', law%smooth%tol
  dlo = 0.0_R8 ; dhi = 0.0_R8
  do k = 0, NK ; do j = 0, NJ
    dlo = min( dlo, law%disp(j,k)%c(1) ) ; dhi = max( dhi, law%disp(j,k)%c(1) )
  enddo ; enddo
  ok = .true.
  do k = 0, NK ; do j = 0, NJ ; do i = 0, NI
    if ( node(i,j,k)%c(1) - node_old(i,j,k)%c(1) < dlo - 1.0e-15_R8 .or. &
         node(i,j,k)%c(1) - node_old(i,j,k)%c(1) > dhi + 1.0e-15_R8 ) ok = .false.
  enddo ; enddo ; enddo
  call check ( ok, 'bump: every node''s displacement lies between the prescribed ones (maximum principle)' )
  dmax = max_diff()
  call check ( dmax > 1.0e-4_R8, 'bump: the smoothed interior differs from the line law (the bump spread into the block)' )
  write(*,'(A,ES10.2,A)') '   max |smooth - line| ', dmax, ' m'
  write(*,'(A,ES23.16,A,ES23.16,A,ES23.16)') '   node(4,3,2): ', node(4,3,2)%c(1), ' ', node(4,3,2)%c(2), ' ', node(4,3,2)%c(3)
  write(*,'(A,ES23.16)') '   checksum (sum of all coordinates): ', checksum()

  ! ---- 4. one sweep allowed: a status, the nodes where they were
  law%smooth%maxit = 1
  call Morph_Update ( node, node_old, dim, t, 1.0_R8, law, geom, st )
  call check ( st%code == MORPH_ERR_SMOOTH, 'maxit = 1 on the bump returns MORPH_ERR_SMOOTH' )
  call check ( index( st%message, 'not converged in 1 sweeps' ) > 0, 'its message says so' )
  write(*,'(A,I0,A,A)') '   status ', st%code, ': ', trim(st%message)
  ok = .true.
  do k = 0, NK ; do j = 0, NJ ; do i = 0, NI
    if ( any( node(i,j,k)%c /= node_old(i,j,k)%c ) ) ok = .false.
  enddo ; enddo ; enddo
  call check ( ok, 'the failed step left every node where it was' )
  law%smooth%maxit = 100000

  ! ---- 5. two neighbouring nodes coincide: a status
  node(4,3,2) = node(5,3,2)
  call Morph_Update ( node, node_old, dim, t, 1.0_R8, law, geom, st )
  call check ( st%code == MORPH_ERR_SMOOTH .and. index( st%message, 'coincide' ) > 0, &
               'two coincident interior nodes return MORPH_ERR_SMOOTH' )
  write(*,'(A,I0,A,A,A,3(I0,1X))') '   status ', st%code, ': ', trim(st%message), ' at ', st%i, st%j, st%k

  ! ---- 6. 2-D block (two node planes in k), a bump in j
  call new_block ( NI, NJ, 1, graded = .false. )
  call set_laws ( side = 1, profile = 2 )
  call run ( 10 )
  call check ( ok, '2-D bump: every update ok, both laws' )
  call check ( worst <= TOL_GCL, '2-D bump: GCL residual' )
  ok = .true.
  do j = 0, NJ ; do i = 0, NI
    if ( any( node(i,j,0)%c(1:2) /= node(i,j,1)%c(1:2) ) ) ok = .false.
  enddo ; enddo
  call check ( ok, '2-D bump: the two node planes hold the same x and y, bit for bit' )
  res = equilibrium()
  call check ( res <= 10.0_R8 * law%smooth%tol, '2-D bump: the last step''s field is the spring equilibrium' )
  dmax = max_diff()
  call check ( dmax > 1.0e-4_R8, '2-D bump: the smoothed interior differs from the line law' )
  write(*,'(A,ES10.2,A,ES10.2,A,I0)') '   max equilibrium residual ', res, ', max |smooth - line| ', dmax, &
        ' m, last step''s face sweeps ', law%smooth_info%sweeps_face
  write(*,'(A,ES23.16)') '   checksum (sum of all coordinates): ', checksum()

  ! ---- 7. zero displacement: bitwise inert, exact zero swept volumes
  node0 = node
  do k = 0, 1 ; do j = 0, NJ
    law%disp(j,k)%c = 0.0_R8
  enddo ; enddo
  call Morph_Update ( node, node_old, dim, t, 1.0_R8, law, geom, st )
  call check ( st%ok(), 'zero displacement: update ok' )
  ok = .true.
  do k = 0, 1 ; do j = 0, NJ ; do i = 0, NI
    if ( any( node(i,j,k)%c /= node0(i,j,k)%c ) ) ok = .false.
  enddo ; enddo ; enddo
  call check ( ok .and. all( geom%dV_swept == 0.0_R8 ), &
               'zero displacement leaves the mesh and the swept volumes exactly unchanged' )

  ! ---- 8. a displacement equal to the line length: the line law's status and message
  do k = 0, 1 ; do j = 0, NJ
    law%disp(j,k)%c  = [ node(NI,j,k)%c(1) - node(0,j,k)%c(1), 0.0_R8, 0.0_R8 ]
    lawl%disp(j,k)%c = law%disp(j,k)%c
  enddo ; enddo
  call Morph_Update ( node,  node_old,  dim, t, 1.0_R8, law,  geom,  st  )
  call Morph_Update ( nodel, nodel_old, dim, t, 1.0_R8, lawl, geoml, stl )
  call check ( st%code == MORPH_ERR_NEGVOL .and. stl%code == MORPH_ERR_NEGVOL, &
               'disp = L_line returns MORPH_ERR_NEGVOL under both interiors' )
  call check ( st%message == stl%message .and. st%i == stl%i .and. st%j == stl%j .and. st%k == stl%k, &
               'with the same message and location' )
  write(*,'(A,I0,A,A)') '   consumed line: status ', st%code, ', ', trim(st%message)

  write(*,*)
  if ( nfail == 0 ) then
    write(*,'(A)') 'ALL SMOOTH CHECKS PASSED'
  else
    write(*,'(I0,A)') nfail, ' SMOOTH CHECK(S) FAILED'
    error stop 1
  endif

contains

  subroutine check ( cond, what )
    logical,          intent(in) :: cond
    character(len=*), intent(in) :: what
    if ( cond ) then
      write(*,'(A,A)') '   PASS ', what
    else
      write(*,'(A,A)') '   FAIL ', what
      nfail = nfail + 1
    endif
  end subroutine check

  !> A fresh block of ni x nj x nk cells, the same mesh in `node` (the smoothed
  !> law) and `nodel` (the line law), both initialised.
  subroutine new_block ( ni, nj, nk, graded )
    integer, intent(in) :: ni, nj, nk
    logical, intent(in) :: graded
    real(R8) :: x(0:ni), r, d1
    integer  :: ii, jj, kk
    dim = [ ni, nj, nk ]
    if ( allocated(node) ) deallocate( node, node_old, node0, nodel, nodel_old )
    allocate( node(0:ni,0:nj,0:nk), node_old(0:ni,0:nj,0:nk), node0(0:ni,0:nj,0:nk), &
              nodel(0:ni,0:nj,0:nk), nodel_old(0:ni,0:nj,0:nk) )
    if ( graded ) then
      r  = 0.1_R8 ** ( 1.0_R8 / ( ni - 1 ) )
      d1 = LX * ( 1.0_R8 - r ) / ( 1.0_R8 - r ** ni )
      x(0) = 0.0_R8
      do ii = 1, ni
        x(ii) = x(ii-1) + d1 * r ** ( ii - 1 )
      enddo
      x(ni) = LX
    else
      do ii = 0, ni
        x(ii) = LX * ii / ni
      enddo
    endif
    do kk = 0, nk ; do jj = 0, nj ; do ii = 0, ni
      node(ii,jj,kk)%c = [ x(ii), LY * jj / nj, LZ * kk / nk ]
    enddo ; enddo ; enddo
    node0 = node
    nodel = node
    call Morph_Init ( node,  dim, geom,  st  )
    call Morph_Init ( nodel, dim, geoml, stl )
    call check ( st%ok() .and. stl%ok(), 'Morph_Init' )
  end subroutine new_block

  !> Both laws on the i face `side`, with the displacement profile over the
  !> face's node indices (j, k): 0 uniform along -x (the i-high face), 1 a
  !> bump 1 + sin(pi j/nj) sin(pi k/nk) along +x, 2 a bump in j only (2-D).
  subroutine set_laws ( side, profile )
    integer, intent(in) :: side, profile
    integer  :: jj, kk
    real(R8) :: g
    law%dir  = 1 ; law%side  = side ; law%interior  = MORPH_INTERIOR_SMOOTH
    lawl%dir = 1 ; lawl%side = side ; lawl%interior = MORPH_INTERIOR_LINE
    if ( allocated(law%disp) )  deallocate( law%disp )
    if ( allocated(lawl%disp) ) deallocate( lawl%disp )
    if ( allocated(law%dsp) )   deallocate( law%dsp )
    allocate( law%disp(0:dim(2), 0:dim(3)), lawl%disp(0:dim(2), 0:dim(3)) )
    do kk = 0, dim(3) ; do jj = 0, dim(2)
      select case ( profile )
        case (0) ; law%disp(jj,kk)%c = [ -DELTA, 0.0_R8, 0.0_R8 ]
        case (1)
          g = sin( PI * jj / dim(2) ) * sin( PI * kk / dim(3) )
          law%disp(jj,kk)%c = [ DELTA * ( 1.0_R8 + g ), 0.0_R8, 0.0_R8 ]
        case default
          g = sin( PI * jj / dim(2) )
          law%disp(jj,kk)%c = [ DELTA * ( 1.0_R8 + g ), 0.0_R8, 0.0_R8 ]
      end select
      lawl%disp(jj,kk)%c = law%disp(jj,kk)%c
    enddo ; enddo
  end subroutine set_laws

  !> `steps` updates of both laws; ok = every status ok, worst = the smoothed
  !> law's largest relative GCL residual.
  subroutine run ( steps )
    integer, intent(in) :: steps
    integer :: ss
    ok = .true. ; worst = 0.0_R8
    do ss = 1, steps
      call Morph_Update ( node,  node_old,  dim, t, 1.0_R8, law,  geom,  st  )
      call Morph_Update ( nodel, nodel_old, dim, t, 1.0_R8, lawl, geoml, stl )
      ok = ok .and. st%ok() .and. stl%ok()
      call Morph_GCL_Residual ( geom, dim, r_abs, r_rel )
      worst = max( worst, r_rel )
    enddo
  end subroutine run

  !> the i face `ii` of the two laws holds the same bits
  logical function same_face ( ii )
    integer, intent(in) :: ii
    integer :: jj, kk
    same_face = .true.
    do kk = 0, dim(3) ; do jj = 0, dim(2)
      if ( any( node(ii,jj,kk)%c /= nodel(ii,jj,kk)%c ) ) same_face = .false.
    enddo ; enddo
  end function same_face

  !> the i face `ii` of the smoothed law holds the bits it started with
  logical function same_face_as_start ( ii )
    integer, intent(in) :: ii
    integer :: jj, kk
    same_face_as_start = .true.
    do kk = 0, dim(3) ; do jj = 0, dim(2)
      if ( any( node(ii,jj,kk)%c /= node0(ii,jj,kk)%c ) ) same_face_as_start = .false.
    enddo ; enddo
  end function same_face_as_start

  !> the largest |node - nodel| over every node and coordinate
  real(R8) function max_diff ()
    integer :: ii, jj, kk
    max_diff = 0.0_R8
    do kk = 0, dim(3) ; do jj = 0, dim(2) ; do ii = 0, dim(1)
      max_diff = max( max_diff, maxval( abs( node(ii,jj,kk)%c - nodel(ii,jj,kk)%c ) ) )
    enddo ; enddo ; enddo
  end function max_diff

  real(R8) function checksum ()
    integer :: ii, jj, kk
    checksum = 0.0_R8
    do kk = 0, dim(3) ; do jj = 0, dim(2) ; do ii = 0, dim(1)
      checksum = checksum + sum( node(ii,jj,kk)%c )
    enddo ; enddo ; enddo
  end function checksum

  !> The spring equilibrium of the last step's field node - node_old, for the
  !> lines along i: at every solved node (the four side faces with their in-face
  !> springs, the interior with all six) |sum_q w (dsp_q - dsp_p)| / sum_q w,
  !> the weights 1/|e| on node_old, maximised over the nodes and the
  !> coordinates and divided by the largest |dsp|.
  real(R8) function equilibrium ()
    integer  :: ii, jj, kk, dd, sgn, q(3)
    real(R8) :: acc(3), sw, w, e(3), dsp_p(3), dsp_q(3), sc
    logical  :: onj, onk, use(3)
    equilibrium = 0.0_R8 ; sc = 0.0_R8
    do kk = 0, dim(3) ; do jj = 0, dim(2) ; do ii = 0, dim(1)
      sc = max( sc, maxval( abs( node(ii,jj,kk)%c - node_old(ii,jj,kk)%c ) ) )
    enddo ; enddo ; enddo
    do kk = 0, dim(3) ; do jj = 0, dim(2) ; do ii = 1, dim(1) - 1
      onj = ( jj == 0 .or. jj == dim(2) )
      onk = ( kk == 0 .or. kk == dim(3) )
      if ( onj .and. onk ) cycle                         ! an edge: prescribed
      use = [ .true., .not. onj, .not. onk ]
      dsp_p = node(ii,jj,kk)%c - node_old(ii,jj,kk)%c
      acc = 0.0_R8 ; sw = 0.0_R8
      do dd = 1, 3
        if ( .not. use(dd) ) cycle
        do sgn = -1, 1, 2
          q = [ ii, jj, kk ] ; q(dd) = q(dd) + sgn
          e = node_old(q(1),q(2),q(3))%c - node_old(ii,jj,kk)%c
          w = 1.0_R8 / sqrt( dot_product( e, e ) )
          dsp_q = node(q(1),q(2),q(3))%c - node_old(q(1),q(2),q(3))%c
          acc = acc + w * ( dsp_q - dsp_p )
          sw  = sw + w
        enddo
      enddo
      equilibrium = max( equilibrium, maxval( abs( acc ) ) / sw )
    enddo ; enddo ; enddo
    equilibrium = equilibrium / sc
  end function equilibrium

end program test_morph_smooth
