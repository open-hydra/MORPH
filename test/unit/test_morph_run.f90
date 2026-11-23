!>@brief Test of Morph_Run, the whole run in one call.
!>
!> It pins:
!>   * the run against the loop it replaces, written out here as the driver had
!>     it: the laws installed by hand, domain%init, then for each step the
!>     nodeshift displacement and domain%advance. The nodes, the worst residuals
!>     and every step's relative residual agree bit for bit, for every moving
!>     law and every prescribed parameter, the taper included;
!>   * the nodeshift face: faces 1 and 6, on the first and on the last block,
!>     give the direction and the side written here by hand;
!>   * the smoothed interior: run%interior and the three smoothing controls
!>     reach the law (against the loop with the law set by hand, bit for bit),
!>     and one sweep allowed stops step 1 with MORPH_ERR_SMOOTH;
!>   * the bump: the face displacement at the centre of the face (velocity
!>     plus bump), at a corner outside the bump (velocity alone) and at a node
!>     in between, with no transverse part;
!>   * the observer: step 0 once the run has started, then each step once it
!>     is taken, the refused one included and marked;
!>   * the tolerance: a residual equal to it is accepted, one above it stops the
!>     run after that step with MORPH_ERR_GCL and the residual's block;
!>   * the failures: in the start (report%step 0, the observer never called)
!>     and in a step (report%step k, no call for it), the block in status%b;
!>   * the input: an unknown law, a block or a face out of range, a smooth
!>     omega of 2 refused with MORPH_ERR_INPUT, the laws left as they were;
!>     the bounds accepted, the largest smooth omega below 2 among them;
!>   * zero steps: the laws installed, the observer called once, at step 0;
!>   * a static run after a moving one, without an observer: nothing moves.
module test_morph_run_recorder
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_API
  implicit none
  private
  public :: recorder_t

  integer(I4), parameter :: MAXCALL = 16

  !> Every call it gets: the step it is told, the domain's own step count, the
  !> verdict and the relative residual.
  type, extends(morph_run_observer_t) :: recorder_t
    integer(I4) :: n = 0
    integer(I4) :: step(MAXCALL)     = -1
    integer(I4) :: dstep(MAXCALL)    = -1
    logical     :: accepted(MAXCALL) = .true.
    real(R8)    :: gcl_rel(MAXCALL)  = -1.0_R8
  contains
    procedure :: on_step => recorder_on_step
    procedure :: reset   => recorder_reset
  end type recorder_t

contains

  subroutine recorder_on_step ( self, domain, info )
    class(recorder_t),      intent(inout) :: self
    type(morph_domain_t),   intent(in)    :: domain
    type(morph_run_step_t), intent(in)    :: info
    self%n = self%n + 1
    if ( self%n > MAXCALL ) return
    self%step(self%n)     = info%step
    self%dstep(self%n)    = domain%step
    self%accepted(self%n) = info%accepted
    self%gcl_rel(self%n)  = info%gcl_rel
  end subroutine recorder_on_step

  subroutine recorder_reset ( self )
    class(recorder_t), intent(inout) :: self
    self%n        = 0
    self%step     = -1
    self%dstep    = -1
    self%accepted = .true.
    self%gcl_rel  = -1.0_R8
  end subroutine recorder_reset

end module test_morph_run_recorder


program test_morph_run
  use iso_fortran_env, only: I4 => int32, I8 => int64, R8 => real64
  use Morph_API
  use test_morph_run_recorder, only: recorder_t
  implicit none
  integer(I4), parameter :: NB = 2
  integer(I4), parameter :: DIMS(3,NB) = reshape( [3,2,2, 2,2,2], [3,NB] )
  type(morph_domain_t)     :: dom, ref
  type(morph_run_t)        :: run, bad
  type(morph_run_report_t) :: rep
  type(morph_status_t)     :: st
  type(recorder_t)         :: rec
  real(R8)    :: worst_rel, worst_abs, gcl(8), centre, corner, edge
  logical     :: ok, installed, transverse

  ok = .true.
  write(*,'(A)') '== Morph one-call run test =='
  write(*,*)

  ! 1. static: time, steps, the observer's calls, the tolerance's bound
  call make_domain ( dom, .false. )
  call make_domain ( ref, .false. )
  run = morph_run_t()
  run%law = 'static' ; run%t0 = 1.0_R8 ; run%dt = 0.5_R8 ; run%steps = 3 ; run%gcl_tol = 0.0_R8
  call rec%reset()
  call Morph_Run ( dom, run, rep, st, rec )
  call expect ( st%ok() .and. rep%step == 3 .and. dom%step == 3, 'static: 3 steps, status ok' )
  call expect ( dom%t == ((1.0_R8 + 0.5_R8) + 0.5_R8) + 0.5_R8, 'static: t is t0 plus dt three times' )
  call expect ( same_nodes(dom, ref), 'static: every node bit for bit as given' )
  call expect ( rep%worst_rel == 0.0_R8 .and. rep%worst_abs == 0.0_R8, 'static: worst residuals exactly 0' )
  call expect ( rec%n == 4 .and. all( rec%step(1:4) == [0,1,2,3] ) .and. all( rec%dstep(1:4) == [0,1,2,3] ), &
                'observer: step 0 once started, then each step once taken' )
  call expect ( all( rec%accepted(1:4) ), 'tolerance: a residual equal to it (0 and 0) is accepted' )

  ! 2. translation, against the loop written out
  run = morph_run_t()
  run%law = 'translation' ; run%dt = 1.0e-3_R8 ; run%steps = 4 ; run%vel = [0.75_R8, -0.5_R8, 0.25_R8]
  call compare ( run, 'translation' )

  ! 3. prescribed, every parameter different, then tapered
  run = morph_run_t()
  run%law = 'prescribed' ; run%t0 = 0.1_R8 ; run%dt = 1.0e-3_R8 ; run%steps = 4
  run%amp = [2.0e-3_R8, 1.5e-3_R8, 1.0e-3_R8] ; run%omega = 50.0_R8
  run%kx  = [10.0_R8, 5.0_R8, 8.0_R8] ; run%ky = [7.0_R8, 12.0_R8, 3.0_R8] ; run%kz = [4.0_R8, 6.0_R8, 9.0_R8]
  call compare ( run, 'prescribed' )
  call expect ( rep%worst_rel > 0.0_R8, 'prescribed: a residual of round-off, not 0, so the worst is tested' )
  run%taper = .true.
  call compare ( run, 'prescribed, tapered' )

  ! 4. nodeshift: the last block's face 6, then block 1's face 1
  run = morph_run_t()
  run%law = 'nodeshift' ; run%dt = 1.0e-3_R8 ; run%steps = 3
  run%shift_block = NB ; run%shift_face = 6 ; run%shift_vel = [0.0_R8, 0.0_R8, -2.5_R8]
  call compare ( run, 'nodeshift, last block, face 6' )
  run%shift_block = 1 ; run%shift_face = 1 ; run%shift_vel = [2.5_R8, 0.0_R8, 0.0_R8]
  call compare ( run, 'nodeshift, block 1, face 1' )

  ! 4b. the same face with the smoothed interior and controls of its own, then
  !     one sweep allowed
  run%interior = MORPH_INTERIOR_SMOOTH
  run%smooth_tol = 1.0e-9_R8 ; run%smooth_maxit = 5000 ; run%smooth_omega = 1.5_R8
  call compare ( run, 'nodeshift, block 1, face 1, smoothed interior' )
  call make_domain ( dom, .false. )
  run%smooth_maxit = 1
  call rec%reset()
  call Morph_Run ( dom, run, rep, st, rec )
  call expect ( st%code == MORPH_ERR_SMOOTH .and. st%b == 1 .and. rep%step == 1 .and. dom%step == 0, &
                'smoothed interior, one sweep allowed: MORPH_ERR_SMOOTH in step 1, block 1' )

  ! 5. a step above the tolerance
  call make_domain ( dom, .false. )
  run = morph_run_t()
  run%law = 'prescribed' ; run%dt = 1.0e-3_R8 ; run%steps = 3 ; run%gcl_tol = -1.0_R8
  run%amp = [2.0e-3_R8, 1.5e-3_R8, 1.0e-3_R8] ; run%omega = 50.0_R8
  run%kx  = [10.0_R8, 5.0_R8, 8.0_R8] ; run%ky = [7.0_R8, 12.0_R8, 3.0_R8] ; run%kz = [4.0_R8, 6.0_R8, 9.0_R8]
  call rec%reset()
  call Morph_Run ( dom, run, rep, st, rec )
  call expect ( st%code == MORPH_ERR_GCL .and. rep%step == 1 .and. dom%step == 1, &
                'tolerance: a step above it stops the run after that step, MORPH_ERR_GCL' )
  call expect ( rep%last%step == 1 .and. .not. rep%last%accepted .and. st%b == rep%last%gcl_block .and. st%b >= 1, &
                'tolerance: the refused step in the report, its block in status%b' )
  call expect ( rec%n == 2 .and. rec%step(2) == 1 .and. rec%accepted(1) .and. .not. rec%accepted(2), &
                'observer: the refused step is shown, marked refused' )

  ! 6. a step that fails: block 2's face 1 asked to recede by more than its lines
  call make_domain ( dom, .false. )
  run = morph_run_t()
  run%law = 'nodeshift' ; run%dt = 1.0e-3_R8 ; run%steps = 3
  run%shift_block = 2 ; run%shift_face = 1 ; run%shift_vel = [500.0_R8, 0.0_R8, 0.0_R8]
  call rec%reset()
  call Morph_Run ( dom, run, rep, st, rec )
  call expect ( st%code == MORPH_ERR_NEGVOL .and. st%b == 2 .and. rep%step == 1 .and. dom%step == 0, &
                'a step that fails: its status, block 2 in status%b, report%step 1' )
  call expect ( rec%n == 1 .and. rec%step(1) == 0, 'observer: not called for a step that failed' )

  ! 7. a start that fails: block 2 inverted
  call make_domain ( dom, .true. )
  run = morph_run_t()
  run%dt = 1.0_R8 ; run%steps = 2
  call rec%reset()
  call Morph_Run ( dom, run, rep, st, rec )
  call expect ( st%code == MORPH_ERR_NEGVOL .and. st%b == 2 .and. rep%step == 0, &
                'a start that fails: its status, block 2 in status%b, report%step 0' )
  call expect ( rec%n == 0, 'observer: never called when the start fails' )

  ! 8. input refused, after a translation run
  call make_domain ( dom, .false. )
  run = morph_run_t()
  run%law = 'translation' ; run%dt = 1.0e-3_R8 ; run%steps = 1 ; run%vel = [1.0_R8, 0.0_R8, 0.0_R8]
  call Morph_Run ( dom, run, rep, st )
  bad = morph_run_t()
  bad%law = 'spin'
  call refused ( bad, 'an unknown law' )
  bad%law = 'nodeshift' ; bad%dt = 1.0e-3_R8 ; bad%steps = 1 ; bad%shift_vel = [2.5_R8, 0.0_R8, 0.0_R8]
  bad%shift_face = 1
  bad%shift_block = 0
  call refused ( bad, 'nodeshift block 0' )
  bad%shift_block = NB + 1
  call refused ( bad, 'nodeshift block nblocks + 1' )
  bad%shift_block = 1
  bad%shift_face = 0
  call refused ( bad, 'nodeshift face 0' )
  bad%shift_face = 7
  call refused ( bad, 'nodeshift face 7' )
  bad%shift_face = 1
  bad%smooth_omega = 2.0_R8
  call refused ( bad, 'smooth omega 2' )

  ! 9. zero steps
  call make_domain ( dom, .false. )
  run = morph_run_t()
  run%law = 'nodeshift' ; run%t0 = 2.0_R8 ; run%dt = 1.0e-3_R8 ; run%steps = 0
  run%shift_block = 1 ; run%shift_face = 1 ; run%shift_vel = [2.5_R8, 0.0_R8, 0.0_R8]
  call rec%reset()
  call Morph_Run ( dom, run, rep, st, rec )
  call expect ( st%ok() .and. rep%step == 0 .and. dom%step == 0 .and. dom%t == 2.0_R8, &
                'zero steps: status ok, the domain started at t0' )
  call expect ( rec%n == 1 .and. rec%step(1) == 0, 'zero steps: the observer called once, at step 0' )
  installed = .false.
  select type ( law => dom%blk(1)%law )
  type is ( morph_motion_nodeshift_t )
    installed = allocated( law%disp )
  end select
  call expect ( installed .and. dom%blk(2)%law%name == 'static', 'zero steps: the laws installed all the same' )

  ! 9b. the bump on block 1's face 1 (3x3 nodes, the centre at (1,1)) and the
  !     largest smooth omega below 2, with zero steps
  call make_domain ( dom, .false. )
  run = morph_run_t()
  run%law = 'nodeshift' ; run%dt = 1.0e-3_R8 ; run%steps = 0
  run%shift_block = 1 ; run%shift_face = 1 ; run%shift_vel = [2.5_R8, 0.0_R8, 0.0_R8]
  run%shift_bump_vel = [2.0_R8, 0.0_R8, 0.0_R8] ; run%shift_bump_width = 0.6_R8
  run%smooth_omega = nearest( 2.0_R8, -1.0_R8 )
  call Morph_Run ( dom, run, rep, st )
  call expect ( st%ok(), 'smooth omega: the largest number below 2 accepted' )
  centre = -1.0_R8 ; corner = -1.0_R8 ; edge = -1.0_R8 ; transverse = .false.
  select type ( law => dom%blk(1)%law )
  type is ( morph_motion_nodeshift_t )
    centre = law%disp(1,1)%c(1)
    corner = law%disp(0,0)%c(1)
    edge   = law%disp(1,0)%c(1)
    transverse = all( law%disp(:,:)%c(2) == 0.0_R8 ) .and. all( law%disp(:,:)%c(3) == 0.0_R8 )
  end select
  call expect ( centre == (2.5_R8 + 2.0_R8) * 1.0e-3_R8, 'bump: the centre of the face by (vel + bump-vel)*dt' )
  call expect ( corner == 2.5_R8 * 1.0e-3_R8, 'bump: a corner, outside the bump, by vel*dt' )
  call expect ( edge > corner .and. edge < centre, 'bump: a node inside the bump, off the centre, in between' )
  call expect ( transverse, 'bump: no transverse displacement' )

  ! 10. static after translation, without an observer
  call make_domain ( dom, .false. )
  call make_domain ( ref, .false. )
  run = morph_run_t()
  run%law = 'translation' ; run%dt = 1.0e-3_R8 ; run%steps = 2 ; run%vel = [0.75_R8, -0.5_R8, 0.25_R8]
  call Morph_Run ( dom, run, rep, st )
  call Morph_Run ( ref, run, rep, st )
  run%law = 'static'
  call Morph_Run ( dom, run, rep, st )
  call expect ( st%ok() .and. dom%step == 2 .and. same_nodes(dom, ref), &
                'static after translation, no observer: nothing moves' )

  write(*,*)
  if ( ok ) then
    write(*,'(A)') '   RESULT: PASS'
  else
    write(*,'(A)') '   RESULT: FAIL'
    error stop 1
  endif

contains

  subroutine expect ( good, what )
    logical,          intent(in) :: good
    character(len=*), intent(in) :: what
    if ( good ) then
      write(*,'(A,A,A)') '   ', what, ': PASS'
    else
      write(*,'(A,A)') '   FAIL: ', what
      ok = .false.
    endif
  end subroutine expect


  !> Morph_Run on a fresh domain against the loop written out on another.
  subroutine compare ( r, what )
    type(morph_run_t), intent(in) :: r
    character(len=*),  intent(in) :: what
    logical     :: good
    integer(I4) :: n

    call make_domain ( dom, .false. )
    call make_domain ( ref, .false. )
    call rec%reset()
    call Morph_Run ( dom, r, rep, st, rec )
    call manual ( ref, r, worst_rel, worst_abs, gcl, good )
    n = r%steps
    call expect ( st%ok() .and. good .and. same_nodes(dom, ref), what//': nodes bit for bit as the loop written out' )
    call expect ( rep%worst_rel == worst_rel .and. rep%worst_abs == worst_abs, what//': worst residuals as the loop' )
    call expect ( rec%n == n + 1 .and. all( rec%gcl_rel(2:n+1) == gcl(1:n) ), &
                  what//': the observer sees every step''s residual' )
  end subroutine compare


  !> The input r refused: MORPH_ERR_INPUT before the start, no observer call,
  !> the translation laws of dom left in place.
  subroutine refused ( r, what )
    type(morph_run_t), intent(in) :: r
    character(len=*),  intent(in) :: what
    call rec%reset()
    call Morph_Run ( dom, r, rep, st, rec )
    call expect ( st%code == MORPH_ERR_INPUT .and. rep%step == 0 .and. rec%n == 0 .and. &
                  dom%blk(1)%law%name == 'translation' .and. dom%blk(2)%law%name == 'translation', &
                  'input refused, the laws left as they were: '//what )
  end subroutine refused


  !> The loop Morph_Run replaces, as the driver had it: the laws installed by
  !> hand (the direction and side of nodeshift faces 1 and 6 written out),
  !> domain%init, then for each step the nodeshift displacement and
  !> domain%advance; the worst residuals kept and each step's relative one in
  !> gcl_step.
  subroutine manual ( d, r, w_rel, w_abs, gcl_step, good )
    type(morph_domain_t), intent(inout) :: d
    type(morph_run_t),    intent(in)    :: r
    real(R8),             intent(out)   :: w_rel, w_abs, gcl_step(:)
    logical,              intent(out)   :: good
    ! Local
    type(morph_status_t) :: s
    integer(I4) :: b, ib, ib_abs, ib_rel, step, dir, side, t1, t2, ob, oi, oj, oc
    real(R8)    :: r_abs, r_rel, ratio

    select case ( trim(r%law) )
    case ( 'translation' )
      do b = 1, d%nblocks
        deallocate( d%blk(b)%law )
        allocate( d%blk(b)%law, source = morph_motion_translation_t( name='translation', vel=r%vel ) )
      enddo
    case ( 'prescribed' )
      do b = 1, d%nblocks
        deallocate( d%blk(b)%law )
        allocate( morph_motion_prescribed_t :: d%blk(b)%law )
        select type ( law => d%blk(b)%law )
        type is ( morph_motion_prescribed_t )
          law%amp   = r%amp
          law%kx    = r%kx
          law%ky    = r%ky
          law%kz    = r%kz
          law%omega = r%omega
          law%taper_to_boundary = r%taper
        end select
      enddo
    case ( 'nodeshift' )
      select case ( r%shift_face )
        case (1) ; dir = 1 ; side = 1 ; t1 = 2 ; t2 = 3
        case (6) ; dir = 3 ; side = 2 ; t1 = 1 ; t2 = 2
        case default ; error stop 'manual: nodeshift faces 1 and 6 only'
      end select
      b = r%shift_block
      deallocate( d%blk(b)%law )
      allocate( d%blk(b)%law, source = morph_motion_nodeshift_t( name='nodeshift', dir=dir, side=side ) )
      select type ( law => d%blk(b)%law )
      type is ( morph_motion_nodeshift_t )
        allocate( law%disp( 0:d%blk(b)%dim(t1), 0:d%blk(b)%dim(t2) ) )
        law%interior     = r%interior
        law%smooth%tol   = r%smooth_tol
        law%smooth%maxit = r%smooth_maxit
        law%smooth%omega = r%smooth_omega
      end select
    case default
      error stop 'manual: translation, prescribed and nodeshift only'
    end select

    call d%init ( r%t0, s, ib )
    good = s%ok()
    w_rel = 0.0_R8
    w_abs = 0.0_R8
    gcl_step = -1.0_R8
    do step = 1, r%steps
      if ( trim(r%law) == 'nodeshift' ) then
        select type ( law => d%blk(r%shift_block)%law )
        type is ( morph_motion_nodeshift_t )
          law%disp = morph_vec3_t( r%shift_vel * r%dt )
        end select
      endif
      call d%advance ( r%dt, s, ib )
      good = good .and. s%ok()
      call d%gcl_residual ( r_abs, r_rel, ib_abs, ib_rel )
      call d%outflow_ratio ( ratio, ob, oi, oj, oc )
      w_rel = max( w_rel, r_rel )
      w_abs = max( w_abs, r_abs )
      gcl_step(step) = r_rel
    enddo
  end subroutine manual


  !> Block 1 on x in [0, 0.3] (3x2x2 cells), block 2 on x in [0.3, 0.5]
  !> (2x2x2), y and z in [0, 0.2], every node moved by a smooth warp so that
  !> no grid line is straight. invert: block 2 mirrored in x, every one of its
  !> cells with a negative volume.
  subroutine make_domain ( d, invert )
    type(morph_domain_t), intent(inout) :: d
    logical,              intent(in)    :: invert
    ! Local
    type(morph_status_t) :: s
    integer(I4) :: b, i, j, k
    real(R8)    :: x, y, z

    call d%allocate_blocks ( DIMS, s )
    if ( .not. s%ok() ) error stop 'make_domain: allocate_blocks failed'
    do b = 1, NB
      do k = 0, DIMS(3,b)
        do j = 0, DIMS(2,b)
          do i = 0, DIMS(1,b)
            x = merge( 0.0_R8, 0.3_R8, b == 1 ) + 0.1_R8*i
            y = 0.1_R8*j
            z = 0.1_R8*k
            d%blk(b)%node(i,j,k)%c = [ x + 0.01_R8*sin(7.0_R8*y + 3.0_R8*z), &
                                       y + 0.01_R8*sin(5.0_R8*x + 2.0_R8*z), &
                                       z + 0.01_R8*sin(4.0_R8*x + 6.0_R8*y) ]
            if ( invert .and. b == 2 ) d%blk(b)%node(i,j,k)%c(1) = -d%blk(b)%node(i,j,k)%c(1)
          enddo
        enddo
      enddo
    enddo
  end subroutine make_domain


  !> Every node of every block holds the same bits in a and b.
  logical function same_nodes ( a, b ) result ( same )
    type(morph_domain_t), intent(in) :: a, b
    ! Local
    integer(I4) :: blk, i, j, k

    same = ( a%nblocks == b%nblocks )
    if ( .not. same ) return
    do blk = 1, a%nblocks
      do k = 0, a%blk(blk)%dim(3)
        do j = 0, a%blk(blk)%dim(2)
          do i = 0, a%blk(blk)%dim(1)
            if ( any( transfer(a%blk(blk)%node(i,j,k)%c, 0_I8, 3) /= &
                      transfer(b%blk(blk)%node(i,j,k)%c, 0_I8, 3) ) ) same = .false.
          enddo
        enddo
      enddo
    enddo
  end function same_nodes

end program test_morph_run
