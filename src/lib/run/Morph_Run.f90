!>@brief One call that moves a multi-block mesh for a number of steps: the law
!>       installed on every block, the domain started, then every step
!>       advanced, checked and shown to the caller.
!>
!> It is the whole computation of the standalone driver (src/app/main.f90 reads
!> the input, makes this call and writes), and any caller that wants MORPH to
!> hold and run the mesh can make the same call. Like the rest of the library it
!> never stops and does no I/O: a failure comes back as the status, and the
!> caller's observer (morph_run_observer_t) sees the mesh once the run has
!> started and after every step, to write what it wants.
!>
!> A step is accepted when its relative discrete-GCL residual, the largest over
!> all cells of all blocks (Morph_Domain_m), is not above run%gcl_tol, as FUSS's
!> Update_Mesh asserts it. A step above it is still shown to the observer, and
!> the run stops with MORPH_ERR_GCL.
!>
!> MORPH must never `use` a FUSS module -- see Morph_Types_m.
module Morph_Run_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m
  use Morph_Domain_m,           only: morph_domain_t
  use Morph_Motion_Static,      only: morph_motion_static_t
  use Morph_Motion_Translation, only: morph_motion_translation_t
  use Morph_Motion_Prescribed,  only: morph_motion_prescribed_t
  use Morph_Motion_Nodeshift,   only: morph_motion_nodeshift_t

  implicit none
  private

  public :: morph_run_t, morph_run_step_t, morph_run_report_t, morph_run_observer_t
  public :: Morph_Run

  !> What to run: one law on every block with its parameters, the steps, and
  !> the largest relative discrete-GCL residual a step may leave. Only numbers
  !> and the law's name; the defaults are the driver's (Register_Morph).
  type :: morph_run_t
    character(len=32) :: law = 'static'       !> static | translation | prescribed | nodeshift
    real(R8)    :: t0      = 0.0_R8           !> time of the mesh as given [s]
    real(R8)    :: dt      = 0.0_R8           !> step [s]
    integer(I4) :: steps   = 0                !> number of steps
    real(R8)    :: gcl_tol = 1.0e-10_R8       !> largest relative discrete-GCL residual accepted in a step
    ! translation
    real(R8)    :: vel(3)  = 0.0_R8           !> velocity of every node [m/s]
    ! prescribed (morph_motion_prescribed_t; its phases are the law's own)
    real(R8)    :: amp(3)  = 0.0_R8           !> displacement amplitude per coordinate [m]
    real(R8)    :: kx(3)   = 0.0_R8           !> wavenumbers multiplying x [1/m]
    real(R8)    :: ky(3)   = 0.0_R8           !> wavenumbers multiplying y [1/m]
    real(R8)    :: kz(3)   = 0.0_R8           !> wavenumbers multiplying z [1/m]
    real(R8)    :: omega   = 0.0_R8           !> angular frequency [1/s]
    logical     :: taper   = .false.          !> taper the motion to zero on the block boundary
    ! nodeshift: one block face recedes at a constant velocity, every other block is static
    integer(I4) :: shift_block  = 0           !> the block, 1..nblocks
    integer(I4) :: shift_face   = 0           !> the face, 1..6 = i-low, i-high, j-low, j-high, k-low, k-high
    real(R8)    :: shift_vel(3) = 0.0_R8      !> velocity of every node of that face [m/s]
  end type morph_run_t

  !> One step as the run saw it. Step 0 is the mesh as started, before the first step.
  type :: morph_run_step_t
    integer(I4) :: step      = 0
    real(R8)    :: gcl_abs   = 0.0_R8         !> largest |(V_new - V_old) - sum_f dV_swept(f)| over all cells [m^3]
    real(R8)    :: gcl_rel   = 0.0_R8         !> the same relative to the cell volume
    integer(I4) :: gcl_block = 0              !> the block of gcl_rel, 0 when it is zero
    real(R8)    :: outflow   = 0.0_R8         !> largest fraction of its own volume a cell swept out
    logical     :: accepted  = .true.         !> .false.: gcl_rel exceeded run%gcl_tol, the run stops after this step
  end type morph_run_step_t

  !> How the run ended. With the status ok, step = run%steps.
  type :: morph_run_report_t
    integer(I4)            :: step      = 0       !> where the run ended: 0 while starting, k in step k
    real(R8)               :: worst_rel = 0.0_R8  !> largest gcl_rel over the steps taken
    real(R8)               :: worst_abs = 0.0_R8  !> largest gcl_abs over the steps taken
    type(morph_run_step_t) :: last                !> the last step that moved the mesh, refused or not (step 0: none)
  end type morph_run_report_t

  !> What the caller does with the mesh while it runs (a log, a snapshot): it
  !> extends this type with its own state, and MORPH calls on_step.
  type, abstract :: morph_run_observer_t
  contains
    procedure(morph_run_on_step_i), deferred :: on_step
  end type morph_run_observer_t

  abstract interface
    !> Called with info%step = 0 once the run has started (the laws installed,
    !> the mesh checked by domain%init), then after every step, the refused one
    !> included, before the run goes on or stops. It cannot stop the run.
    subroutine morph_run_on_step_i ( self, domain, info )
      import :: morph_run_observer_t, morph_domain_t, morph_run_step_t
      class(morph_run_observer_t), intent(inout) :: self
      type(morph_domain_t),        intent(in)    :: domain
      type(morph_run_step_t),      intent(in)    :: info
    end subroutine morph_run_on_step_i
  end interface

contains

  !> Install run%law on every block, start the domain at run%t0 and advance it
  !> run%steps times by run%dt, showing each step to the observer.
  !>
  !> On a failure the run stops; the status says why, and status%b the block:
  !>   MORPH_ERR_INPUT  the law or its parameters; nothing was installed and
  !>                    report%step = 0
  !>   MORPH_ERR_GCL    step report%step was refused: its gcl_rel exceeded
  !>                    run%gcl_tol; the mesh is the one that step made
  !>   any other code   domain%init (report%step = 0) or domain%advance in step
  !>                    report%step, which leaves the domain part-way through
  !>                    that step (Morph_Domain_m)
  subroutine Morph_Run ( domain, run, report, status, observer )
    type(morph_domain_t),        intent(inout)           :: domain
    type(morph_run_t),           intent(in)              :: run
    type(morph_run_report_t),    intent(out)             :: report
    type(morph_status_t),        intent(out)             :: status
    class(morph_run_observer_t), intent(inout), optional :: observer
    ! Local
    type(morph_run_step_t)      :: info
    integer(I4)                 :: step, ib, ib_abs, ob, oi, oj, oc
    character(len=MORPH_MSGLEN) :: msg

    call status%clear()

    call install_laws ( domain, run, status )
    if ( .not. status%ok() ) return

    call domain%init ( run%t0, status, ib )
    if ( .not. status%ok() ) then
      status%b = ib
      return
    endif

    if ( present(observer) ) call observer%on_step ( domain, info )

    do step = 1, run%steps
      report%step = step

      call domain%advance ( run%dt, status, ib )
      if ( .not. status%ok() ) then
        status%b = ib
        return
      endif

      info%step = step
      call domain%gcl_residual ( info%gcl_abs, info%gcl_rel, ib_abs, info%gcl_block )
      call domain%outflow_ratio ( info%outflow, ob, oi, oj, oc )
      info%accepted = .not. ( info%gcl_rel > run%gcl_tol )
      report%worst_rel = max( report%worst_rel, info%gcl_rel )
      report%worst_abs = max( report%worst_abs, info%gcl_abs )
      report%last      = info

      if ( present(observer) ) call observer%on_step ( domain, info )

      if ( .not. info%accepted ) then
        write(msg,'(A,I0,A,ES10.3,A,I0,A,ES10.3)') 'run: step ', step, ': discrete GCL residual ', &
              info%gcl_rel, ' (relative, block ', info%gcl_block, ') exceeds the tolerance ', run%gcl_tol
        call status%fail ( MORPH_ERR_GCL, msg, b=info%gcl_block )
        return
      endif
    enddo

  end subroutine Morph_Run


  !> run%law on every block, one instance per block (the prescribed law caches
  !> the mesh of the block it first moves). The law and its parameters are
  !> checked first, so a refused run leaves the domain's laws as they were.
  !> Under nodeshift every other block is static, and the receding face gets
  !> its displacement once: run%shift_vel*run%dt at every node, the same every
  !> step, which the law reads and never consumes.
  subroutine install_laws ( domain, run, status )
    type(morph_domain_t), intent(inout) :: domain
    type(morph_run_t),    intent(in)    :: run
    type(morph_status_t), intent(inout) :: status
    ! Local
    integer(I4)                 :: b, d, t1, t2, a, c
    character(len=MORPH_MSGLEN) :: msg

    select case ( trim(run%law) )
    case ( 'static', 'translation', 'prescribed' )
    case ( 'nodeshift' )
      if ( run%shift_block < 1 .or. run%shift_block > domain%nblocks ) then
        write(msg,'(A,I0,A,I0)') 'run: nodeshift block ', run%shift_block, &
              ' is not a block of the domain, 1..', domain%nblocks
        call status%fail ( MORPH_ERR_INPUT, msg )
        return
      endif
      if ( run%shift_face < 1 .or. run%shift_face > 6 ) then
        write(msg,'(A,I0,A)') 'run: nodeshift face ', run%shift_face, &
              ' is not a face, 1..6 (i-low, i-high, j-low, j-high, k-low, k-high)'
        call status%fail ( MORPH_ERR_INPUT, msg )
        return
      endif
    case default
      call status%fail ( MORPH_ERR_INPUT, 'run: law '//trim(run%law)// &
                         ' is not known (static, translation, prescribed, nodeshift)' )
      return
    end select

    do b = 1, domain%nblocks
      if ( allocated(domain%blk(b)%law) ) deallocate( domain%blk(b)%law )
      select case ( trim(run%law) )
      case ( 'translation' )
        allocate( domain%blk(b)%law, source = morph_motion_translation_t( name='translation', vel=run%vel ) )
      case ( 'prescribed' )
        allocate( morph_motion_prescribed_t :: domain%blk(b)%law )
        select type ( law => domain%blk(b)%law )
        type is ( morph_motion_prescribed_t )
          law%name  = 'prescribed'
          law%amp   = run%amp
          law%kx    = run%kx
          law%ky    = run%ky
          law%kz    = run%kz
          law%omega = run%omega
          law%taper_to_boundary = run%taper
        end select
      case default
        ! static, and every block but the receding one under nodeshift
        allocate( morph_motion_static_t :: domain%blk(b)%law )
        domain%blk(b)%law%name = 'static'
      end select
    enddo

    if ( trim(run%law) /= 'nodeshift' ) return

    ! the BC 306 convention of FUSS: faces 1/2 -> i, 3/4 -> j, 5/6 -> k; odd = the low-index face
    b = run%shift_block
    d = ( run%shift_face + 1 ) / 2
    deallocate( domain%blk(b)%law )
    allocate( domain%blk(b)%law, &
              source = morph_motion_nodeshift_t( name='nodeshift', dir=d, side=2-mod(run%shift_face,2) ) )
    select case ( d )
      case (1) ; t1 = 2 ; t2 = 3
      case (2) ; t1 = 1 ; t2 = 3
      case default ; t1 = 1 ; t2 = 2
    end select
    select type ( law => domain%blk(b)%law )
    type is ( morph_motion_nodeshift_t )
      allocate( law%disp( 0:domain%blk(b)%dim(t1), 0:domain%blk(b)%dim(t2) ) )
      do c = 0, domain%blk(b)%dim(t2)
        do a = 0, domain%blk(b)%dim(t1)
          law%disp(a,c)%c = run%shift_vel * run%dt
        enddo
      enddo
    end select

  end subroutine install_laws

end module Morph_Run_m
