!>@brief The MORPH driver: read a multi-block mesh, move it under a motion law
!>       for a number of steps, write it (plan 09 section 1, item 3).
!>
!> Run from a case directory holding input.ini ([MORPH-IO], [MORPH-Time],
!> [MORPH-Motion]; docs: README.md). The computation is one library call,
!> Morph_Run: it installs the law on every block, moves the mesh step by step
!> and refuses a step whose discrete-GCL residual exceeds gcl-tolerance, as
!> FUSS's Update_Mesh asserts it. This program reads the input, makes the call
!> and writes: the per-step log, the step lines and the snapshots through its
!> observer (io/Morph_Output.f90), the mesh at the end.
!>
!> MORPH the library never stops; this program is its caller and does, with a
!> status that says why:
!>   1  input (input.ini, the mesh file, the motion parameters, writing the mesh)
!>   2  MORPH returned an error (inverted cell, singular metric, line consumed,
!>      smoothing not converged)
!>   3  the discrete GCL residual of a step exceeded gcl-tolerance
program MORPH_program
#if defined (_OPENMP)
  use omp_lib
#endif
  use iso_fortran_env,      only: I4 => int32
  use Morph_API
  use Morph_Config_Types_m, only: obj_io, obj_time, obj_motion
  use Morph_Read_Ini_m,     only: Read_Inifile
  use Morph_IO_m,           only: Read_Mesh, Write_Mesh
  use Morph_Output_m,       only: morph_output_t, Fail
  implicit none
  type(morph_domain_t)     :: domain
  type(morph_run_t)        :: run
  type(morph_run_report_t) :: report
  type(morph_status_t)     :: status
  type(morph_output_t)     :: output
  character(len=4096)      :: message
  integer(I4)              :: nthreads
  logical                  :: smooth

#if defined (_OPENMP)
  !$omp parallel
  !$omp master
  nthreads = omp_get_num_threads()
  !$omp end master
  !$omp end parallel
  write(*,'(A)')    ' MORPH: parallel execution'
  write(*,'(A,I4)') ' -  Number of threads --> ', nthreads
#else
  nthreads = 1
  write(*,'(A)')    ' MORPH: serial execution'
#endif

  ! ---- Input
  call Read_Inifile ( 'input.ini', message )
  if ( len_trim(message) > 0 ) call Fail ( 1, message )

  call Read_Mesh ( obj_io%mesh_in, obj_io%mesh_format, domain, message )
  if ( len_trim(message) > 0 ) call Fail ( 1, message )
  write(*,'(A,I0,A,A)') ' -  Mesh              --> ', domain%nblocks, ' block(s) from ', trim(obj_io%mesh_in)
  write(*,'(A,A)')     ' -  Motion law        --> ', trim(obj_motion%law)
  smooth = ( trim(obj_motion%law) == 'nodeshift' .and. trim(obj_motion%interior) == 'smooth' )
  if ( trim(obj_motion%law) == 'nodeshift' ) write(*,'(A,A)') ' -  Interior          --> ', trim(obj_motion%interior)
  if ( smooth ) write(*,'(A,ES9.2,A,I0,A,ES9.2,A)') ' -  Smoothing         --> tolerance ', obj_motion%smooth_tol, &
                      ', maxit ', obj_motion%smooth_maxit, ', omega ', obj_motion%smooth_omega, ' (0: 2/(1 + sin(pi/N)))'
  write(*,'(A,I0,A,ES12.5,A)') ' -  Steps             --> ', obj_time%steps, ' of ', obj_time%dt, ' s'

  run%law         = obj_motion%law(1:len(run%law))   ! the names Register_Morph allows all fit
  run%t0          = obj_time%t0
  run%dt          = obj_time%dt
  run%steps       = obj_time%steps
  run%gcl_tol     = obj_motion%gcl_tol
  run%vel         = obj_motion%vel
  run%amp         = obj_motion%amp
  run%kx          = obj_motion%kx
  run%ky          = obj_motion%ky
  run%kz          = obj_motion%kz
  run%omega       = obj_motion%omega
  run%taper       = obj_motion%taper
  run%shift_block = obj_motion%shift_block
  run%shift_face  = obj_motion%shift_face
  run%shift_vel   = obj_motion%shift_vel
  run%shift_bump_vel   = obj_motion%shift_bump_vel
  run%shift_bump_width = obj_motion%shift_bump_width
  run%interior         = MORPH_INTERIOR_LINE              ! the registry allows line and smooth only
  if ( trim(obj_motion%interior) == 'smooth' ) run%interior = MORPH_INTERIOR_SMOOTH
  run%smooth_tol       = obj_motion%smooth_tol
  run%smooth_maxit     = obj_motion%smooth_maxit
  run%smooth_omega     = obj_motion%smooth_omega

  ! The observer reports the smoothing of the receding block under each step line
  output%smooth      = smooth
  output%shift_block = obj_motion%shift_block

  ! ---- The run: one call
  call Morph_Run ( domain, run, report, status, output )
  call output%close_log()

  if ( .not. status%ok() ) then
    select case ( status%code )
    case ( MORPH_ERR_INPUT )
      call Fail ( 1, '[ERROR] [MORPH-Motion] '//trim(status%message) )
    case ( MORPH_ERR_GCL )
      write(message,'(A,I0,A,ES10.3,A,I0,A,ES10.3,A)') '[ERROR] step ', report%last%step, &
            ': discrete GCL residual ', report%last%gcl_rel, ' (relative, block ', report%last%gcl_block, &
            ') exceeds gcl-tolerance = ', run%gcl_tol, ': the swept volumes and the volume change disagree'
      call Fail ( 3, message )
    case default
      if ( report%step == 0 ) then
        call Fail_Status ( 'setup', status )
      else
        call Fail_Status ( 'step', status )
      endif
    end select
  endif

  ! ---- Output
  call Write_Mesh ( obj_io%mesh_out, obj_io%mesh_format, domain, message )
  if ( len_trim(message) > 0 ) call Fail ( 1, message )

  write(*,'(A,I0,A,ES14.7,A)') ' MORPH: ', domain%step, ' steps, t = ', domain%t, ' s'
  write(*,'(A,ES10.3,A,ES10.3,A)') ' MORPH: largest discrete GCL residual ', report%worst_rel, ' relative, ', &
                                   report%worst_abs, ' m^3 absolute'
  write(*,'(A,A)') ' MORPH: mesh written to ', trim(obj_io%mesh_out)
  write(*,'(A)')   ' MORPH: done'

  call domain%free()

contains

  !> A MORPH status, in its block status%b. The status's (i,j,k) is a cell,
  !> except for the nodeshift law's consumed line (the line's two transverse
  !> node indices and its direction) and the smoothing's MORPH_ERR_SMOOTH (the
  !> first node of the set that did not converge, or the node whose neighbour
  !> coincides with it); the message says which.
  subroutine Fail_Status ( where, st )
    character(len=*),     intent(in) :: where
    type(morph_status_t), intent(in) :: st
    character(len=1024) :: msg
    write(msg,'(A,I0,A,I0,A,3(I0,A),A)') '[ERROR] MORPH '//where//' failed in block ', st%b, &
          ', status ', st%code, ' at (', st%i, ',', st%j, ',', st%k, '): ', trim(st%message)
    call Fail ( 2, msg )
  end subroutine Fail_Status

end program MORPH_program
