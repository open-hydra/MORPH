!>@brief What the driver writes while MORPH runs: the per-step log of the
!>       discrete-GCL residual, the step lines on the screen and the snapshots.
!>
!> Morph_Run calls on_step with step 0 once the mesh is started, then after
!> every step, the refused one included; the library does no I/O, so all of it
!> is here. The settings are the driver's [MORPH-IO] and [MORPH-Time].
module Morph_Output_m
  use iso_fortran_env,      only: I4 => int32, output_unit
  use Morph_API,            only: morph_domain_t, morph_run_observer_t, morph_run_step_t
  use Morph_Config_Types_m, only: obj_io, obj_time, llen
  use Morph_IO_m,           only: Write_Mesh

  implicit none
  private

  public :: morph_output_t, Fail

  !> The driver's observer: the gcl-log is open from step 0 to close_log.
  type, extends(morph_run_observer_t) :: morph_output_t
    integer(I4) :: ulog   = 0
    logical     :: opened = .false.
  contains
    procedure :: on_step   => output_on_step
    procedure :: close_log => output_close_log
  end type morph_output_t

contains

  !> Step 0: open gcl-log and write its header. Step k: its log row, a step
  !> line every shell-diter steps and at the last, and a snapshot every
  !> out-diter steps before the last, unless MORPH refused the step.
  subroutine output_on_step ( self, domain, info )
    class(morph_output_t),  intent(inout) :: self
    type(morph_domain_t),   intent(in)    :: domain
    type(morph_run_step_t), intent(in)    :: info
    ! Local
    character(len=4096) :: message
    character(len=llen) :: snapshot
    integer(I4)         :: ios

    if ( info%step == 0 ) then
      open( newunit=self%ulog, file=trim(obj_io%gcl_log), status='replace', action='write', iostat=ios )
      if ( ios /= 0 ) then
        call Make_Parent ( obj_io%gcl_log )
        open( newunit=self%ulog, file=trim(obj_io%gcl_log), status='replace', action='write', iostat=ios )
        if ( ios /= 0 ) call Fail ( 1, '[ERROR] cannot open gcl-log '//trim(obj_io%gcl_log) )
      endif
      self%opened = .true.
      write(self%ulog,'(A)') '# MORPH per-step discrete-GCL residual: max over all cells of all blocks of'
      write(self%ulog,'(A)') '#   |(V_new - V_old) - sum_f dV_swept(f)|, absolute [m^3] and relative to the cell volume,'
      write(self%ulog,'(A)') '#   the block of the relative maximum, and the largest volume fraction a cell swept out'
      write(self%ulog,'(A)') '# step  t_end  gcl_abs  gcl_rel  block  outflow_ratio'
      return
    endif

    write(self%ulog,'(I8,1X,ES24.16E3,2(1X,ES12.5E3),1X,I5,1X,ES12.5E3)') info%step, domain%t, &
          info%gcl_abs, info%gcl_rel, info%gcl_block, info%outflow
    if ( obj_io%shell_diter > 0 ) then
      if ( mod(info%step, obj_io%shell_diter) == 0 .or. info%step == obj_time%steps ) &
        write(*,'(A,I8,A,ES14.7,A,ES10.3,A,ES10.3)') ' step ', info%step, '  t = ', domain%t, &
              '  GCL rel = ', info%gcl_rel, '  outflow = ', info%outflow
    endif

    ! The run stops after a refused step: no snapshot of it
    if ( .not. info%accepted ) return

    if ( obj_io%out_diter > 0 ) then
      if ( mod(info%step, obj_io%out_diter) == 0 .and. info%step < obj_time%steps ) then
        snapshot = Snapshot_Name ( obj_io%mesh_out, info%step )
        call Write_Mesh ( snapshot, obj_io%mesh_format, domain, message )
        if ( len_trim(message) > 0 ) call Fail ( 1, message )
      endif
    endif

  end subroutine output_on_step


  subroutine output_close_log ( self )
    class(morph_output_t), intent(inout) :: self
    if ( self%opened ) close( self%ulog )
    self%opened = .false.
  end subroutine output_close_log


  !> mesh-out with _<step> before the extension.
  function Snapshot_Name ( base, n ) result ( name )
    character(len=*), intent(in) :: base
    integer(I4),      intent(in) :: n
    character(len=llen)          :: name
    ! Local
    integer(I4)       :: dot, slash
    character(len=16) :: sn

    write(sn,'(I0)') n
    dot   = index( trim(base), '.', back=.true. )
    slash = index( trim(base), '/', back=.true. )
    if ( dot > slash + 1 ) then
      name = base(1:dot-1)//'_'//trim(sn)//base(dot:len_trim(base))
    else
      name = trim(base)//'_'//trim(sn)
    endif

  end function Snapshot_Name


  subroutine Make_Parent ( filename )
    character(len=*), intent(in) :: filename
    integer(I4) :: slash, stat
    slash = index( trim(filename), '/', back=.true. )
    if ( slash > 1 ) call execute_command_line( 'mkdir -p "'//filename(1:slash-1)//'"', exitstat=stat )
  end subroutine Make_Parent


  !> Print msg and stop with the driver's exit status code (1, 2 or 3).
  subroutine Fail ( code, msg )
    integer(I4),      intent(in) :: code
    character(len=*), intent(in) :: msg
    write(*,'(A)') trim(msg)
    ! before the runtime's ERROR STOP line on stderr, in a log that has both
    flush(output_unit)
    select case ( code )
      case (1) ; error stop 1
      case (2) ; error stop 2
      case default ; error stop 3
    end select
  end subroutine Fail

end module Morph_Output_m
