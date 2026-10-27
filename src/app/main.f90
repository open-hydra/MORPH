!>@brief The MORPH driver: read a multi-block mesh, move it under a motion law
!>       for a number of steps, write it (plan 09 section 1, item 3).
!>
!> Run from a case directory holding input.ini ([MORPH-IO], [MORPH-Time],
!> [MORPH-Motion]; docs: README.md). Each step moves every block and rebuilds its
!> interior geometry (Morph_Domain_m); the discrete-GCL residual of the step is
!> logged and asserted against gcl-tolerance, as FUSS's Update_Mesh asserts it.
!>
!> MORPH the library never stops; this program is its caller and does, with a
!> status that says why:
!>   1  input (input.ini, the mesh file, writing the mesh)
!>   2  MORPH returned an error (inverted cell, singular metric, line consumed)
!>   3  the discrete GCL residual of a step exceeded gcl-tolerance
program MORPH_program
#if defined (_OPENMP)
  use omp_lib
#endif
  use iso_fortran_env,      only: I4 => int32, R8 => real64, output_unit
  use Morph_API
  use Morph_Domain_m,       only: morph_domain_t
  use Morph_Config_Types_m, only: obj_io, obj_time, obj_motion, llen
  use Morph_Read_Ini_m,     only: Read_Inifile
  use Morph_IO_m,           only: Read_Mesh, Write_Mesh
  implicit none
  type(morph_domain_t) :: domain
  type(morph_status_t) :: status
  character(len=4096)  :: message
  character(len=llen)  :: snapshot
  integer(I4)          :: nthreads, step, ib, ib_abs, ib_rel, ob, oi, oj, oc, ios, ulog
  real(R8)             :: r_abs, r_rel, ratio, worst_rel, worst_abs
  logical              :: smooth

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

  call Set_Laws ( message )
  if ( len_trim(message) > 0 ) call Fail ( 1, message )
  write(*,'(A,A)')     ' -  Motion law        --> ', trim(obj_motion%law)
  smooth = ( trim(obj_motion%law) == 'nodeshift' .and. trim(obj_motion%interior) == 'smooth' )
  if ( trim(obj_motion%law) == 'nodeshift' ) write(*,'(A,A)') ' -  Interior          --> ', trim(obj_motion%interior)
  if ( smooth ) write(*,'(A,ES9.2,A,I0,A,ES9.2,A)') ' -  Smoothing         --> tolerance ', obj_motion%smooth_tol, &
                      ', maxit ', obj_motion%smooth_maxit, ', omega ', obj_motion%smooth_omega, ' (0: 2/(1 + sin(pi/N)))'
  write(*,'(A,I0,A,ES12.5,A)') ' -  Steps             --> ', obj_time%steps, ' of ', obj_time%dt, ' s'

  call domain%init ( obj_time%t0, status, ib )
  if ( .not. status%ok() ) call Fail_Status ( 'setup', ib, status )

  open( newunit=ulog, file=trim(obj_io%gcl_log), status='replace', action='write', iostat=ios )
  if ( ios /= 0 ) then
    call Make_Parent ( obj_io%gcl_log )
    open( newunit=ulog, file=trim(obj_io%gcl_log), status='replace', action='write', iostat=ios )
    if ( ios /= 0 ) call Fail ( 1, '[ERROR] cannot open gcl-log '//trim(obj_io%gcl_log) )
  endif
  write(ulog,'(A)') '# MORPH per-step discrete-GCL residual: max over all cells of all blocks of'
  write(ulog,'(A)') '#   |(V_new - V_old) - sum_f dV_swept(f)|, absolute [m^3] and relative to the cell volume,'
  write(ulog,'(A)') '#   the block of the relative maximum, and the largest volume fraction a cell swept out'
  write(ulog,'(A)') '# step  t_end  gcl_abs  gcl_rel  block  outflow_ratio'

  ! ---- Steps
  worst_rel = 0.0_R8
  worst_abs = 0.0_R8
  do step = 1, obj_time%steps

    if ( trim(obj_motion%law) == 'nodeshift' ) call Set_Shift ( obj_time%dt )

    call domain%advance ( obj_time%dt, status, ib )
    if ( .not. status%ok() ) call Fail_Status ( 'step', ib, status )

    call domain%gcl_residual ( r_abs, r_rel, ib_abs, ib_rel )
    call domain%outflow_ratio ( ratio, ob, oi, oj, oc )
    worst_rel = max( worst_rel, r_rel )
    worst_abs = max( worst_abs, r_abs )

    write(ulog,'(I8,1X,ES24.16E3,2(1X,ES12.5E3),1X,I5,1X,ES12.5E3)') step, domain%t, r_abs, r_rel, ib_rel, ratio
    if ( obj_io%shell_diter > 0 ) then
      if ( mod(step, obj_io%shell_diter) == 0 .or. step == obj_time%steps ) then
        write(*,'(A,I8,A,ES14.7,A,ES10.3,A,ES10.3)') ' step ', step, '  t = ', domain%t, &
              '  GCL rel = ', r_rel, '  outflow = ', ratio
        if ( smooth ) call Report_Smooth ()
      endif
    endif

    if ( r_rel > obj_motion%gcl_tol ) then
      write(message,'(A,I0,A,ES10.3,A,I0,A,ES10.3,A)') '[ERROR] step ', step, ': discrete GCL residual ', r_rel, &
            ' (relative, block ', ib_rel, ') exceeds gcl-tolerance = ', obj_motion%gcl_tol, &
            ': the swept volumes and the volume change disagree'
      close(ulog)
      call Fail ( 3, message )
    endif

    if ( obj_io%out_diter > 0 ) then
      if ( mod(step, obj_io%out_diter) == 0 .and. step < obj_time%steps ) then
        snapshot = Snapshot_Name ( obj_io%mesh_out, step )
        call Write_Mesh ( snapshot, obj_io%mesh_format, domain, message )
        if ( len_trim(message) > 0 ) call Fail ( 1, message )
      endif
    endif

  enddo
  close(ulog)

  ! ---- Output
  call Write_Mesh ( obj_io%mesh_out, obj_io%mesh_format, domain, message )
  if ( len_trim(message) > 0 ) call Fail ( 1, message )

  write(*,'(A,I0,A,ES14.7,A)') ' MORPH: ', domain%step, ' steps, t = ', domain%t, ' s'
  write(*,'(A,ES10.3,A,ES10.3,A)') ' MORPH: largest discrete GCL residual ', worst_rel, ' relative, ', &
                                   worst_abs, ' m^3 absolute'
  write(*,'(A,A)') ' MORPH: mesh written to ', trim(obj_io%mesh_out)
  write(*,'(A)')   ' MORPH: done'

  call domain%free()

contains

  !> One law instance per block, from [MORPH-Motion]. The names are those of
  !> Register_Morph's allowed list; FUSS's adapter builds its laws the same way
  !> (one prescribed instance per block, Adapter_Morph.f90:154-176 at 2943354).
  subroutine Set_Laws ( msg )
    character(len=*), intent(out) :: msg
    ! Local
    integer(I4) :: b, f, d
    character(len=16) :: sb

    msg = ''

    select case ( trim(obj_motion%law) )

    case ( 'static' )
      ! domain%allocate_blocks put the static law on every block

    case ( 'translation' )
      do b = 1, domain%nblocks
        deallocate( domain%blk(b)%law )
        allocate( domain%blk(b)%law, source = morph_motion_translation_t( name='translation', vel=obj_motion%vel ) )
      enddo

    case ( 'prescribed' )
      do b = 1, domain%nblocks
        deallocate( domain%blk(b)%law )
        allocate( morph_motion_prescribed_t :: domain%blk(b)%law )
        select type ( law => domain%blk(b)%law )
        type is ( morph_motion_prescribed_t )
          law%name  = 'prescribed'
          law%amp   = obj_motion%amp
          law%kx    = obj_motion%kx
          law%ky    = obj_motion%ky
          law%kz    = obj_motion%kz
          law%omega = obj_motion%omega
          law%taper_to_boundary = obj_motion%taper
        end select
      enddo

    case ( 'nodeshift' )
      b = obj_motion%shift_block
      f = obj_motion%shift_face
      if ( b < 1 .or. b > domain%nblocks ) then
        write(sb,'(I0)') domain%nblocks
        msg = '[ERROR] [MORPH-Motion] nodeshift-block must name a block, 1..'//trim(sb)
        return
      endif
      if ( f < 1 .or. f > 6 ) then
        msg = '[ERROR] [MORPH-Motion] nodeshift-face must be 1..6 (i-low, i-high, j-low, j-high, k-low, k-high)'
        return
      endif
      if ( obj_motion%smooth_omega >= 2.0_R8 ) then
        msg = '[ERROR] [MORPH-Motion] smooth-omega must be below 2 (0 selects 2/(1 + sin(pi/N)))'
        return
      endif
      ! the BC 306 convention of FUSS (plan 08 section 2.3): faces 1/2 -> i,
      ! 3/4 -> j, 5/6 -> k; odd = the low-index face
      d = (f + 1) / 2
      deallocate( domain%blk(b)%law )
      allocate( domain%blk(b)%law, source = morph_motion_nodeshift_t( name='nodeshift', dir=d, side=2-mod(f,2) ) )
      ! the interior: the registry has already refused a name but line and smooth
      select type ( law => domain%blk(b)%law )
      type is ( morph_motion_nodeshift_t )
        law%interior = MORPH_INTERIOR_LINE
        if ( trim(obj_motion%interior) == 'smooth' ) law%interior = MORPH_INTERIOR_SMOOTH
        law%smooth%tol   = obj_motion%smooth_tol
        law%smooth%maxit = obj_motion%smooth_maxit
        law%smooth%omega = obj_motion%smooth_omega
      end select

    case default
      msg = '[ERROR] [MORPH-Motion] law '//trim(obj_motion%law)//' is not known'

    end select

  end subroutine Set_Laws


  !> This step's displacement of the receding face: nodeshift-vel * dt at every
  !> node of the face (FUSS forms it from the ablating wall's recession rate),
  !> plus, when nodeshift-bump-vel is set, nodeshift-bump-vel * g * dt with the
  !> bump g = (1 - r^2/width^2)^2 inside r < width and 0 outside, r the distance
  !> of the node from the centre of the face mapped to the unit square (a
  !> transverse direction of one cell, a 2-D mesh, does not count). verify.py
  !> repeats this arithmetic; without a bump the expression is the old one.
  subroutine Set_Shift ( dt )
    real(R8), intent(in) :: dt
    ! Local
    integer(I4) :: b, t1, t2, a, c, n1, n2
    real(R8)    :: xi, eta, r2, g, w
    logical     :: bump

    b    = obj_motion%shift_block
    bump = any( obj_motion%shift_bump_vel /= 0.0_R8 )
    w    = obj_motion%shift_bump_width
    select type ( law => domain%blk(b)%law )
    type is ( morph_motion_nodeshift_t )
      select case ( law%dir )
        case (1) ; t1 = 2 ; t2 = 3
        case (2) ; t1 = 1 ; t2 = 3
        case default ; t1 = 1 ; t2 = 2
      end select
      n1 = domain%blk(b)%dim(t1)
      n2 = domain%blk(b)%dim(t2)
      if ( .not. allocated(law%disp) ) allocate( law%disp( 0:n1, 0:n2 ) )
      do c = 0, n2
        do a = 0, n1
          if ( bump ) then
            r2 = 0.0_R8
            if ( n1 > 1 ) then
              xi = real(a,R8) / real(n1,R8) - 0.5_R8
              r2 = r2 + xi * xi
            endif
            if ( n2 > 1 ) then
              eta = real(c,R8) / real(n2,R8) - 0.5_R8
              r2  = r2 + eta * eta
            endif
            g = max( 0.0_R8, 1.0_R8 - r2 / ( w * w ) )
            g = g * g
            law%disp(a,c)%c = ( obj_motion%shift_vel + obj_motion%shift_bump_vel * g ) * dt
          else
            law%disp(a,c)%c = obj_motion%shift_vel * dt
          endif
        enddo
      enddo
    end select

  end subroutine Set_Shift


  !> What the step's smoothing did, from the law of the receding block.
  subroutine Report_Smooth ()
    select type ( law => domain%blk(obj_motion%shift_block)%law )
    type is ( morph_motion_nodeshift_t )
      write(*,'(A,I0,A,I0,A,ES9.2,A,F6.3,A,F6.3)') '          smoothing: sweeps ', law%smooth_info%sweeps_face, &
            ' (faces) + ', law%smooth_info%sweeps_interior, ' (interior), last update ', &
            max( law%smooth_info%update_face, law%smooth_info%update_interior ), ' of the displacement, omega ', &
            law%smooth_info%omega_face, ' / ', law%smooth_info%omega_interior
    end select
  end subroutine Report_Smooth


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


  !> A MORPH status, with the block MORPH does not know it is in. The status's
  !> (i,j,k) is a cell, except for the nodeshift law's consumed line (the
  !> line's two transverse node indices and its direction) and the smoothing's
  !> status 6 (the first node of the set that did not converge, or the node
  !> whose neighbour coincides with it); the message says which.
  subroutine Fail_Status ( where, b, st )
    character(len=*),     intent(in) :: where
    integer(I4),          intent(in) :: b
    type(morph_status_t), intent(in) :: st
    character(len=1024) :: msg
    write(msg,'(A,I0,A,I0,A,3(I0,A),A)') '[ERROR] MORPH '//where//' failed in block ', b, &
          ', status ', st%code, ' at (', st%i, ',', st%j, ',', st%k, '): ', trim(st%message)
    call Fail ( 2, msg )
  end subroutine Fail_Status


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

end program MORPH_program
