!>@brief Every input.ini key of the MORPH driver, registered as FUSS registers
!>       its own (FUSS src/lib/config/Register_*.f90): section, name, target,
!>       default, description, allowed values, required.
!>
!> The registry assigns the default when a key is added, Load_Ini overwrites
!> what the file sets, and Validate_Registry checks the required keys and the
!> allowed lists. A key in a [MORPH-*] section that is not registered here is
!> refused (Check_Unknown_Options), so a misspelt key cannot silently fall back
!> to its default.
module Morph_Register_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Config_Types_m
  use Morph_Input_Registry

  implicit none
  private
  public :: Register_Morph

contains

  subroutine Register_Morph ()
    implicit none
    ! Local
    character(len=llen) :: section

    !! ------------------------------------------------------
    !! Input and output -------------------------------------
    !! ------------------------------------------------------
    section = codename//'-IO'

    call reg%add( trim(section), 'mesh-in', obj_io%mesh_in, '', &
                  'Mesh read at start: PLOT3D ASCII multiblock (p3d) or Tecplot ASCII .tec / .dat (tecplot)', '', .true. )
    call reg%add( trim(section), 'mesh-out', obj_io%mesh_out, 'OUTPUT/mesh.p3d', &
                  'Mesh written at the end of the run (its directory is created)', '', .false. )
    ! p3d is exact (17 significant digits); the Tecplot writer keeps 15.
    call reg%add( trim(section), 'mesh-format', obj_io%mesh_format, 'p3d', &
                  'Format of mesh-in and mesh-out', 'p3d, tecplot', .false. )
    call reg%add( trim(section), 'gcl-log', obj_io%gcl_log, 'OUTPUT/gcl.dat', &
                  'Per-step log of the discrete-GCL residual and the swept-outflow ratio', '', .false. )
    call reg%add( trim(section), 'out-diter', obj_io%out_diter, '0', &
                  'Also write the mesh every out-diter steps, as mesh-out with _<step> before the extension (0: only at the end)', &
                  '>= 0', .false. )
    call reg%add( trim(section), 'shell-diter', obj_io%shell_diter, '1', &
                  'Print a step line every shell-diter steps (0: never)', '>= 0', .false. )

    !! ------------------------------------------------------
    !! Time -------------------------------------------------
    !! ------------------------------------------------------
    section = codename//'-Time'

    call reg%add( trim(section), 't0', obj_time%t0, '0.0', &
                  'Time of the mesh as read [s] (the prescribed law is a function of time)', 'real', .false. )
    call reg%add( trim(section), 'dt', obj_time%dt, '0.0', 'Step [s]', '> 0', .true. )
    call reg%add( trim(section), 'steps', obj_time%steps, '0', 'Number of steps', '>= 0', .true. )

    !! ------------------------------------------------------
    !! Motion: the keys of FUSS's [FUSS-MeshMotion] ---------
    !! ------------------------------------------------------
    section = codename//'-Motion'

    ! The allowed list is load-bearing: the driver selects the law by name, and
    ! Validate_Registry refuses a name not listed here. Keep it in sync with
    ! Set_Laws in main.f90.
    call reg%add( trim(section), 'law', obj_motion%law, 'static', &
                  'Mesh motion law', 'static ,  prescribed ,  translation ,  nodeshift', .false. )

    call reg%add( trim(section), 'vel', obj_motion%vel, '0.0', &
                  'Rigid-translation velocity [m/s]', 'real', .false. )

    call reg%add( trim(section), 'amp', obj_motion%amp, '0.0', &
                  'Prescribed-motion displacement amplitude per coordinate [m]', '>= 0', .false. )
    call reg%add( trim(section), 'kx', obj_motion%kx, '0.0', &
                  'Prescribed-motion wavenumbers multiplying x [1/m]', 'real', .false. )
    call reg%add( trim(section), 'ky', obj_motion%ky, '0.0', &
                  'Prescribed-motion wavenumbers multiplying y [1/m]', 'real', .false. )
    call reg%add( trim(section), 'kz', obj_motion%kz, '0.0', &
                  'Prescribed-motion wavenumbers multiplying z [1/m]', 'real', .false. )
    call reg%add( trim(section), 'omega', obj_motion%omega, '0.0', &
                  'Prescribed-motion angular frequency [1/s]', 'real', .false. )
    call reg%add( trim(section), 'taper-to-boundary', obj_motion%taper, '.false.', &
                  'Taper prescribed motion to zero on the block boundary', 'logical', .false. )

    ! Asserted every step, as FUSS's Update_Mesh does: a residual above it means
    ! the swept volumes and the volume change disagree.
    call reg%add( trim(section), 'gcl-tolerance', obj_motion%gcl_tol, '1.0e-10', &
                  'Max relative discrete-GCL residual before the run is stopped', '> 0', .false. )

    ! nodeshift: in FUSS the displacement comes from the ablating wall; here it
    ! is a constant velocity of every node of one block face.
    call reg%add( trim(section), 'nodeshift-block', obj_motion%shift_block, '0', &
                  'nodeshift: the block whose face recedes (every other block is static)', '>= 0', .false. )
    call reg%add( trim(section), 'nodeshift-face', obj_motion%shift_face, '0', &
                  'nodeshift: the receding face, 1..6 = i-low, i-high, j-low, j-high, k-low, k-high', '>= 0', .false. )
    call reg%add( trim(section), 'nodeshift-vel', obj_motion%shift_vel, '0.0', &
                  'nodeshift: velocity of the nodes of the receding face [m/s]', 'real', .false. )
    ! A displacement that varies across the face, to exercise the interior:
    ! vel + bump-vel * (1 - r^2/width^2)^2 inside r < width, r the distance
    ! from the face's centre over its unit square (Set_Shift in main.f90).
    call reg%add( trim(section), 'nodeshift-bump-vel', obj_motion%shift_bump_vel, '0.0', &
                  'nodeshift: velocity added at the centre of the receding face, times a bump of radius nodeshift-bump-width [m/s]', &
                  'real', .false. )
    call reg%add( trim(section), 'nodeshift-bump-width', obj_motion%shift_bump_width, '0.25', &
                  'nodeshift: radius of the bump over the face mapped to the unit square', '> 0', .false. )

    ! The interior of the receding block: Morph_Motion_Nodeshift's two options.
    ! The default is the line law, which FUSS's regression references hold.
    call reg%add( trim(section), 'interior', obj_motion%interior, 'line', &
                  'nodeshift: how the nodes between the receding and the pinned face move: '// &
                  'line (each grid line compressed uniformly) or smooth (the lineal-spring equilibrium of the face displacement)', &
                  'line ,  smooth', .false. )
    call reg%add( trim(section), 'smooth-tolerance', obj_motion%smooth_tol, '1.0e-12', &
                  'interior = smooth: a set of nodes is converged when a sweep moves none of them by more than '// &
                  'this fraction of the largest face displacement', '> 0', .false. )
    call reg%add( trim(section), 'smooth-maxit', obj_motion%smooth_maxit, '100000', &
                  'interior = smooth: sweeps allowed per set of nodes before MORPH returns status MORPH_ERR_SMOOTH', &
                  '>= 1', .false. )
    call reg%add( trim(section), 'smooth-omega', obj_motion%smooth_omega, '0.0', &
                  'interior = smooth: the over-relaxation factor, below 2; 0 = 2/(1 + sin(pi/N)) for the extent N of the set', &
                  '>= 0', .false. )

  end subroutine Register_Morph

end module Morph_Register_m
