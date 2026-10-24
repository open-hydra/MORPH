!>@brief The driver's input, one object per input.ini section (Register_Morph
!>       binds every key to a component below).
!>
!> The [MORPH-Motion] keys are those of FUSS's [FUSS-MeshMotion] section
!> (FUSS Register_Numerics.f90:78-130 at 2943354), with the same meaning and
!> defaults, so a FUSS case's motion runs unchanged here; the nodeshift keys
!> stand in for FUSS's recession, which hands the law its displacement.
module Morph_Config_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private

  integer,          parameter, public :: llen     = 256
  character(len=*), parameter, public :: codename = 'MORPH'

  !> [MORPH-IO]
  type, public :: io_t
    character(len=llen) :: mesh_in        !> mesh read at start
    character(len=llen) :: mesh_out       !> mesh written at the end
    character(len=llen) :: mesh_format    !> p3d | tecplot
    character(len=llen) :: gcl_log        !> per-step residual log
    integer             :: out_diter      !> also write the mesh every out_diter steps (0: never)
    integer             :: shell_diter    !> print a step line every shell_diter steps (0: never)
  end type io_t

  !> [MORPH-Time]
  type, public :: time_t
    real(R8)            :: t0             !> time of the mesh as read [s]
    real(R8)            :: dt             !> step [s]
    integer             :: steps          !> number of steps
  end type time_t

  !> [MORPH-Motion]
  type, public :: motion_t
    character(len=llen) :: law            !> static | prescribed | translation | nodeshift
    real(R8)            :: vel(3)         !> rigid translation velocity [m/s]
    real(R8)            :: amp(3)         !> prescribed: displacement amplitude per coordinate [m]
    real(R8)            :: kx(3)          !> prescribed: wavenumbers multiplying x [1/m]
    real(R8)            :: ky(3)
    real(R8)            :: kz(3)
    real(R8)            :: omega          !> prescribed: angular frequency [1/s]
    logical             :: taper          !> prescribed: taper the motion to zero on the block boundary
    real(R8)            :: gcl_tol        !> largest relative discrete-GCL residual accepted in a step
    integer             :: shift_block    !> nodeshift: the block whose face recedes
    integer             :: shift_face     !> nodeshift: that face, 1..6 (i-low, i-high, j-low, j-high, k-low, k-high)
    real(R8)            :: shift_vel(3)   !> nodeshift: velocity of every node of that face [m/s]
  end type motion_t

  type(io_t),     public, target :: obj_io
  type(time_t),   public, target :: obj_time
  type(motion_t), public, target :: obj_motion

end module Morph_Config_Types_m
