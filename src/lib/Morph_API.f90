!>@brief The MORPH public façade. A solver calls this module and nothing else.
!>
!> Keeping the surface this narrow is what makes plan 09's extraction checklist
!> cheap: when MORPH becomes a standalone hydra solver, everything below stays
!> and only the caller changes.
!>
!> Typical use from a solver:
!>
!>   call Morph_Init   ( node, dim, geom, status )              ! once, at setup
!>   ...
!>   call Morph_Update ( node, node_old, dim, t, dt, law, geom, status )   ! each step
!>   call Morph_GCL_Residual ( geom, dim, r_abs, r_rel )        ! assert r_rel ~ 1e-15
!>
!> Or, to let MORPH hold the mesh (morph_domain_t) and run it, one call: the
!> law installed, every step taken and checked, each one shown to the caller's
!> observer (src/lib/run/Morph_Run.f90):
!>
!>   call Morph_Run ( domain, run, report, status, observer )
!>
!> MORPH must never `use` a FUSS module -- see Morph_Types_m.
module Morph_API
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m
  use Morph_Metrics,           only: Morph_Mesh_Type, Morph_Geom_Allocate
  use Morph_GCL,               only: Morph_GCL_Residual, Morph_GCL_Volume_Gap, &
                                     Morph_Swept_Outflow_Ratio
  use Morph_Quality,           only: Morph_Quality_Check
  use Morph_Motion_m,          only: morph_motion_t
  use Morph_Motion_Static,     only: morph_motion_static_t
  use Morph_Motion_Prescribed, only: morph_motion_prescribed_t
  use Morph_Motion_Translation,only: morph_motion_translation_t
  use Morph_Motion_Nodeshift,  only: morph_motion_nodeshift_t
  use Morph_Smooth,            only: morph_smooth_ctl_t, morph_smooth_info_t, &
                                     MORPH_INTERIOR_LINE, MORPH_INTERIOR_SMOOTH
  use Morph_Motion_Given,      only: morph_motion_given_t            ! positions decided elsewhere
  use Morph_Events                                                   ! event codes and lists
  use Morph_Block_m,           only: Morph_Init, Morph_Update        ! the block core
  use Morph_Domain_m,          only: morph_domain_t, morph_block_t   ! the multi-block domain
  use Morph_Surface_Types,     only: morph_conn_t, morph_bface_t, morph_plane_t, morph_mask_t, morph_surface_t, &
                                     MORPH_BND_OTHER, MORPH_BND_RECEDE, MORPH_BND_IMPOSED, MORPH_BND_FIXEDSURF, &
                                     MORPH_BND_PLANE, MORPH_BND_CONN, MORPH_BND_PERIODIC, MORPH_BND_CHIMERA
  use Morph_Surface,           only: morph_surface_stats_t, Morph_Surface_Build, Morph_Surface_Scatter
  use Morph_Fixed,             only: Morph_Fixed_Masks                ! the wall and the fixed nodes
  use Morph_Run_m,             only: morph_run_t, morph_run_step_t, morph_run_report_t, &
                                     morph_run_observer_t, Morph_Run   ! the whole run in one call

  implicit none
  private

  ! Types
  public :: morph_vec3_t, morph_tens3_t, morph_face_t, morph_dir_t
  public :: morph_geom_t, morph_status_t
  public :: morph_motion_t, morph_motion_static_t, morph_motion_prescribed_t
  public :: morph_motion_translation_t, morph_motion_nodeshift_t
  public :: morph_smooth_ctl_t, morph_smooth_info_t        ! the nodeshift law's smoothed interior
  public :: MORPH_INTERIOR_LINE, MORPH_INTERIOR_SMOOTH
  public :: morph_domain_t, morph_block_t
  public :: morph_run_t, morph_run_step_t, morph_run_report_t, morph_run_observer_t
  public :: morph_conn_t, morph_bface_t, morph_plane_t, morph_mask_t, morph_surface_t, morph_surface_stats_t
  public :: MORPH_BND_OTHER, MORPH_BND_RECEDE, MORPH_BND_IMPOSED, MORPH_BND_FIXEDSURF
  public :: MORPH_BND_PLANE, MORPH_BND_CONN, MORPH_BND_PERIODIC, MORPH_BND_CHIMERA
  public :: morph_motion_given_t, morph_nodefield_t
  public :: morph_event_t, morph_event_list_t, MORPH_EV_STOP_LIMIT
  public :: MORPH_EV_FACE_INERT, MORPH_EV_SPAN_BLOCKEDGE, MORPH_EV_PATCH_CONSUMED, MORPH_EV_FRONTS_MERGE
  public :: MORPH_EV_LINE_CONSUMED, MORPH_EV_BREAKTHROUGH, MORPH_EV_STEP_MIDWALL, MORPH_EV_QUALITY
  public :: MORPH_EV_MATERIAL, MORPH_EV_THIN, MORPH_EV_ANCHOR_BORN

  ! Status codes
  public :: MORPH_OK, MORPH_ERR_SINGULAR, MORPH_ERR_NEGVOL, MORPH_ERR_JACOBIAN
  public :: MORPH_ERR_SHAPE, MORPH_ERR_GCL
  public :: MORPH_ERR_INPUT, MORPH_ERR_STATE, MORPH_ERR_EXCHANGE, MORPH_ERR_SMOOTH

  ! Entry points
  public :: Morph_Init
  public :: Morph_Update
  public :: Morph_Mesh_Type
  public :: Morph_GCL_Residual
  public :: Morph_GCL_Volume_Gap
  public :: Morph_Swept_Outflow_Ratio
  public :: Morph_Quality_Check
  public :: Morph_Geom_Allocate
  public :: Morph_Run
  public :: Morph_Surface_Build, Morph_Surface_Scatter, Morph_Fixed_Masks


end module Morph_API
