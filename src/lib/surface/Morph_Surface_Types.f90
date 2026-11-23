!>@brief Types of the surface layer: boundary kinds and connections of every
!> boundary cell face (from the MORPH boundary-condition file written by
!> ATLAS, or from a coupler), the planes, the wall as one surface across the
!> blocks, the fixed-node masks.
!>
!> Face-local axes (m, n) of a block face: faces 1, 2 -> (j, k); 3, 4 -> (i, k);
!> 5, 6 -> (i, j), the indices of the B<b>F<f> zones of the wall files.
!>
!> MORPH must never `use` a FUSS module -- see Morph_Types_m.
module Morph_Surface_Types
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m

  implicit none
  public

  !> Kind of a boundary cell face. No solver code enters MORPH: the codes
  !> 701-703 of the file are MORPH's own.
  integer(I4), parameter :: MORPH_BND_OTHER     = 0   !> 0 (null face): the volume law keeps it as given
  integer(I4), parameter :: MORPH_BND_RECEDE    = 1   !> 701: rate from the wall files or the coupler
  integer(I4), parameter :: MORPH_BND_IMPOSED   = 2   !> 702: imposed rate, 0 = inhibited
  integer(I4), parameter :: MORPH_BND_FIXEDSURF = 3   !> 703: fixed surface the contact nodes slide on
  integer(I4), parameter :: MORPH_BND_PLANE     = 4   !> 200, 300: plane fitted from its nodes
  integer(I4), parameter :: MORPH_BND_CONN      = 5   !> 101: the nodes are copies of the neighbour's
  integer(I4), parameter :: MORPH_BND_PERIODIC  = 6   !> 201: refused in a moving block (first version)
  integer(I4), parameter :: MORPH_BND_CHIMERA   = 7   !> 102: refused in a moving block

  !> One connected boundary cell face: one 101 record.
  type :: morph_conn_t
    integer(I4) :: b = 0, face = 0, ijk(3) = 0      !> our cell (1-based cell indices) and its face (1..6)
    integer(I4) :: b2 = 0, face2 = 0, ijk2(3) = 0   !> the neighbour's cell and face
    !> orientation map, column by column as MOSE reads it (Lib_Ghost.f90:536-583, MOSE 1e5cc8d):
    !> (c1 c2 c3 c4) = (C(1,1) C(2,1) C(1,2) C(2,2)); one step along our face-local axis 1 is
    !> (C(1,1), C(2,1)) on the partner face, along axis 2 it is (C(1,2), C(2,2))
    integer(I4) :: C(2,2) = 0
  end type morph_conn_t

  !> One block face: one value per cell face, face-local cell indices (m, n).
  type :: morph_bface_t
    integer(I4), allocatable :: kind(:,:)   !> MORPH_BND_*
    integer(I4), allocatable :: mat(:,:)    !> material (701, 702), 0 elsewhere
    real(R8),    allocatable :: rate(:,:)   !> imposed rate [m/s] (702), 0 elsewhere
  end type morph_bface_t

  !> One plane: a connected set of 200/300 cell faces.
  type :: morph_plane_t
    real(R8) :: p0(3) = 0.0_R8, n(3) = 0.0_R8   !> mean of its nodes and unit normal (least squares)
  end type morph_plane_t

  !> Fixed nodes of one block: .true. = zero displacement.
  type :: morph_mask_t
    logical, allocatable :: m(:,:,:)            !> (0:dim(1), 0:dim(2), 0:dim(3))
  end type morph_mask_t

  !> The wall as one surface across the blocks.
  type :: morph_surface_t
    integer(I4) :: nnode = 0, nface = 0
    !> global wall nodes and their copies (CSR): copies of node g are first(g) .. first(g+1)-1
    integer(I4), allocatable :: first(:)            !> (nnode+1)
    integer(I4), allocatable :: copy_b(:), copy_ijk(:,:)   !> (ncopy), (3, ncopy)
    integer(I4), allocatable :: owner(:)            !> (nnode): index of the copy whose position is the master
    !> wall faces (one per wall cell face) and their four global nodes, counter-clockwise seen from the fluid
    integer(I4), allocatable :: face_node(:,:)      !> (4, nface)
    integer(I4), allocatable :: face_cell(:,:)      !> (4, nface): b, face, m, n = the cell of the wall-file zone B<b>F<face>
    integer(I4), allocatable :: face_kind(:)        !> MORPH_BND_RECEDE or MORPH_BND_IMPOSED
    integer(I4), allocatable :: face_mat(:)         !> material at t = 0
    !> boundary edges of the wall: the two global nodes and what the wall touches there
    integer(I4) :: nedge = 0
    integer(I4), allocatable :: edge_node(:,:)      !> (2, nedge)
    integer(I4), allocatable :: edge_kind(:)        !> MORPH_BND_* of the non-wall cell face across the edge
    integer(I4), allocatable :: edge_ref(:)         !> plane id (PLANE) or fixed-surface component id (FIXEDSURF)
  end type morph_surface_t

end module Morph_Surface_Types
