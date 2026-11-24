!>@brief The wall as one surface across the blocks, built from the per-cell
!> connections and boundary kinds held by the domain.
!>
!>  * Node classes: the four corners of every 101 record are joined to the
!>    partner's through the orientation map C, never through coordinates
!>    (in a production motor mesh, 295 shared wall-node copies are not
!>    bit-identical).
!>    Joined nodes must coincide within JOIN_TOL, the tolerance ATLAS joins
!>    faces with (docs/user-guide/bcb/connectivity.md, ATLAS 3696d55).
!>  * Orientation map, as MOSE reads it (Lib_Ghost.f90:536-583, MOSE 1e5cc8d):
!>    one step along our face-local axis 1 is (C(1,1), C(2,1)) on the partner
!>    face, along axis 2 it is (C(1,2), C(2,2)); a corner at offset
!>    (2a-1, 2b-1) from the face centre goes to C*(2a-1, 2b-1).
!>  * Planes: 200/300 cell faces sharing an edge are joined unless they belong to
!>    the same block and to different block faces. Structural rule, no angle
!>    threshold: in a production motor wedge the tip and web planes meet along
!>    the axis at 16.36 degrees inside one block; the two sides of an
!>    axisymmetric wedge are faces 5 and 6 of one block. The worst outcome is one
!>    plane split into two coplanar records. Then least squares (eigenvector of
!>    the smallest eigenvalue of the covariance) and planarity within plane_tol.
!>
!> MORPH must never `use` a FUSS module -- see Morph_Types_m.
module Morph_Surface
  use iso_fortran_env, only: I4 => int32, I8 => int64, R8 => real64
  use Morph_Types_m
  use Morph_Surface_Types
  use Morph_Domain_m,  only: morph_domain_t

  implicit none
  private

  real(R8), parameter :: JOIN_TOL = 1.0e-7_R8   !> [mesh units] ATLAS joins faces whose centres are this close
  real(R8), parameter :: PI = 3.14159265358979323846_R8

  !> Counts and diagnostics of a build, for the tests and the driver's log.
  type, public :: morph_surface_stats_t
    integer(I4) :: edges_by_faces(3) = 0                 !> wall edges with 1, 2, >= 3 wall faces
    integer(I4) :: nplane = 0, nfixed = 0, naxis = 0     !> plane components, fixed-surface components, collinear (axis) ones
    real(R8)    :: max_joined_deg = 0.0_R8               !> largest angle between normals of joined plane faces
    real(R8)    :: min_split_deg = 180.0_R8              !> smallest angle between plane faces sharing an edge, not joined
    real(R8)    :: max_plane_dev = 0.0_R8                !> largest distance of a plane node from its fitted plane
    integer(I4), allocatable :: wall_nodes_on_plane(:)   !> wall nodes on each plane component
    integer(I4), allocatable :: wall_nodes_on_fixed(:)   !> wall nodes on each fixed-surface component
    integer(I4), allocatable :: edges_on_plane(:), edges_on_fixed(:)   !> wall boundary edges per component
  end type morph_surface_stats_t

  public :: Morph_Surface_Build, Morph_Surface_Scatter, Morph_Node_Classes
  public :: Morph_Node_Offsets, Morph_Corner_IJK, Morph_Face_Cell_Local

contains

  ! ===========================================================================
  ! Indices
  ! ===========================================================================

  !> off(b) = number of nodes of the blocks before b; node (b,i,j,k) has id
  !> off(b) + 1 + i + (d1+1)*(j + (d2+1)*k).
  subroutine Morph_Node_Offsets ( domain, off, ntot )
    type(morph_domain_t),     intent(in)  :: domain
    integer(I4), allocatable, intent(out) :: off(:)
    integer(I4),              intent(out) :: ntot
    integer(I4) :: b
    allocate ( off(domain%nblocks) )
    ntot = 0
    do b = 1, domain%nblocks
      off(b) = ntot
      ntot = ntot + product( domain%blk(b)%dim + 1 )
    enddo
  end subroutine Morph_Node_Offsets


  pure integer(I4) function node_id ( off, dim, ijk ) result ( id )
    integer(I4), intent(in) :: off, dim(3), ijk(3)
    id = off + 1 + ijk(1) + (dim(1)+1) * ( ijk(2) + (dim(2)+1) * ijk(3) )
  end function node_id


  !> Face-local cell indices (m, n) of the cell ijk on face f.
  pure subroutine Morph_Face_Cell_Local ( f, ijk, m, n )
    integer(I4), intent(in)  :: f, ijk(3)
    integer(I4), intent(out) :: m, n
    select case ( f )
    case ( 1, 2 ) ; m = ijk(2) ; n = ijk(3)
    case ( 3, 4 ) ; m = ijk(1) ; n = ijk(3)
    case default  ; m = ijk(1) ; n = ijk(2)
    end select
  end subroutine Morph_Face_Cell_Local


  !> Node indices (0-based) of corner (a, b) in {0,1}^2 of the face-local cell
  !> (m, n) of face f of a block of dims dim.
  pure subroutine Morph_Corner_IJK ( f, dim, m, n, a, b, ijk )
    integer(I4), intent(in)  :: f, dim(3), m, n, a, b
    integer(I4), intent(out) :: ijk(3)
    integer(I4) :: p, q, nrm
    p = m - 1 + a ; q = n - 1 + b
    select case ( f )
    case ( 1 ) ; nrm = 0
    case ( 2 ) ; nrm = dim(1)
    case ( 3 ) ; nrm = 0
    case ( 4 ) ; nrm = dim(2)
    case ( 5 ) ; nrm = 0
    case default ; nrm = dim(3)
    end select
    select case ( f )
    case ( 1, 2 ) ; ijk = [ nrm, p, q ]
    case ( 3, 4 ) ; ijk = [ p, nrm, q ]
    case default  ; ijk = [ p, q, nrm ]
    end select
  end subroutine Morph_Corner_IJK


  ! ===========================================================================
  ! Union-find (path halving; the smaller index becomes the root: deterministic)
  ! ===========================================================================

  subroutine uf_find ( parent, x, r )
    integer(I4), intent(inout) :: parent(:)
    integer(I4), intent(in)    :: x
    integer(I4), intent(out)   :: r
    r = x
    do while ( parent(r) /= r )
      parent(r) = parent(parent(r))
      r = parent(r)
    enddo
  end subroutine uf_find


  subroutine uf_union ( parent, x, y )
    integer(I4), intent(inout) :: parent(:)
    integer(I4), intent(in)    :: x, y
    integer(I4) :: rx, ry
    call uf_find ( parent, x, rx )
    call uf_find ( parent, y, ry )
    if ( rx < ry ) then
      parent(ry) = rx
    elseif ( ry < rx ) then
      parent(rx) = ry
    endif
  end subroutine uf_union

  ! ===========================================================================
  ! Node classes from the connections
  ! ===========================================================================

  !> parent(id) after the joins of every 101 record; root(id) = the class of node id.
  subroutine Morph_Node_Classes ( domain, parent, status )
    type(morph_domain_t),     intent(in)  :: domain
    integer(I4), allocatable, intent(out) :: parent(:)
    type(morph_status_t),     intent(out) :: status
    integer(I4), allocatable :: off(:)
    integer(I4) :: ntot, ic, i, a, b, m, n, m2, n2, s1, s2, ijk(3), ijk2(3), d(2)
    real(R8)    :: x(3), y(3)

    call status%clear()
    call Morph_Node_Offsets ( domain, off, ntot )
    allocate ( parent(ntot) )
    do i = 1, ntot
      parent(i) = i
    enddo
    if ( .not. allocated(domain%conns) ) return
    do ic = 1, size(domain%conns)
      associate ( c => domain%conns(ic) )
        call Morph_Face_Cell_Local ( c%face,  c%ijk,  m,  n )
        call Morph_Face_Cell_Local ( c%face2, c%ijk2, m2, n2 )
        do b = 0, 1
          do a = 0, 1
            s1 = 2*a - 1 ; s2 = 2*b - 1
            d(1) = c%C(1,1)*s1 + c%C(1,2)*s2
            d(2) = c%C(2,1)*s1 + c%C(2,2)*s2
            if ( abs(d(1)) /= 1 .or. abs(d(2)) /= 1 ) then
              call status%fail ( MORPH_ERR_INPUT, 'connection: the orientation map is not a signed permutation', &
                                 c%ijk(1), c%ijk(2), c%ijk(3), c%b )
              return
            endif
            call Morph_Corner_IJK ( c%face,  domain%blk(c%b)%dim,  m,  n,  a, b, ijk )
            call Morph_Corner_IJK ( c%face2, domain%blk(c%b2)%dim, m2, n2, (d(1)+1)/2, (d(2)+1)/2, ijk2 )
            x = domain%blk(c%b)%node(ijk(1), ijk(2), ijk(3))%c
            y = domain%blk(c%b2)%node(ijk2(1), ijk2(2), ijk2(3))%c
            if ( norm2(x - y) > JOIN_TOL ) then
              call status%fail ( MORPH_ERR_INPUT, &
                'connection: joined nodes farther apart than 1e-7 (wrong orientation map or partner cell)', &
                ijk(1), ijk(2), ijk(3), c%b )
              return
            endif
            call uf_union ( parent, node_id(off(c%b), domain%blk(c%b)%dim, ijk), &
                                    node_id(off(c%b2), domain%blk(c%b2)%dim, ijk2) )
          enddo
        enddo
      end associate
    enddo
  end subroutine Morph_Node_Classes

  ! ===========================================================================
  ! Build
  ! ===========================================================================

  !> Fill domain%surf and domain%planes. MORPH_ERR_INPUT on: a 201
  !> or 102 face in a block that is not frozen; joined nodes farther apart than
  !> JOIN_TOL; a wall edge shared by three or more wall faces; a wall boundary
  !> edge with no non-wall neighbour face or with more than one; a plane
  !> component farther than plane_tol from its fitted plane.
  subroutine Morph_Surface_Build ( domain, plane_tol, status, stats )
    type(morph_domain_t),                  intent(inout) :: domain
    real(R8),                              intent(in)    :: plane_tol   !> [m]
    type(morph_status_t),                  intent(out)   :: status
    type(morph_surface_stats_t), optional, intent(out)   :: stats
    ! Local
    type(morph_surface_stats_t) :: st
    integer(I4), parameter :: CA(4) = [0, 1, 1, 0], CB(4) = [0, 0, 1, 1]
    integer(I4), allocatable :: parent(:), off(:), groot(:), ncp(:), fill(:)
    integer(I4), allocatable :: pf(:,:), pf_r(:,:), pf_comp(:), pf_par(:), comp_of(:)
    integer(I4), allocatable :: xf(:,:), xf_r(:,:), xf_comp(:), xf_par(:)
    real(R8),    allocatable :: pf_nrm(:,:)
    integer(I8), allocatable :: key(:), ckey(:)
    integer(I4), allocatable :: perm(:), ckind(:), cref(:), cperm(:)
    integer(I4) :: nb, ntot, b, f, m, n, i, j, k, r, g, nface, ncopy, nplf, nxf, e, ne, ncand, ncomp, nx
    integer(I4) :: corner(4), cijk(3), ijk(3), u, v, p, i0, i1, nm, kind
    real(R8)    :: xx(3,4), nrm(3), din(3)

    call status%clear()
    nb = domain%nblocks
    if ( .not. allocated(domain%bfaces) ) then
      call status%fail ( MORPH_ERR_INPUT, 'Morph_Surface_Build: the boundary kinds (bfaces) are not set' )
      return
    endif
    do b = 1, nb                                 ! periodic and chimera faces only in frozen blocks
      if ( is_frozen(b) ) cycle
      do f = 1, 6
        if ( any( domain%bfaces(b,f)%kind == MORPH_BND_PERIODIC .or. domain%bfaces(b,f)%kind == MORPH_BND_CHIMERA ) ) then
          call status%fail ( MORPH_ERR_INPUT, 'Morph_Surface_Build: a periodic or chimera face in a block that moves', b = b )
          return
        endif
      enddo
    enddo

    ! ---- node classes
    call Morph_Node_Classes ( domain, parent, status )
    if ( .not. status%ok() ) return
    call Morph_Node_Offsets ( domain, off, ntot )

    ! ---- wall classes: roots of the corners of the wall cell faces
    allocate ( groot(ntot) ) ; groot = 0
    nface = 0
    do b = 1, nb
      do f = 1, 6
        do n = 1, size(domain%bfaces(b,f)%kind, 2)
          do m = 1, size(domain%bfaces(b,f)%kind, 1)
            if ( .not. is_wall(domain%bfaces(b,f)%kind(m,n)) ) cycle
            nface = nface + 1
            do i = 1, 4
              call Morph_Corner_IJK ( f, domain%blk(b)%dim, m, n, CA(i), CB(i), cijk )
              call uf_find ( parent, node_id(off(b), domain%blk(b)%dim, cijk), r )
              groot(r) = -1
            enddo
          enddo
        enddo
      enddo
    enddo

    ! ---- global numbering at first meeting (blocks in order, then k, j, i with i fastest); all copies of a class
    call clear_surface ( domain%surf )
    allocate ( ncp(ntot) ) ; ncp = 0
    do b = 1, nb
      do k = 0, domain%blk(b)%dim(3) ; do j = 0, domain%blk(b)%dim(2) ; do i = 0, domain%blk(b)%dim(1)
        call uf_find ( parent, node_id(off(b), domain%blk(b)%dim, [i,j,k]), r )
        if ( groot(r) == -1 ) then
          domain%surf%nnode = domain%surf%nnode + 1
          groot(r) = domain%surf%nnode
        endif
        if ( groot(r) > 0 ) ncp(groot(r)) = ncp(groot(r)) + 1
      enddo ; enddo ; enddo
    enddo
    associate ( s => domain%surf )
    allocate ( s%first(s%nnode+1), s%owner(s%nnode) )
    s%first(1) = 1
    do g = 1, s%nnode
      s%first(g+1) = s%first(g) + ncp(g)
    enddo
    ncopy = s%first(s%nnode+1) - 1
    allocate ( s%copy_b(ncopy), s%copy_ijk(3,ncopy), fill(s%nnode) )
    fill = s%first(1:s%nnode)
    do b = 1, nb
      do k = 0, domain%blk(b)%dim(3) ; do j = 0, domain%blk(b)%dim(2) ; do i = 0, domain%blk(b)%dim(1)
        call uf_find ( parent, node_id(off(b), domain%blk(b)%dim, [i,j,k]), r )
        g = groot(r)
        if ( g <= 0 ) cycle
        s%copy_b(fill(g)) = b ; s%copy_ijk(:,fill(g)) = [i,j,k]
        fill(g) = fill(g) + 1
      enddo ; enddo ; enddo
    enddo
    s%owner = s%first(1:s%nnode)                 ! the first copy met is the master

    ! ---- wall faces: blocks, block faces, cells (m fastest); counter-clockwise seen from the fluid
    s%nface = nface
    allocate ( s%face_node(4,nface), s%face_cell(4,nface), s%face_kind(nface), s%face_mat(nface) )
    e = 0
    do b = 1, nb
      do f = 1, 6
        do n = 1, size(domain%bfaces(b,f)%kind, 2)
          do m = 1, size(domain%bfaces(b,f)%kind, 1)
            if ( .not. is_wall(domain%bfaces(b,f)%kind(m,n)) ) cycle
            e = e + 1
            do i = 1, 4
              call Morph_Corner_IJK ( f, domain%blk(b)%dim, m, n, CA(i), CB(i), cijk )
              call uf_find ( parent, node_id(off(b), domain%blk(b)%dim, cijk), r )
              corner(i) = groot(r)
              xx(:,i) = domain%blk(b)%node(cijk(1), cijk(2), cijk(3))%c
            enddo
            call Morph_Corner_IJK ( f, domain%blk(b)%dim, m, n, 0, 0, ijk )   ! one layer inside, from corner 1
            select case ( f )
            case ( 1 ) ; ijk(1) = 1
            case ( 2 ) ; ijk(1) = domain%blk(b)%dim(1) - 1
            case ( 3 ) ; ijk(2) = 1
            case ( 4 ) ; ijk(2) = domain%blk(b)%dim(2) - 1
            case ( 5 ) ; ijk(3) = 1
            case default ; ijk(3) = domain%blk(b)%dim(3) - 1
            end select
            din = domain%blk(b)%node(ijk(1), ijk(2), ijk(3))%c - xx(:,1)
            nrm = cross( xx(:,3) - xx(:,1), xx(:,4) - xx(:,2) )
            if ( dot_product(nrm, din) >= 0.0_R8 ) then
              s%face_node(:,e) = corner
            else
              s%face_node(:,e) = [ corner(1), corner(4), corner(3), corner(2) ]
            endif
            s%face_cell(:,e) = [ b, f, m, n ]
            s%face_kind(e) = domain%bfaces(b,f)%kind(m,n)
            s%face_mat(e)  = domain%bfaces(b,f)%mat(m,n)
          enddo
        enddo
      enddo
    enddo

    ! ---- wall edges and their multiplicity; the boundary edges are those with one wall face
    allocate ( key(4*nface) )
    ne = 0
    do e = 1, nface
      do i = 1, 4
        u = s%face_node(i,e) ; v = s%face_node(mod(i,4)+1,e)
        if ( u == v ) cycle                      ! degenerate edge (both ends one node): skipped
        ne = ne + 1
        key(ne) = pair_key( u, v, s%nnode )
      enddo
    enddo
    allocate ( perm(ne) )
    call sort_i8 ( ne, key, perm )
    s%nedge = 0
    i = 1
    do while ( i <= ne )
      call run_end ( ne, key, perm, i, i1 )
      nm = i1 - i + 1
      st%edges_by_faces(min(nm,3)) = st%edges_by_faces(min(nm,3)) + 1
      if ( nm == 1 ) s%nedge = s%nedge + 1
      i = i1 + 1
    enddo
    if ( st%edges_by_faces(3) > 0 ) then
      call status%fail ( MORPH_ERR_INPUT, 'Morph_Surface_Build: a wall edge is shared by three or more wall faces' )
      return
    endif
    allocate ( s%edge_node(2,s%nedge), s%edge_kind(s%nedge), s%edge_ref(s%nedge) )
    e = 0 ; i = 1
    do while ( i <= ne )
      call run_end ( ne, key, perm, i, i1 )
      if ( i1 == i ) then
        e = e + 1
        s%edge_node(1,e) = int( key(perm(i)) / int(s%nnode+1,I8), I4 )
        s%edge_node(2,e) = int( mod( key(perm(i)), int(s%nnode+1,I8) ), I4 )
      endif
      i = i1 + 1
    enddo

    ! ---- plane faces: (b, f, m, n), corner roots, unit normal
    nplf = count_kind( MORPH_BND_PLANE )
    allocate ( pf(4,nplf), pf_r(4,nplf), pf_nrm(3,nplf), pf_par(nplf), pf_comp(nplf) )
    call list_faces ( MORPH_BND_PLANE, pf, pf_r )
    do p = 1, nplf
      do i = 1, 4
        call Morph_Corner_IJK ( pf(2,p), domain%blk(pf(1,p))%dim, pf(3,p), pf(4,p), CA(i), CB(i), cijk )
        xx(:,i) = domain%blk(pf(1,p))%node(cijk(1), cijk(2), cijk(3))%c
      enddo
      pf_nrm(:,p) = cross( xx(:,3) - xx(:,1), xx(:,4) - xx(:,2) )
      if ( norm2(pf_nrm(:,p)) > 0.0_R8 ) pf_nrm(:,p) = pf_nrm(:,p) / norm2(pf_nrm(:,p))
      pf_par(p) = p
    enddo
    call join_by_edges ( nplf, pf, pf_r, pf_nrm, ntot, pf_par, .true., st )
    call number_components ( nplf, pf_par, pf_comp, ncomp )
    allocate ( comp_of(ncomp) ) ; comp_of = 0
    call fit_planes ( domain, nplf, pf, pf_r, pf_comp, ncomp, ntot, plane_tol, comp_of, st, status )
    if ( .not. status%ok() ) return

    ! ---- fixed-surface components: faces joined through any shared edge
    nxf = count_kind( MORPH_BND_FIXEDSURF )
    allocate ( xf(4,nxf), xf_r(4,nxf), xf_par(nxf), xf_comp(nxf) )
    call list_faces ( MORPH_BND_FIXEDSURF, xf, xf_r )
    do p = 1, nxf
      xf_par(p) = p
    enddo
    call join_by_edges ( nxf, xf, xf_r, pf_nrm, ntot, xf_par, .false., st )
    call number_components ( nxf, xf_par, xf_comp, nx )
    st%nfixed = nx

    ! ---- kind of every wall boundary edge: the one non-wall, non-connection face that has it
    ncand = 0
    do b = 1, nb
      do f = 1, 6
        ncand = ncand + 4 * count( .not. ( domain%bfaces(b,f)%kind == MORPH_BND_RECEDE .or. &
                                           domain%bfaces(b,f)%kind == MORPH_BND_IMPOSED .or. &
                                           domain%bfaces(b,f)%kind == MORPH_BND_CONN ) )
      enddo
    enddo
    allocate ( ckey(ncand), ckind(ncand), cref(ncand) )
    ncand = 0
    p = 0 ; i0 = 0
    do b = 1, nb
      do f = 1, 6
        do n = 1, size(domain%bfaces(b,f)%kind, 2)
          do m = 1, size(domain%bfaces(b,f)%kind, 1)
            kind = domain%bfaces(b,f)%kind(m,n)
            if ( is_wall(kind) .or. kind == MORPH_BND_CONN ) cycle
            if ( kind == MORPH_BND_PLANE ) p = p + 1
            if ( kind == MORPH_BND_FIXEDSURF ) i0 = i0 + 1
            do i = 1, 4
              call Morph_Corner_IJK ( f, domain%blk(b)%dim, m, n, CA(i), CB(i), cijk )
              call uf_find ( parent, node_id(off(b), domain%blk(b)%dim, cijk), r )
              corner(i) = groot(r)
            enddo
            do i = 1, 4
              u = corner(i) ; v = corner(mod(i,4)+1)
              if ( u <= 0 .or. v <= 0 .or. u == v ) cycle
              ncand = ncand + 1
              ckey(ncand) = pair_key( u, v, s%nnode )
              ckind(ncand) = kind ; cref(ncand) = 0
              if ( kind == MORPH_BND_PLANE ) then
                ckind(ncand) = MORPH_BND_OTHER
                if ( comp_of(pf_comp(p)) > 0 ) then
                  ckind(ncand) = MORPH_BND_PLANE ; cref(ncand) = comp_of(pf_comp(p))
                endif
              endif
              if ( kind == MORPH_BND_FIXEDSURF ) cref(ncand) = xf_comp(i0)
            enddo
          enddo
        enddo
      enddo
    enddo
    allocate ( cperm(ncand) )
    call sort_i8 ( ncand, ckey, cperm )
    allocate ( st%edges_on_plane(st%nplane), st%edges_on_fixed(nx) )
    st%edges_on_plane = 0 ; st%edges_on_fixed = 0
    do e = 1, s%nedge
      call find_range ( ncand, ckey, cperm, pair_key( s%edge_node(1,e), s%edge_node(2,e), s%nnode ), i0, i1 )
      if ( i1 - i0 + 1 /= 1 ) then
        if ( i1 < i0 ) then
          call status%fail ( MORPH_ERR_INPUT, 'Morph_Surface_Build: a wall boundary edge touches no other boundary' )
        else
          call status%fail ( MORPH_ERR_INPUT, 'Morph_Surface_Build: a wall boundary edge touches two other boundaries' )
        endif
        return
      endif
      s%edge_kind(e) = ckind(cperm(i0))
      s%edge_ref(e)  = cref(cperm(i0))
      if ( s%edge_kind(e) == MORPH_BND_PLANE ) st%edges_on_plane(s%edge_ref(e)) = st%edges_on_plane(s%edge_ref(e)) + 1
      if ( s%edge_kind(e) == MORPH_BND_FIXEDSURF ) st%edges_on_fixed(s%edge_ref(e)) = st%edges_on_fixed(s%edge_ref(e)) + 1
    enddo

    ! ---- wall nodes on each component (diagnostics)
    allocate ( st%wall_nodes_on_plane(st%nplane), st%wall_nodes_on_fixed(nx) )
    call wall_nodes_on ( nplf, pf_r, pf_comp, comp_of, st%nplane, groot, st%wall_nodes_on_plane )
    call wall_nodes_on ( nxf, xf_r, xf_comp, [( i, i = 1, nx )], nx, groot, st%wall_nodes_on_fixed )
    end associate
    if ( present(stats) ) stats = st

  contains

    logical function is_frozen ( bb )
      integer(I4), intent(in) :: bb
      is_frozen = .false.
      if ( allocated(domain%frozen) ) is_frozen = domain%frozen(bb)
    end function is_frozen

    pure logical function is_wall ( kk )
      integer(I4), intent(in) :: kk
      is_wall = ( kk == MORPH_BND_RECEDE .or. kk == MORPH_BND_IMPOSED )
    end function is_wall

    integer(I4) function count_kind ( kk ) result ( c )
      integer(I4), intent(in) :: kk
      integer(I4) :: bb, ff
      c = 0
      do bb = 1, nb
        do ff = 1, 6
          c = c + count( domain%bfaces(bb,ff)%kind == kk )
        enddo
      enddo
    end function count_kind

    !> The faces of one kind, in block / face / cell order: (b, f, m, n) and the roots of the four corners.
    subroutine list_faces ( kk, fc, rr )
      integer(I4), intent(in)  :: kk
      integer(I4), intent(out) :: fc(:,:), rr(:,:)
      integer(I4) :: bb, ff, mm, nn, t, c, q(3)
      t = 0
      do bb = 1, nb
        do ff = 1, 6
          do nn = 1, size(domain%bfaces(bb,ff)%kind, 2)
            do mm = 1, size(domain%bfaces(bb,ff)%kind, 1)
              if ( domain%bfaces(bb,ff)%kind(mm,nn) /= kk ) cycle
              t = t + 1
              fc(:,t) = [ bb, ff, mm, nn ]
              do c = 1, 4
                call Morph_Corner_IJK ( ff, domain%blk(bb)%dim, mm, nn, CA(c), CB(c), q )
                call uf_find ( parent, node_id(off(bb), domain%blk(bb)%dim, q), rr(c,t) )
              enddo
            enddo
          enddo
        enddo
      enddo
    end subroutine list_faces

  end subroutine Morph_Surface_Build


  pure integer(I8) function pair_key ( u, v, nn ) result ( kv )
    integer(I4), intent(in) :: u, v, nn
    kv = int(min(u,v),I8) * int(nn+1,I8) + int(max(u,v),I8)
  end function pair_key


  !> i1 = last sorted position whose key equals that at position i.
  subroutine run_end ( n, key, perm, i, i1 )
    integer(I4), intent(in)  :: n, perm(:), i
    integer(I8), intent(in)  :: key(:)
    integer(I4), intent(out) :: i1
    i1 = i
    do while ( i1 < n )
      if ( key(perm(i1+1)) /= key(perm(i)) ) exit
      i1 = i1 + 1
    enddo
  end subroutine run_end


  !> Join faces that share an edge (corner roots rr). With `structural`, two
  !> faces of one block on different block faces are not joined, and the angles
  !> between normals are recorded as diagnostics.
  subroutine join_by_edges ( nfc, fc, rr, nrm, nn, par, structural, st )
    integer(I4),                 intent(in)    :: nfc, fc(:,:), rr(:,:), nn
    real(R8),                    intent(in)    :: nrm(:,:)
    integer(I4),                 intent(inout) :: par(:)
    logical,                     intent(in)    :: structural
    type(morph_surface_stats_t), intent(inout) :: st
    integer(I8), allocatable :: kk(:)
    integer(I4), allocatable :: ow(:), pp(:)
    integer(I4) :: t, c, ia, ib, ja, jb, nk, aa, bb, fa, fb
    real(R8)    :: an
    allocate ( kk(4*nfc), ow(4*nfc) )
    nk = 0
    do t = 1, nfc
      do c = 1, 4
        aa = rr(c,t) ; bb = rr(mod(c,4)+1,t)
        if ( aa == bb ) cycle
        nk = nk + 1
        kk(nk) = pair_key( aa, bb, nn )
        ow(nk) = t
      enddo
    enddo
    allocate ( pp(nk) )
    call sort_i8 ( nk, kk, pp )
    ia = 1
    do while ( ia <= nk )
      call run_end ( nk, kk, pp, ia, ib )
      do ja = ia, ib
        do jb = ja + 1, ib
          fa = ow(pp(ja)) ; fb = ow(pp(jb))
          if ( structural ) then
            an = angle_deg( nrm(:,fa), nrm(:,fb) )
            if ( fc(1,fa) == fc(1,fb) .and. fc(2,fa) /= fc(2,fb) ) then
              st%min_split_deg = min( st%min_split_deg, an )
              cycle
            endif
            st%max_joined_deg = max( st%max_joined_deg, an )
          endif
          call uf_union ( par, fa, fb )
        enddo
      enddo
      ia = ib + 1
    enddo
  end subroutine join_by_edges


  subroutine number_components ( nfc, par, comp, nc )
    integer(I4), intent(in)    :: nfc
    integer(I4), intent(inout) :: par(:)
    integer(I4), intent(out)   :: comp(:), nc
    integer(I4), allocatable :: lab(:)
    integer(I4) :: t, rt
    allocate ( lab(max(nfc,1)) ) ; lab = 0
    nc = 0
    do t = 1, nfc
      call uf_find ( par, t, rt )
      if ( lab(rt) == 0 ) then
        nc = nc + 1 ; lab(rt) = nc
      endif
      comp(t) = lab(rt)
    enddo
  end subroutine number_components


  !> Least-squares plane of every component over its distinct nodes (faces
  !> visited component by component, so a node shared by two components, as on
  !> the wedge axis, counts once in each). A collinear component is an axis: no
  !> plane, comp_of = 0. Then the planarity check.
  subroutine fit_planes ( domain, nfc, fc, rr, comp, ncomp, nn, tol, comp_of, st, status )
    type(morph_domain_t),        intent(inout) :: domain
    integer(I4),                 intent(in)    :: nfc, fc(:,:), rr(:,:), comp(:), ncomp, nn
    real(R8),                    intent(in)    :: tol
    integer(I4),                 intent(inout) :: comp_of(:)
    type(morph_surface_stats_t), intent(inout) :: st
    type(morph_status_t),        intent(inout) :: status
    integer(I4), parameter :: CA(4) = [0, 1, 1, 0], CB(4) = [0, 0, 1, 1]
    integer(I4), allocatable :: stamp(:), order(:), start(:)
    real(R8),    allocatable :: sx(:,:), cnt(:), cov(:,:,:), dev(:)
    real(R8) :: x(3), y(3), w(3), vec(3,3)
    integer(I4) :: t, it, c, ic, pass, q(3), r
    allocate ( stamp(nn), order(nfc), start(ncomp+1), sx(3,ncomp), cnt(ncomp), cov(3,3,ncomp), dev(ncomp) )
    start = 0                                   ! counting sort of the faces by component
    do t = 1, nfc
      start(comp(t)+1) = start(comp(t)+1) + 1
    enddo
    start(1) = 1
    do c = 1, ncomp
      start(c+1) = start(c+1) + start(c)
    enddo
    stamp(1:ncomp) = start(1:ncomp)            ! reused as fill pointer
    do t = 1, nfc
      order(stamp(comp(t))) = t
      stamp(comp(t)) = stamp(comp(t)) + 1
    enddo
    sx = 0 ; cnt = 0 ; cov = 0 ; dev = 0
    do pass = 1, 3                              ! 1: centroid, 2: covariance, 3: distance from the plane
      stamp = 0
      do it = 1, nfc
        t = order(it) ; c = comp(t)
        do ic = 1, 4
          r = rr(ic,t)
          if ( stamp(r) == c ) cycle
          stamp(r) = c
          call Morph_Corner_IJK ( fc(2,t), domain%blk(fc(1,t))%dim, fc(3,t), fc(4,t), CA(ic), CB(ic), q )
          x = domain%blk(fc(1,t))%node(q(1), q(2), q(3))%c
          select case ( pass )
          case ( 1 )
            sx(:,c) = sx(:,c) + x ; cnt(c) = cnt(c) + 1.0_R8
          case ( 2 )
            y = x - sx(:,c)
            cov(:,:,c) = cov(:,:,c) + spread(y,2,3) * spread(y,1,3)
          case default
            if ( comp_of(c) > 0 ) dev(c) = max( dev(c), abs( dot_product( x - domain%planes(comp_of(c))%p0, &
                                                                           domain%planes(comp_of(c))%n ) ) )
          end select
        enddo
      enddo
      if ( pass == 1 ) then
        do c = 1, ncomp
          sx(:,c) = sx(:,c) / cnt(c)
        enddo
      elseif ( pass == 2 ) then
        st%nplane = 0 ; st%naxis = 0
        do c = 1, ncomp
          call eig3 ( cov(:,:,c) / cnt(c), w, vec )
          if ( w(2) < tol**2 ) then             ! collinear nodes: an axis, not a plane
            st%naxis = st%naxis + 1
          else
            st%nplane = st%nplane + 1
            comp_of(c) = st%nplane
          endif
        enddo
        if ( allocated(domain%planes) ) deallocate ( domain%planes )
        allocate ( domain%planes(st%nplane) )
        do c = 1, ncomp
          if ( comp_of(c) == 0 ) cycle
          call eig3 ( cov(:,:,c) / cnt(c), w, vec )
          domain%planes(comp_of(c))%p0 = sx(:,c)
          domain%planes(comp_of(c))%n  = vec(:,1)
        enddo
      endif
    enddo
    do c = 1, ncomp
      if ( comp_of(c) == 0 ) cycle
      st%max_plane_dev = max( st%max_plane_dev, dev(c) )
      if ( dev(c) > tol ) then
        t = order(start(c))
        call status%fail ( MORPH_ERR_INPUT, 'Morph_Surface_Build: a plane component is not planar within plane_tol', &
                           fc(3,t), fc(4,t), fc(2,t), fc(1,t) )
        return
      endif
    enddo
  end subroutine fit_planes


  !> Distinct wall nodes (groot > 0) on the faces of each component; lab maps a
  !> component to its reported index (0 = not reported).
  subroutine wall_nodes_on ( nfc, rr, comp, lab, nlab, groot, cnts )
    integer(I4), intent(in)  :: nfc, rr(:,:), comp(:), lab(:), nlab, groot(:)
    integer(I4), intent(out) :: cnts(:)
    integer(I8), allocatable :: kk(:)
    integer(I4), allocatable :: pp(:)
    integer(I4) :: t, c, nk, i, i1, g
    allocate ( kk(4*max(nfc,1)) )
    nk = 0
    do t = 1, nfc
      if ( lab(comp(t)) == 0 ) cycle
      do c = 1, 4
        g = groot(rr(c,t))
        if ( g <= 0 ) cycle
        nk = nk + 1
        kk(nk) = int(g,I8) * int(nlab+1,I8) + int(lab(comp(t)),I8)
      enddo
    enddo
    allocate ( pp(max(nk,1)) )
    call sort_i8 ( nk, kk, pp )
    cnts = 0
    i = 1
    do while ( i <= nk )
      call run_end ( nk, kk, pp, i, i1 )
      t = int( mod( kk(pp(i)), int(nlab+1,I8) ), I4 )
      cnts(t) = cnts(t) + 1
      i = i1 + 1
    enddo
  end subroutine wall_nodes_on

  subroutine clear_surface ( s )
    type(morph_surface_t), intent(inout) :: s
    if ( allocated(s%first) )     deallocate( s%first )
    if ( allocated(s%copy_b) )    deallocate( s%copy_b )
    if ( allocated(s%copy_ijk) )  deallocate( s%copy_ijk )
    if ( allocated(s%owner) )     deallocate( s%owner )
    if ( allocated(s%face_node) ) deallocate( s%face_node )
    if ( allocated(s%face_cell) ) deallocate( s%face_cell )
    if ( allocated(s%face_kind) ) deallocate( s%face_kind )
    if ( allocated(s%face_mat) )  deallocate( s%face_mat )
    if ( allocated(s%edge_node) ) deallocate( s%edge_node )
    if ( allocated(s%edge_kind) ) deallocate( s%edge_kind )
    if ( allocated(s%edge_ref) )  deallocate( s%edge_ref )
    s%nnode = 0 ; s%nface = 0 ; s%nedge = 0
  end subroutine clear_surface

  ! ===========================================================================
  ! Small numerics
  ! ===========================================================================

  pure function cross ( a, b ) result ( c )
    real(R8), intent(in) :: a(3), b(3)
    real(R8) :: c(3)
    c = [ a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1) ]
  end function cross


  pure real(R8) function angle_deg ( a, b ) result ( d )
    real(R8), intent(in) :: a(3), b(3)
    real(R8) :: c
    c = max( -1.0_R8, min( 1.0_R8, abs( dot_product(a, b) ) ) )
    d = acos(c) * 180.0_R8 / PI
  end function angle_deg


  !> Eigenvalues (ascending) and eigenvectors (columns) of a symmetric 3x3
  !> matrix, by cyclic Jacobi rotations.
  pure subroutine eig3 ( a_in, w, v )
    real(R8), intent(in)  :: a_in(3,3)
    real(R8), intent(out) :: w(3), v(3,3)
    real(R8) :: a(3,3), th, t, c, s, tau, tmp(3)
    integer(I4) :: sweep, p, q, r, i
    a = a_in
    v = 0.0_R8
    do i = 1, 3
      v(i,i) = 1.0_R8
    enddo
    do sweep = 1, 50
      if ( abs(a(1,2)) + abs(a(1,3)) + abs(a(2,3)) == 0.0_R8 ) exit
      do p = 1, 2
        do q = p + 1, 3
          if ( a(p,q) == 0.0_R8 ) cycle
          th = ( a(q,q) - a(p,p) ) / ( 2.0_R8 * a(p,q) )
          t = sign(1.0_R8, th) / ( abs(th) + sqrt(th*th + 1.0_R8) )
          c = 1.0_R8 / sqrt(t*t + 1.0_R8) ; s = t*c ; tau = s / (1.0_R8 + c)
          a(p,p) = a(p,p) - t*a(p,q)
          a(q,q) = a(q,q) + t*a(p,q)
          a(p,q) = 0.0_R8 ; a(q,p) = 0.0_R8
          do r = 1, 3
            if ( r == p .or. r == q ) cycle
            tmp(1) = a(r,p) ; tmp(2) = a(r,q)
            a(r,p) = tmp(1) - s*(tmp(2) + tau*tmp(1)) ; a(p,r) = a(r,p)
            a(r,q) = tmp(2) + s*(tmp(1) - tau*tmp(2)) ; a(q,r) = a(r,q)
          enddo
          do r = 1, 3
            tmp(1) = v(r,p) ; tmp(2) = v(r,q)
            v(r,p) = tmp(1) - s*(tmp(2) + tau*tmp(1))
            v(r,q) = tmp(2) + s*(tmp(1) - tau*tmp(2))
          enddo
        enddo
      enddo
    enddo
    do i = 1, 3
      w(i) = a(i,i)
    enddo
    ! sort ascending
    do p = 1, 2
      do q = p + 1, 3
        if ( w(q) < w(p) ) then
          th = w(p) ; w(p) = w(q) ; w(q) = th
          tmp = v(:,p) ; v(:,p) = v(:,q) ; v(:,q) = tmp
        endif
      enddo
    enddo
  end subroutine eig3


  !> Stable merge sort: perm orders key ascending.
  subroutine sort_i8 ( n, key, perm )
    integer(I4), intent(in)  :: n
    integer(I8), intent(in)  :: key(:)
    integer(I4), intent(out) :: perm(:)
    integer(I4), allocatable :: tmp(:)
    integer(I4) :: width, lo, mid, hi, i, j, k
    do i = 1, n
      perm(i) = i
    enddo
    if ( n < 2 ) return
    allocate ( tmp(n) )
    width = 1
    do while ( width < n )
      lo = 1
      do while ( lo <= n )
        mid = min( lo + width - 1, n ) ; hi = min( lo + 2*width - 1, n )
        i = lo ; j = mid + 1 ; k = lo
        do while ( i <= mid .and. j <= hi )
          if ( key(perm(j)) < key(perm(i)) ) then
            tmp(k) = perm(j) ; j = j + 1
          else
            tmp(k) = perm(i) ; i = i + 1
          endif
          k = k + 1
        enddo
        do while ( i <= mid )
          tmp(k) = perm(i) ; i = i + 1 ; k = k + 1
        enddo
        do while ( j <= hi )
          tmp(k) = perm(j) ; j = j + 1 ; k = k + 1
        enddo
        lo = lo + 2*width
      enddo
      perm(1:n) = tmp(1:n)
      width = 2*width
    enddo
  end subroutine sort_i8


  !> The range [i0, i1] of sorted positions whose key equals kv (i1 < i0 if none).
  subroutine find_range ( n, key, perm, kv, i0, i1 )
    integer(I4), intent(in)  :: n, perm(:)
    integer(I8), intent(in)  :: key(:), kv
    integer(I4), intent(out) :: i0, i1
    integer(I4) :: lo, hi, mid
    lo = 1 ; hi = n
    do while ( lo <= hi )                       ! first position with key >= kv
      mid = (lo + hi) / 2
      if ( key(perm(mid)) < kv ) then
        lo = mid + 1
      else
        hi = mid - 1
      endif
    enddo
    i0 = lo ; i1 = lo - 1
    do while ( i1 + 1 <= n )
      if ( key(perm(i1+1)) /= kv ) exit
      i1 = i1 + 1
    enddo
  end subroutine find_range

  ! ===========================================================================
  ! Scatter
  ! ===========================================================================

  !> Give every copy of a wall node the bits of its owner copy. Called after
  !> every update of the wall nodes.
  subroutine Morph_Surface_Scatter ( domain, status )
    type(morph_domain_t), intent(inout) :: domain
    type(morph_status_t), intent(out)   :: status
    ! Local
    integer(I4) :: g, c, o
    real(R8)    :: x(3)

    call status%clear()
    associate ( s => domain%surf )
      if ( s%nnode == 0 ) return
      if ( .not. allocated(s%first) .or. .not. allocated(s%owner) .or. &
           .not. allocated(s%copy_b) .or. .not. allocated(s%copy_ijk) ) then
        call status%fail ( MORPH_ERR_SHAPE, 'Morph_Surface_Scatter: the surface is not built' )
        return
      endif
      do g = 1, s%nnode
        o = s%owner(g)
        x = domain%blk(s%copy_b(o))%node(s%copy_ijk(1,o), s%copy_ijk(2,o), s%copy_ijk(3,o))%c
        do c = s%first(g), s%first(g+1) - 1
          if ( c == o ) cycle
          domain%blk(s%copy_b(c))%node(s%copy_ijk(1,c), s%copy_ijk(2,c), s%copy_ijk(3,c))%c = x
        enddo
      enddo
    end associate
  end subroutine Morph_Surface_Scatter

end module Morph_Surface
