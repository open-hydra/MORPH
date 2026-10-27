!>@brief Surface recession: the nodes of one block face move by a given
!> displacement, the opposite face stays put, and the nodes between them follow
!> the deformation: along their grid line (the line law, plan 08 section 2.3)
!> or through the smoothing of Morph_Smooth, chosen by `interior`.
!>
!> The block face that moves is named by `dir` (the grid direction of the lines
!> that end on it: 1 = i, 2 = j, 3 = k) and `side` (1 = the low-index face,
!> 2 = the high-index face). `disp(a, b)` is the displacement of the surface
!> node of the line at transverse node indices (a, b): (j, k) for dir 1,
!> (i, k) for dir 2, (i, j) for dir 3. The caller (FUSS's recession, from the
!> ablating wall's recession rate) supplies it every step; the law displaces
!> INCREMENTALLY from node_old, like the translation law, so a restart needs
!> nothing beyond the mesh in the solution file.
!>
!> THE LINE LAW (interior = MORPH_INTERIOR_LINE, the default). A node at arc
!> length s from the pinned end of its line, on a line of length L, moves by
!> disp * s/L. Every cell of the line therefore shrinks by the same fraction
!> |disp|/L and a graded mesh keeps its grading; a cell cannot invert before
!> the whole line is consumed. s and L are summed from node_old along the
!> line, the same sum for both, so the surface node moves by exactly disp
!> (s/L = 1 bitwise) and the pinned node by exactly nothing. Each line moves on
!> its own: a displacement that varies across the face is imprinted, scaled by
!> s/L, on every node plane between the two faces.
!>
!> THE SMOOTHED INTERIOR (interior = MORPH_INTERIOR_SMOOTH). The two faces keep
!> their values (disp and 0, bit for bit as under the line law) and the four
!> edges of the block that join them keep the line law's closed form, evaluated
!> the same way in every block that holds the edge, so two blocks sharing one
!> agree bit for bit (the edges are what a neighbouring block sees in a 2-D
!> mesh; against the line path's own arithmetic, which the compiler may
!> contract into an FMA, one rounding can differ). Every other node takes the
!> displacement of the lineal-spring equilibrium Morph_Smooth solves, so a
!> non-uniform face displacement spreads into the block and dies away towards
!> the pinned face instead of being confined to each line. On a block whose
!> lines are straight and uniformly displaced the two interiors coincide to the
!> solver's tolerance. The solve starts from the previous step's field (`dsp`,
!> kept in the law; the line law on the first step), which makes the cost of a
!> steady recession a few sweeps per step; the field is therefore reproducible
!> only to the tolerance across a restart. The line-consumption check is the
!> same, with the same status and message, for both interiors.
!>
!> A line whose length would drop to |disp| or less returns MORPH_ERR_NEGVOL:
!> a status, never a stop (MORPH guard (c)); the caller decides. A zero
!> displacement everywhere leaves the nodes bit-identical, which the swept-
!> volume fill needs to give exact zeros.
module Morph_Motion_Nodeshift
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m
  use Morph_Motion_m
  use Morph_Smooth, only: morph_smooth_ctl_t, morph_smooth_info_t, Morph_Smooth_Lineal, &
                          MORPH_INTERIOR_LINE, MORPH_INTERIOR_SMOOTH

  implicit none
  private

  public :: morph_motion_nodeshift_t

  type, extends(morph_motion_t) :: morph_motion_nodeshift_t
    integer(I4) :: dir  = 0                          !> 1..3, the direction of the receding lines
    integer(I4) :: side = 0                          !> 1 = low-index face recedes, 2 = high-index face
    type(morph_vec3_t), allocatable :: disp(:,:)     !> surface-node displacement over the transverse node indices
    integer(I4) :: interior = MORPH_INTERIOR_LINE    !> how the nodes between the faces move (module header)
    type(morph_smooth_ctl_t)  :: smooth              !> MORPH_INTERIOR_SMOOTH: the solver's controls
    type(morph_smooth_info_t) :: smooth_info         !> MORPH_INTERIOR_SMOOTH: what the last step's solve did
    type(morph_vec3_t), allocatable :: dsp(:,:,:)    !> MORPH_INTERIOR_SMOOTH: the last step's displacement of every node
  contains
    procedure :: apply => nodeshift_apply
  end type morph_motion_nodeshift_t

contains

  subroutine nodeshift_apply ( self, node, node_old, dim, t, dt, status )
    class(morph_motion_nodeshift_t), intent(inout) :: self
    type(morph_vec3_t),              intent(inout) :: node(0:,0:,0:)
    type(morph_vec3_t),              intent(in)    :: node_old(0:,0:,0:)
    integer(I4),                     intent(in)    :: dim(3)
    real(R8),                        intent(in)    :: t, dt
    type(morph_status_t),            intent(out)   :: status
    ! Local
    integer(I4) :: t1, t2, a, b, m, n, bad_a, bad_b
    integer(I4) :: idx(3)
    real(R8)    :: d(3), Lline, s, seg(3)
    real(R8)    :: arc(0:maxval(dim))
    logical     :: neg, moving

    call status%clear()

    if ( self%dir < 1 .or. self%dir > 3 .or. self%side < 1 .or. self%side > 2 ) then
      call status%fail ( MORPH_ERR_SHAPE, 'nodeshift: dir must be 1..3 and side 1 or 2' )
      return
    endif
    select case ( self%dir )
      case (1) ; t1 = 2 ; t2 = 3
      case (2) ; t1 = 1 ; t2 = 3
      case default ; t1 = 1 ; t2 = 2
    end select
    if ( .not. allocated(self%disp) ) then
      call status%fail ( MORPH_ERR_SHAPE, 'nodeshift: no displacement' )
      return
    endif
    if ( lbound(self%disp,1) /= 0 .or. ubound(self%disp,1) /= dim(t1) .or. &
         lbound(self%disp,2) /= 0 .or. ubound(self%disp,2) /= dim(t2) ) then
      call status%fail ( MORPH_ERR_SHAPE, 'nodeshift: the displacement array does not match the receding face' )
      return
    endif
    if ( self%interior /= MORPH_INTERIOR_LINE .and. self%interior /= MORPH_INTERIOR_SMOOTH ) then
      call status%fail ( MORPH_ERR_SHAPE, 'nodeshift: interior must be MORPH_INTERIOR_LINE or MORPH_INTERIOR_SMOOTH' )
      return
    endif

    ! node already holds node_old (Morph_Update copies before calling apply).
    ! No displacement anywhere: leave every node bit-identical.
    moving = .false.
    do b = 0, dim(t2)
      do a = 0, dim(t1)
        if ( any( self%disp(a,b)%c /= 0.0_R8 ) ) moving = .true.
      enddo
    enddo
    if ( .not. moving ) return

    n = dim(self%dir)
    neg = .false. ; bad_a = 0 ; bad_b = 0

    if ( self%interior == MORPH_INTERIOR_SMOOTH ) then
      call apply_smooth ()
      return
    endif

    !$omp parallel do collapse(2) private(a,b,m,idx,d,Lline,s,seg,arc) reduction(.or.:neg)
    do b = 0, dim(t2)
    do a = 0, dim(t1)
      d = self%disp(a,b)%c
      if ( all( d == 0.0_R8 ) ) cycle
      ! arc(m) = arc length from the pinned end to the node m steps away
      ! from it, summed over node_old
      arc(0) = 0.0_R8
      do m = 1, n
        seg = point( node_old, line_idx( m, a, b ) ) - point( node_old, line_idx( m - 1, a, b ) )
        arc(m) = arc(m-1) + sqrt( dot_product( seg, seg ) )
      enddo
      Lline = arc(n)
      if ( Lline - sqrt( dot_product( d, d ) ) <= 0.0_R8 ) then
        neg = .true.
      else
        do m = 1, n
          s = arc(m)
          idx = line_idx( m, a, b )
          node(idx(1),idx(2),idx(3))%c = node_old(idx(1),idx(2),idx(3))%c + d * ( s / Lline )
        enddo
      endif
    enddo
    enddo
    !$omp end parallel do

    if ( neg ) call fail_consumed ()

  contains

    !> (i, j, k) of the node m steps from the pinned end on line (a, b)
    pure function line_idx ( mm, aa, bb ) result ( ijk )
      integer(I4), intent(in) :: mm, aa, bb
      integer(I4) :: ijk(3)
      integer(I4) :: l
      l = mm
      if ( self%side == 1 ) l = n - mm     ! the high-index face is pinned
      ijk(self%dir) = l
      ijk(t1) = aa
      ijk(t2) = bb
    end function line_idx

    pure function point ( nd, ijk ) result ( x )
      type(morph_vec3_t), intent(in) :: nd(0:,0:,0:)
      integer(I4),        intent(in) :: ijk(3)
      real(R8) :: x(3)
      x = nd(ijk(1),ijk(2),ijk(3))%c
    end function point

    !> The first offending line (serial, only on failure), then the status:
    !> its location is the line's two transverse node indices and its direction
    subroutine fail_consumed ()
      outer: do b = 0, dim(t2)
        do a = 0, dim(t1)
          d = self%disp(a,b)%c
          if ( all( d == 0.0_R8 ) ) cycle
          Lline = 0.0_R8
          do m = 1, n
            seg = point( node_old, line_idx( m, a, b ) ) - point( node_old, line_idx( m - 1, a, b ) )
            Lline = Lline + sqrt( dot_product( seg, seg ) )
          enddo
          if ( Lline - sqrt( dot_product( d, d ) ) <= 0.0_R8 ) then
            bad_a = a ; bad_b = b
            exit outer
          endif
        enddo
      enddo outer
      call status%fail ( MORPH_ERR_NEGVOL, 'nodeshift: the displacement consumes a whole grid line', &
                         bad_a, bad_b, self%dir )
    end subroutine fail_consumed

    !> MORPH_INTERIOR_SMOOTH. The displacement field `dsp` of every node is
    !> built as Morph_Smooth_Lineal expects it (its header): the line law on
    !> the receding face, the pinned face and the four edges, every step; the
    !> line law everywhere on a first call (the start of the solve) and the
    !> previous step's field otherwise. Then the solve, then the move. The nodes
    !> are written only once everything succeeded, so a failed step leaves
    !> `node` as `node_old`.
    subroutine apply_smooth ()
      integer(I4) :: i, j, k, aa, bb, mm, ijk(3)
      real(R8)    :: dd(3), sg(3), LL, sc
      real(R8)    :: ar(0:n)
      logical     :: cold, edge

      cold = .not. allocated( self%dsp )
      if ( .not. cold ) cold = any( shape( self%dsp ) /= dim + 1 )
      if ( cold ) then
        if ( allocated( self%dsp ) ) deallocate( self%dsp )
        allocate( self%dsp( 0:dim(1), 0:dim(2), 0:dim(3) ) )
      endif

      sc = 0.0_R8
      !$omp parallel do collapse(2) private(aa,bb,mm,ijk,dd,sg,LL,ar,edge) reduction(.or.:neg) reduction(max:sc)
      do bb = 0, dim(t2)
      do aa = 0, dim(t1)
        dd   = self%disp(aa,bb)%c
        sc   = max( sc, maxval( abs( dd ) ) )
        edge = ( aa == 0 .or. aa == dim(t1) ) .and. ( bb == 0 .or. bb == dim(t2) )
        ijk  = line_idx( 0, aa, bb )
        self%dsp(ijk(1),ijk(2),ijk(3))%c = 0.0_R8                     ! the pinned node
        if ( all( dd == 0.0_R8 ) ) then
          ! a line whose surface node does not move: the line law is zero
          do mm = 1, n
            if ( mm == n .or. edge .or. cold ) then
              ijk = line_idx( mm, aa, bb )
              self%dsp(ijk(1),ijk(2),ijk(3))%c = 0.0_R8
            endif
          enddo
        else
          ar(0) = 0.0_R8
          do mm = 1, n
            sg = point( node_old, line_idx( mm, aa, bb ) ) - point( node_old, line_idx( mm - 1, aa, bb ) )
            ar(mm) = ar(mm-1) + sqrt( dot_product( sg, sg ) )
          enddo
          LL = ar(n)
          if ( LL - sqrt( dot_product( dd, dd ) ) <= 0.0_R8 ) then
            neg = .true.
          else
            do mm = 1, n
              if ( mm == n .or. edge .or. cold ) then
                ijk = line_idx( mm, aa, bb )
                self%dsp(ijk(1),ijk(2),ijk(3))%c = dd * ( ar(mm) / LL )
              endif
            enddo
          endif
        endif
      enddo
      enddo
      !$omp end parallel do

      if ( neg ) then
        call fail_consumed ()
        return
      endif

      call Morph_Smooth_Lineal ( node_old, dim, self%dir, self%side, sc, self%dsp, &
                                 self%smooth, self%smooth_info, status )
      if ( .not. status%ok() ) return

      !$omp parallel do collapse(3) private(i,j,k)
      do k = 0, dim(3)
      do j = 0, dim(2)
      do i = 0, dim(1)
        node(i,j,k)%c = node_old(i,j,k)%c + self%dsp(i,j,k)%c
      enddo; enddo; enddo
      !$omp end parallel do
    end subroutine apply_smooth

  end subroutine nodeshift_apply

end module Morph_Motion_Nodeshift
