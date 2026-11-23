!>@brief Test of the `given` law: the copy rule and the `pending` order.
!>
!> It pins:
!>   * the copy rule: a node not set keeps node_old bit for bit, a negative
!>     zero included, and a block with nothing set sweeps exactly zero volume;
!>   * a set node takes its target position;
!>   * `fill` hands the arrays over without copying them;
!>   * the order: an update without a new fill is refused, and so is a second
!>     fill before the update, which leaves the first positions in place;
!>   * the shapes: positions and mask must agree, and must match the block.
program test_morph_given
  use iso_fortran_env, only: I4 => int32, I8 => int64, R8 => real64
  use Morph_API
  implicit none
  integer(I4), parameter :: N = 2
  type(morph_vec3_t), allocatable :: node(:,:,:), node_old(:,:,:), pos(:,:,:)
  logical,            allocatable :: set(:,:,:)
  type(morph_geom_t)         :: geom
  type(morph_status_t)       :: st
  type(morph_motion_given_t) :: law
  integer(I4) :: dim(3), i, j, k
  logical     :: ok, same

  ok = .true.
  write(*,'(A)') '== Morph given-law test =='
  write(*,*)

  dim = [N, N, N]
  allocate( node(0:N,0:N,0:N), node_old(0:N,0:N,0:N) )
  do k = 0, N; do j = 0, N; do i = 0, N
    node(i,j,k)%c = [0.1_R8*i, 0.1_R8*j, 0.1_R8*k]
  enddo; enddo; enddo
  node(0,0,0)%c(1) = -0.0_R8                        ! a negative zero
  call Morph_Init ( node, dim, geom, st )
  call expect ( st%ok(), 'Morph_Init' )

  ! 1. nothing set
  call new_target ()
  call law%fill ( pos, set, st )
  call expect ( st%ok() .and. law%pending, 'fill accepted, positions pending' )
  call expect ( .not. allocated(pos) .and. .not. allocated(set), 'fill takes the arrays without copying them' )
  call Morph_Update ( node, node_old, dim, 0.0_R8, 1.0_R8, law, geom, st )
  call expect ( st%ok() .and. .not. law%pending, 'update applies the positions and lowers the flag' )
  same = .true.
  do k = 0, N; do j = 0, N; do i = 0, N
    if ( any( transfer(node(i,j,k)%c, 0_I8, 3) /= transfer(node_old(i,j,k)%c, 0_I8, 3) ) ) same = .false.
  enddo; enddo; enddo
  call expect ( same, 'nothing set: every node bit-identical to node_old' )
  call expect ( sign(1.0_R8, node(0,0,0)%c(1)) < 0.0_R8, 'a negative zero stays negative (copy, not old + 0)' )
  call expect ( all( geom%dV_swept == 0.0_R8 ), 'nothing set: swept volumes exactly zero' )

  ! 2. one node set
  call new_target ()
  set(2,2,2) = .true.
  pos(2,2,2)%c = [0.25_R8, 0.21_R8, 0.22_R8]
  call law%fill ( pos, set, st )
  call Morph_Update ( node, node_old, dim, 1.0_R8, 1.0_R8, law, geom, st )
  call expect ( st%ok() .and. all( node(2,2,2)%c == [0.25_R8, 0.21_R8, 0.22_R8] ), 'a set node takes its position' )

  ! 3. an update without a new fill
  call Morph_Update ( node, node_old, dim, 2.0_R8, 1.0_R8, law, geom, st )
  call expect ( st%code == MORPH_ERR_STATE, 'update without a new fill refused' )

  ! 4. two fills before an update
  call new_target ()
  set(1,1,1) = .true.
  pos(1,1,1)%c = [0.11_R8, 0.12_R8, 0.13_R8]
  call law%fill ( pos, set, st )
  call new_target ()
  call law%fill ( pos, set, st )
  call expect ( st%code == MORPH_ERR_STATE .and. allocated(pos) .and. allocated(set), &
                'second fill before the update refused, arrays left to the caller' )
  call Morph_Update ( node, node_old, dim, 3.0_R8, 1.0_R8, law, geom, st )
  call expect ( st%ok() .and. all( node(1,1,1)%c == [0.11_R8, 0.12_R8, 0.13_R8] ), &
                'the first fill is the one applied' )

  ! 5. positions and mask of different shape
  deallocate( set )
  allocate( set(0:N-1,0:N,0:N) )
  set = .false.
  call law%fill ( pos, set, st )
  call expect ( st%code == MORPH_ERR_SHAPE .and. .not. law%pending, 'positions and mask of different shape refused' )

  ! 6. target arrays that do not match the block
  deallocate( pos, set )
  allocate( pos(0:N+1,0:N,0:N), set(0:N+1,0:N,0:N) )
  do k = 0, N; do j = 0, N; do i = 0, N+1
    pos(i,j,k)%c = 0.0_R8
  enddo; enddo; enddo
  set = .false.
  call law%fill ( pos, set, st )
  call Morph_Update ( node, node_old, dim, 4.0_R8, 1.0_R8, law, geom, st )
  call expect ( st%code == MORPH_ERR_SHAPE, 'target arrays that do not match the block refused' )

  write(*,*)
  if ( ok ) then
    write(*,'(A)') '   RESULT: PASS'
  else
    write(*,'(A)') '   RESULT: FAIL'
    error stop 1
  endif

contains

  subroutine expect ( good, what )
    logical,          intent(in) :: good
    character(len=*), intent(in) :: what
    if ( good ) then
      write(*,'(A,A,A)') '   ', what, ': PASS'
    else
      write(*,'(A,A)') '   FAIL: ', what
      ok = .false.
    endif
  end subroutine expect

  !> Fresh target arrays for the block: positions = the current nodes, nothing set.
  subroutine new_target ()
    if ( allocated(pos) ) deallocate( pos )
    if ( allocated(set) ) deallocate( set )
    allocate( pos(0:N,0:N,0:N), set(0:N,0:N,0:N) )
    pos = node
    set = .false.
  end subroutine new_target

end program test_morph_given
