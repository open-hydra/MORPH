!>@brief N simulated ranks in one process, for the tests. Each rank is a
!> thread of an OpenMP parallel region of N threads that the test opens, and
!> the collectives meet at the team's barriers. The level-set blocks are split
!> among the ranks as a run on N processes would split them, and the values
!> pass by copy through mailboxes the ranks share, so one machine shows that
!> the bits do not depend on the number of ranks.
!>
!> USE. One world, made outside the parallel region, and one instance of
!> morph_exchange_sim_t per rank-thread, pointing to it:
!>     call world%init ( n, status )
!>     !$omp parallel num_threads(n)
!>       ... an instance ex local to the thread:
!>       ex%me = omp_get_thread_num()  ;  ex%world => world
!>       ... MORPH's calls with ex ...
!>     !$omp end parallel
!>
!> A COLLECTIVE is: thread t writes mailbox t (its data, and whether its own
!> arguments are consistent), barrier, every thread reads every mailbox,
!> barrier (nobody writes the next call's data before all have read). A
!> thread that returned before the first barrier would leave the others
!> waiting there for ever, so every refusal is taken by all the ranks:
!>  * the team must have n threads: OMP_DYNAMIC, OMP_THREAD_LIMIT or a call
!>    outside the region give fewer, and the missing ranks would leave their
!>    mailboxes stale. The team is the same for all its threads: they all
!>    refuse before writing anything;
!>  * thread t must be rank t, and its counts must fit its buffers: written
!>    in mailbox t, read by every rank after the barrier;
!>  * what two ranks must agree on (the counts one sends and the other
!>    expects, the counts of an allgather, the owners of the blocks) is
!>    compared after the barrier, by every rank on the same mailboxes.
!> Either every rank returns MORPH_ERR_EXCHANGE or none does. The world must
!> be the same for all the instances: one without it refuses alone.
!>
!> THREAD SAFETY. With these ranks, whatever MORPH runs between two exchanges
!> runs on n threads at once: the library writes no module variable and keeps
!> no SAVE'd local (a local declared with an initial value is SAVE'd).
!>
!> Without OpenMP there is one thread: a world of one rank works, a larger one
!> is refused by the team check.
module Morph_Exchange_Sim
  use iso_fortran_env, only: I4 => int32, R8 => real64
  !$ use omp_lib, only: omp_get_num_threads, omp_get_thread_num
  use Morph_Types_m
  use Morph_Exchange,  only: morph_exchange_t

  implicit none
  private

  !> What one rank posts in a collective.
  type, public :: morph_sim_box_t
    logical :: ok = .false.                  !> its own arguments are consistent
    integer(I4), allocatable :: cnt(:)       !> swap: counts it sends; allgather: the counts (0:n-1)
    integer(I4), allocatable :: rcnt(:)      !> swap: counts it expects (0:n-1)
    integer(I4), allocatable :: own(:)       !> ordered sum: the owner of each block
    real(R8),    allocatable :: v(:)         !> its values
  end type morph_sim_box_t

  !> The mailboxes shared by the n ranks: box(t) is written by rank t only.
  type, public :: morph_sim_world_t
    integer(I4) :: n = 0
    type(morph_sim_box_t), allocatable :: box(:)     !> (0:n-1)
  contains
    procedure :: init => world_init
  end type morph_sim_world_t

  type, extends(morph_exchange_t), public :: morph_exchange_sim_t
    integer(I4) :: me = 0                                  !> this rank
    type(morph_sim_world_t), pointer :: world => null()    !> shared by the n instances
  contains
    procedure :: rank        => sim_rank
    procedure :: nranks      => sim_nranks
    procedure :: swap        => sim_swap
    procedure :: allgather   => sim_allgather
    procedure :: ordered_sum => sim_ordered_sum
  end type morph_exchange_sim_t

contains

  !> A world of n ranks, made outside the parallel region.
  subroutine world_init ( self, n, status )
    class(morph_sim_world_t), intent(inout) :: self
    integer(I4),              intent(in)    :: n
    type(morph_status_t),     intent(out)   :: status

    call status%clear()
    if ( n < 1 ) then
      call status%fail ( MORPH_ERR_EXCHANGE, 'simulated ranks: a world needs at least one rank' )
      return
    endif
    if ( allocated(self%box) ) deallocate( self%box )
    allocate( self%box(0:n-1) )
    self%n = n
  end subroutine world_init


  integer(I4) function sim_rank ( self ) result ( r )
    class(morph_exchange_sim_t), intent(in) :: self
    r = self%me
  end function sim_rank


  !> The number of ranks of the world; 0 without one.
  integer(I4) function sim_nranks ( self ) result ( r )
    class(morph_exchange_sim_t), intent(in) :: self
    r = 0
    if ( associated(self%world) ) r = self%world%n
  end function sim_nranks


  !> The checks that every thread of the team decides alike: there is a world
  !> and the team has one thread per rank. t = this thread's number.
  logical function team_ok ( self, what, t, status ) result ( ok )
    class(morph_exchange_sim_t), intent(in)    :: self
    character(len=*),            intent(in)    :: what
    integer(I4),                 intent(out)   :: t
    type(morph_status_t),        intent(inout) :: status
    integer(I4) :: nteam
    character(len=MORPH_MSGLEN) :: msg

    t = 0
    nteam = 1
    !$ t = omp_get_thread_num()
    !$ nteam = omp_get_num_threads()
    ok = associated(self%world)
    if ( .not. ok ) then
      call status%fail ( MORPH_ERR_EXCHANGE, 'simulated '//what//': the instance has no world' )
      return
    endif
    ok = nteam == self%world%n
    if ( .not. ok ) then
      write(msg,'(A,A,A,I0,A,I0,A)') 'simulated ', what, ': a team of ', nteam, ' threads for ', &
            self%world%n, ' ranks'
      call status%fail ( MORPH_ERR_EXCHANGE, msg )
    endif
  end function team_ok


  !> Every rank posted consistent arguments (read after the first barrier).
  logical function all_posted ( world ) result ( ok )
    type(morph_sim_world_t), intent(in) :: world
    integer(I4) :: p
    ok = .true.
    do p = 0, world%n - 1
      if ( .not. world%box(p)%ok ) ok = .false.
    enddo
  end function all_posted


  !> dst(lo:lo+size(src)-1) = src, allocated anew: the bounds never depend on
  !> an earlier call. (An explicit allocation: gfortran 13.3 does not
  !> reallocate on assignment a component reached through the world pointer
  !> of a polymorphic argument, and writes through a null address.)
  subroutine put ( dst, src, lo )
    integer(I4), allocatable, intent(inout) :: dst(:)
    integer(I4),              intent(in)    :: src(:)
    integer(I4),              intent(in)    :: lo
    if ( allocated(dst) ) deallocate( dst )
    allocate( dst(lo:lo+size(src)-1) )
    dst = src
  end subroutine put


  !> dst(1:size(src)) = src, allocated anew (see put).
  subroutine put_v ( dst, src )
    real(R8), allocatable, intent(inout) :: dst(:)
    real(R8),              intent(in)    :: src(:)
    if ( allocated(dst) ) deallocate( dst )
    allocate( dst(size(src)) )
    dst = src
  end subroutine put_v


  subroutine sim_swap ( self, scount, sbuf, rcount, rbuf, status )
    class(morph_exchange_sim_t), intent(inout) :: self
    integer(I4),                 intent(in)    :: scount(0:)
    real(R8),                    intent(in)    :: sbuf(:)
    integer(I4),                 intent(in)    :: rcount(0:)
    real(R8),                    intent(out)   :: rbuf(:)
    type(morph_status_t),        intent(out)   :: status
    integer(I4) :: t, n, p, q, ks, kr, m
    logical :: ok

    call status%clear()
    if ( .not. team_ok ( self, 'swap', t, status ) ) return
    n = self%world%n
    ok = self%me == t .and. size(scount) == n .and. size(rcount) == n
    if ( ok ) ok = minval(scount) >= 0 .and. sum(scount) <= size(sbuf) .and. sum(rcount) <= size(rbuf)
    if ( ok ) then
      call put ( self%world%box(t)%cnt,  scount, 0 )
      call put ( self%world%box(t)%rcnt, rcount, 0 )
      call put_v ( self%world%box(t)%v, sbuf(1:sum(scount)) )
    endif
    self%world%box(t)%ok = ok
    !$omp barrier
    ok = all_posted ( self%world )
    if ( ok ) then                                   ! p sends q what q expects from p
      do q = 0, n-1
        do p = 0, n-1
          if ( self%world%box(p)%cnt(q) /= self%world%box(q)%rcnt(p) ) ok = .false.
        enddo
      enddo
    endif
    if ( ok ) then
      kr = 0
      do p = 0, n-1                                  ! from each rank, in rank order
        ks = sum( self%world%box(p)%cnt(0:t-1) )     ! where p's values for t start
        m  = self%world%box(p)%cnt(t)
        rbuf(kr+1:kr+m) = self%world%box(p)%v(ks+1:ks+m)
        kr = kr + m
      enddo
    endif
    !$omp barrier
    if ( .not. ok ) call status%fail ( MORPH_ERR_EXCHANGE, &
         'simulated swap: the ranks, their counts or their buffers do not match' )
  end subroutine sim_swap


  subroutine sim_allgather ( self, mine, all, counts, status )
    class(morph_exchange_sim_t), intent(inout) :: self
    real(R8),                    intent(in)    :: mine(:)
    real(R8),                    intent(out)   :: all(:)
    integer(I4),                 intent(in)    :: counts(0:)
    type(morph_status_t),        intent(out)   :: status
    integer(I4) :: t, n, p, k, m
    logical :: ok

    call status%clear()
    if ( .not. team_ok ( self, 'allgather', t, status ) ) return
    n = self%world%n
    ok = self%me == t .and. size(counts) == n
    if ( ok ) ok = counts(t) == size(mine) .and. sum(counts) <= size(all)
    if ( ok ) then
      call put ( self%world%box(t)%cnt, counts, 0 )
      call put_v ( self%world%box(t)%v, mine )
    endif
    self%world%box(t)%ok = ok
    !$omp barrier
    ok = all_posted ( self%world )
    if ( ok ) then                                   ! the same counts on every rank
      do p = 1, n-1
        if ( any( self%world%box(p)%cnt /= self%world%box(0)%cnt ) ) ok = .false.
      enddo
    endif
    if ( ok ) then
      k = 0
      do p = 0, n-1                                  ! in rank order
        m = self%world%box(p)%cnt(p)
        all(k+1:k+m) = self%world%box(p)%v(1:m)
        k = k + m
      enddo
    endif
    !$omp barrier
    if ( .not. ok ) call status%fail ( MORPH_ERR_EXCHANGE, &
         'simulated allgather: the ranks, their counts or their buffers do not match' )
  end subroutine sim_allgather


  !> Rank p posts the partials of its blocks in block order; every rank sums,
  !> in block order, the k-th partial of owner(b): the additions of the serial
  !> exchange, one by one, so the same bits.
  subroutine sim_ordered_sum ( self, nblk, owner, partial, total, status )
    class(morph_exchange_sim_t), intent(inout) :: self
    integer(I4),                 intent(in)    :: nblk
    integer(I4),                 intent(in)    :: owner(:)
    real(R8),                    intent(in)    :: partial(:)
    real(R8),                    intent(out)   :: total
    type(morph_status_t),        intent(out)   :: status
    integer(I4), allocatable :: k(:)
    integer(I4) :: t, n, p, b
    logical :: ok

    call status%clear()
    total = 0.0_R8
    if ( .not. team_ok ( self, 'ordered sum', t, status ) ) return
    n = self%world%n
    ok = self%me == t .and. size(owner) == nblk
    if ( ok ) ok = all( owner >= 0 .and. owner < n ) .and. count( owner == t ) == size(partial)
    if ( ok ) then
      call put ( self%world%box(t)%own, owner, 1 )
      call put_v ( self%world%box(t)%v, partial )
    endif
    self%world%box(t)%ok = ok
    !$omp barrier
    ok = all_posted ( self%world )
    if ( ok ) then                                   ! the same owners on every rank
      do p = 1, n-1
        if ( size(self%world%box(p)%own) /= size(self%world%box(0)%own) ) then
          ok = .false.
        else if ( any( self%world%box(p)%own /= self%world%box(0)%own ) ) then
          ok = .false.
        endif
      enddo
    endif
    if ( ok ) then
      allocate( k(0:n-1) )
      k = 0
      do b = 1, nblk
        p = owner(b)
        k(p) = k(p) + 1
        total = total + self%world%box(p)%v(k(p))
      enddo
    endif
    !$omp barrier
    if ( .not. ok ) call status%fail ( MORPH_ERR_EXCHANGE, &
         'simulated ordered sum: the ranks, their owners or their partials do not match' )
  end subroutine sim_ordered_sum

end module Morph_Exchange_Sim
