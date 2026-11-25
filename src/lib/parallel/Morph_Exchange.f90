!>@brief How MORPH talks to the other ranks without knowing MPI: an abstract
!> exchange that the library calls and the caller implements, with its own
!> MPI or with the serial one below (the simulated ranks of the tests are in
!> Morph_Exchange_Sim). The library itself stays free of MPI.
!>
!> THE CONTRACT.
!>  * swap: the exchange of packed buffers of MPI_Alltoallv: scount(q) values
!>    go to rank q, taken from sbuf in rank order; rcount(q) values come from
!>    rank q, stored in rbuf in rank order.
!>  * allgather: the values of every rank, concatenated in rank order, on
!>    every rank; counts(q) is the number rank q gives.
!>  * ordered_sum: one partial value per block of the level-set grid, owner(b)
!>    the rank of block b, each rank giving the partials of its own blocks in
!>    block order; the total is summed in block order, so it has the same bits
!>    on every rank and for any number of ranks.
!>  * MORPH calls these outside its OpenMP regions only (one thread per rank
!>    talks to the others, as MPI_THREAD_FUNNELED requires) and fills the
!>    buffers in an order it fixes (block, then node), so the result does not
!>    depend on the implementation.
!>  * Every rank makes the same calls in the same order: a rank that skips one
!>    leaves the others waiting. MORPH decides anything that changes its calls
!>    from data that every rank shares.
!>
!> Never stops, no I/O: a call that does not match returns MORPH_ERR_EXCHANGE.
module Morph_Exchange
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Morph_Types_m

  implicit none
  private

  type, abstract, public :: morph_exchange_t
  contains
    procedure(rank_i), deferred :: rank          !> this rank, 0 .. nranks-1
    procedure(rank_i), deferred :: nranks
    procedure(swap_i), deferred :: swap
    procedure(allg_i), deferred :: allgather
    procedure(osum_i), deferred :: ordered_sum
  end type morph_exchange_t

  abstract interface
    integer(I4) function rank_i ( self )
      import :: morph_exchange_t, I4
      class(morph_exchange_t), intent(in) :: self
    end function rank_i
    subroutine swap_i ( self, scount, sbuf, rcount, rbuf, status )
      import :: morph_exchange_t, I4, R8, morph_status_t
      class(morph_exchange_t), intent(inout) :: self
      integer(I4),             intent(in)    :: scount(0:)   !> values for each rank, in rank order
      real(R8),                intent(in)    :: sbuf(:)
      integer(I4),             intent(in)    :: rcount(0:)   !> values from each rank, in rank order
      real(R8),                intent(out)   :: rbuf(:)
      type(morph_status_t),    intent(out)   :: status
    end subroutine swap_i
    subroutine allg_i ( self, mine, all, counts, status )
      import :: morph_exchange_t, I4, R8, morph_status_t
      class(morph_exchange_t), intent(inout) :: self
      real(R8),                intent(in)    :: mine(:)
      real(R8),                intent(out)   :: all(:)
      integer(I4),             intent(in)    :: counts(0:)   !> values of each rank, in rank order
      type(morph_status_t),    intent(out)   :: status
    end subroutine allg_i
    subroutine osum_i ( self, nblk, owner, partial, total, status )
      import :: morph_exchange_t, I4, R8, morph_status_t
      class(morph_exchange_t), intent(inout) :: self
      integer(I4),             intent(in)    :: nblk
      integer(I4),             intent(in)    :: owner(:)     !> rank of each block
      real(R8),                intent(in)    :: partial(:)   !> this rank's blocks' partials, in block order
      real(R8),                intent(out)   :: total
      type(morph_status_t),    intent(out)   :: status
    end subroutine osum_i
  end interface

  !> One rank: swap and allgather copy, ordered_sum sums in block order.
  type, extends(morph_exchange_t), public :: morph_exchange_serial_t
  contains
    procedure :: rank        => serial_rank
    procedure :: nranks      => serial_nranks
    procedure :: swap        => serial_swap
    procedure :: allgather   => serial_allgather
    procedure :: ordered_sum => serial_ordered_sum
  end type morph_exchange_serial_t

contains

  integer(I4) function serial_rank ( self ) result ( r )
    class(morph_exchange_serial_t), intent(in) :: self
    r = 0
  end function serial_rank


  integer(I4) function serial_nranks ( self ) result ( r )
    class(morph_exchange_serial_t), intent(in) :: self
    r = 1
  end function serial_nranks


  subroutine serial_swap ( self, scount, sbuf, rcount, rbuf, status )
    class(morph_exchange_serial_t), intent(inout) :: self
    integer(I4),                    intent(in)    :: scount(0:)
    real(R8),                       intent(in)    :: sbuf(:)
    integer(I4),                    intent(in)    :: rcount(0:)
    real(R8),                       intent(out)   :: rbuf(:)
    type(morph_status_t),           intent(out)   :: status

    call status%clear()
    if ( size(scount) /= 1 .or. size(rcount) /= 1 ) then
      call status%fail ( MORPH_ERR_EXCHANGE, 'serial exchange: swap counts for more than one rank' )
      return
    endif
    if ( scount(0) /= rcount(0) .or. scount(0) < 0 .or. scount(0) > size(sbuf) .or. rcount(0) > size(rbuf) ) then
      call status%fail ( MORPH_ERR_EXCHANGE, 'serial exchange: swap counts that do not match the buffers' )
      return
    endif
    rbuf(1:rcount(0)) = sbuf(1:scount(0))
  end subroutine serial_swap


  subroutine serial_allgather ( self, mine, all, counts, status )
    class(morph_exchange_serial_t), intent(inout) :: self
    real(R8),                       intent(in)    :: mine(:)
    real(R8),                       intent(out)   :: all(:)
    integer(I4),                    intent(in)    :: counts(0:)
    type(morph_status_t),           intent(out)   :: status

    call status%clear()
    if ( size(counts) /= 1 ) then
      call status%fail ( MORPH_ERR_EXCHANGE, 'serial exchange: allgather counts for more than one rank' )
      return
    endif
    if ( counts(0) /= size(mine) .or. size(all) < counts(0) ) then
      call status%fail ( MORPH_ERR_EXCHANGE, 'serial exchange: allgather counts that do not match the buffers' )
      return
    endif
    all(1:counts(0)) = mine
  end subroutine serial_allgather


  subroutine serial_ordered_sum ( self, nblk, owner, partial, total, status )
    class(morph_exchange_serial_t), intent(inout) :: self
    integer(I4),                    intent(in)    :: nblk
    integer(I4),                    intent(in)    :: owner(:)
    real(R8),                       intent(in)    :: partial(:)
    real(R8),                       intent(out)   :: total
    type(morph_status_t),           intent(out)   :: status
    integer(I4) :: b

    call status%clear()
    total = 0.0_R8
    if ( size(owner) /= nblk .or. size(partial) /= nblk ) then
      call status%fail ( MORPH_ERR_EXCHANGE, 'serial exchange: ordered sum of partials that do not match the blocks' )
      return
    endif
    if ( any( owner /= 0 ) ) then
      call status%fail ( MORPH_ERR_EXCHANGE, 'serial exchange: a block on a rank other than 0' )
      return
    endif
    do b = 1, nblk                      ! block order, one at a time: the bits are fixed
      total = total + partial(b)
    enddo
  end subroutine serial_ordered_sum

end module Morph_Exchange
