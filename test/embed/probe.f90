!>@brief MORPH as a solver sees it: one block, Morph_Domain_m and Morph_API
!>       only, the unit cube translated by a quarter of its side.
program probe
  use iso_fortran_env, only: R8 => real64
  use Morph_API,       only: morph_status_t, morph_motion_translation_t
  use Morph_Domain_m,  only: morph_domain_t
  implicit none
  type(morph_domain_t) :: d
  type(morph_status_t) :: st
  integer  :: ib, ia, ir, i, j, k
  real(R8) :: r_abs, r_rel

  call d%allocate_blocks ( reshape([1, 1, 1], [3, 1]), st )
  if ( .not. st%ok() ) error stop 1
  do k = 0, 1; do j = 0, 1; do i = 0, 1
    d%blk(1)%node(i,j,k)%c = [real(i,R8), real(j,R8), real(k,R8)]
  enddo; enddo; enddo
  call d%init ( 0.0_R8, st, ib )
  if ( .not. st%ok() ) error stop 1
  if ( d%blk(1)%geom%vol(1,1,1) /= 1.0_R8 ) error stop 2

  deallocate( d%blk(1)%law )
  allocate( d%blk(1)%law, source = morph_motion_translation_t( name='translation', vel=[0.25_R8, 0.0_R8, 0.0_R8] ) )
  call d%advance ( 1.0_R8, st, ib )
  if ( .not. st%ok() ) error stop 3
  call d%gcl_residual ( r_abs, r_rel, ia, ir )
  if ( r_rel > 1.0e-14_R8 .or. d%blk(1)%geom%vol(1,1,1) /= 1.0_R8 ) error stop 4

  write(*,'(A,ES10.3)') 'probe: unit cube translated, volume 1 exactly, GCL residual ', r_rel
end program probe
