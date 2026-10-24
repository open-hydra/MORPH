!>@brief Mesh input and output of the driver, through ORION: the only module
!>       that names both ORION's and MORPH's types.
!>
!> Two formats, both ASCII (the driver builds ORION without TecIO):
!>   p3d      PLOT3D multiblock. ORION writes 17 significant digits, so a mesh
!>            survives a write and a read bit for bit.
!>   tecplot  Tecplot structured zones (.tec / .dat), e.g. the INPUT/ic.tec of a
!>            FUSS case; the variables are read and dropped, the mesh is written
!>            alone. ORION writes 15 significant digits (E23.15E3): not exact.
!>
!> ORION stops the program with exit status 0 on input it cannot handle (a
!> binary Tecplot file without TecIO, a binary PLOT3D write), and a stop with
!> status 0 looks like success to a script. Everything ORION would stop on is
!> therefore refused here first, with a message and a non-empty `message`.
module Morph_IO_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use Lib_ORION_data,  only: orion_data
  use Lib_PLOT3D,      only: p3d_read_multiblock, p3d_write_multiblock
  use Lib_Tecplot,     only: tec_read_structured_multiblock, tec_write_structured_multiblock
  use Morph_API,       only: morph_status_t
  use Morph_Domain_m,  only: morph_domain_t

  implicit none
  private
  public :: Read_Mesh, Write_Mesh

contains

  !> Read a multi-block structured mesh into `domain` (blocks allocated, nodes
  !> filled, the static law on every block).
  subroutine Read_Mesh ( filename, format, domain, message )
    character(len=*),     intent(in)    :: filename
    character(len=*),     intent(in)    :: format
    type(morph_domain_t), intent(inout) :: domain
    character(len=*),     intent(out)   :: message
    ! Local
    type(orion_data)         :: orion
    type(morph_status_t)     :: status
    integer(I4), allocatable :: dims(:,:)
    integer(I4)              :: b, i, j, k, d, err
    logical                  :: present
    character(len=16)        :: sb

    message = ''

    inquire( file=trim(filename), exist=present )
    if ( .not. present ) then
      message = '[ERROR] mesh-in '//trim(filename)//' not found'
      return
    endif

    select case ( trim(format) )
    case ( 'p3d' )
      err = p3d_read_multiblock( orion, trim(filename) )
    case ( 'tecplot' )
      if ( .not. ascii_tecplot_name(filename) ) then
        message = '[ERROR] mesh-in '//trim(filename)//': a tecplot mesh must be ASCII, .tec or .dat '// &
                  '(binary Tecplot needs TecIO, which this build does not link)'
        return
      endif
      err = tec_read_structured_multiblock( orion=orion, filename=trim(filename) )
    case default
      message = '[ERROR] mesh-format '//trim(format)//' is not p3d or tecplot'
      return
    end select

    if ( err /= 0 ) then
      write(sb,'(I0)') err
      message = '[ERROR] reading mesh-in '//trim(filename)//' failed (ORION error '//trim(sb)//')'
      return
    endif
    if ( .not. allocated(orion%block) ) then
      message = '[ERROR] mesh-in '//trim(filename)//' holds no block'
      return
    endif
    if ( size(orion%block) < 1 ) then
      message = '[ERROR] mesh-in '//trim(filename)//' holds no block'
      return
    endif

    allocate( dims(3, size(orion%block)) )
    do b = 1, size(orion%block)
      write(sb,'(I0)') b
      if ( .not. allocated(orion%block(b)%mesh) ) then
        message = '[ERROR] mesh-in: block '//trim(sb)//' has no coordinates'
        return
      endif
      ! MORPH works on hexahedra: a 2-D mesh carries two node planes in k with
      ! three coordinates, as FUSS's meshes do. ORION reads a single node plane
      ! (PLOT3D 2-D, a Tecplot zone with K = 1) as two coordinates and Nk = 0.
      if ( size(orion%block(b)%mesh,1) /= 3 .or. orion%block(b)%Nk < 1 ) then
        message = '[ERROR] mesh-in: block '//trim(sb)//' is a single node plane; MORPH needs hexahedra '// &
                  '(a 2-D mesh with two node planes in k and three coordinates, as in FUSS)'
        return
      endif
      if ( orion%block(b)%Ni < 1 .or. orion%block(b)%Nj < 1 ) then
        message = '[ERROR] mesh-in: block '//trim(sb)//' has no cell in i or j'
        return
      endif
      dims(:,b) = [ orion%block(b)%Ni, orion%block(b)%Nj, orion%block(b)%Nk ]
    enddo

    call domain%allocate_blocks ( dims, status )
    if ( .not. status%ok() ) then
      message = '[ERROR] '//trim(status%message)
      return
    endif

    do b = 1, domain%nblocks
      do k = 0, dims(3,b)
      do j = 0, dims(2,b)
      do i = 0, dims(1,b)
        do d = 1, 3
          domain%blk(b)%node(i,j,k)%c(d) = orion%block(b)%mesh(d,i,j,k)
        enddo
      enddo; enddo; enddo
    enddo

  end subroutine Read_Mesh


  !> Write the domain's current nodes, mesh only, in `format`. The directory of
  !> `filename` is created when it does not exist.
  subroutine Write_Mesh ( filename, format, domain, message )
    character(len=*),     intent(in)  :: filename
    character(len=*),     intent(in)  :: format
    type(morph_domain_t), intent(in)  :: domain
    character(len=*),     intent(out) :: message
    ! Local
    type(orion_data)  :: orion
    integer(I4)       :: b, i, j, k, d, err, slash, stat
    character(len=16) :: sb

    message = ''

    select case ( trim(format) )
    case ( 'p3d' )
    case ( 'tecplot' )
      if ( .not. ascii_tecplot_name(filename) ) then
        message = '[ERROR] mesh-out '//trim(filename)//': a tecplot mesh is written as ASCII, .tec or .dat'
        return
      endif
    case default
      message = '[ERROR] mesh-format '//trim(format)//' is not p3d or tecplot'
      return
    end select

    slash = index( trim(filename), '/', back=.true. )
    if ( slash > 1 ) then
      call execute_command_line( 'mkdir -p "'//filename(1:slash-1)//'"', exitstat=stat )
      if ( stat /= 0 ) then
        message = '[ERROR] cannot create the directory of '//trim(filename)
        return
      endif
    endif

    allocate( orion%block(domain%nblocks) )
    orion%solutiontime = domain%t
    do b = 1, domain%nblocks
      write(sb,'(I0)') b
      orion%block(b)%name = 'block'//trim(sb)
      orion%block(b)%Ni   = domain%blk(b)%dim(1)
      orion%block(b)%Nj   = domain%blk(b)%dim(2)
      orion%block(b)%Nk   = domain%blk(b)%dim(3)
      allocate( orion%block(b)%mesh(1:3, 0:domain%blk(b)%dim(1), 0:domain%blk(b)%dim(2), 0:domain%blk(b)%dim(3)) )
      do k = 0, domain%blk(b)%dim(3)
      do j = 0, domain%blk(b)%dim(2)
      do i = 0, domain%blk(b)%dim(1)
        do d = 1, 3
          orion%block(b)%mesh(d,i,j,k) = domain%blk(b)%node(i,j,k)%c(d)
        enddo
      enddo; enddo; enddo
    enddo

    select case ( trim(format) )
    case ( 'p3d' )
      orion%p3d%format = 'ascii'
      err = p3d_write_multiblock( orion, trim(filename) )
    case default
      orion%tec%format = 'ascii'
      err = tec_write_structured_multiblock( orion=orion, filename=trim(filename) )
    end select

    if ( err /= 0 ) then
      write(sb,'(I0)') err
      message = '[ERROR] writing '//trim(filename)//' failed (ORION error '//trim(sb)//')'
    endif

  end subroutine Write_Mesh


  !> .tec or .dat at the end of the name: what ORION reads and writes as ASCII.
  logical function ascii_tecplot_name ( filename ) result ( ok )
    character(len=*), intent(in) :: filename
    integer :: n
    n  = len_trim(filename)
    ok = .false.
    if ( n > 4 ) ok = ( filename(n-3:n) == '.tec' .or. filename(n-3:n) == '.dat' )
  end function ascii_tecplot_name

end module Morph_IO_m
