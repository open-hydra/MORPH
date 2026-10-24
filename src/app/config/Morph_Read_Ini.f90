!>@brief input.ini -> registry -> validated input objects (FUSS's
!>       Read_Inifile, src/lib/config/Read_Ini.f90, in one pass).
module Morph_Read_Ini_m

  implicit none
  private
  public :: Read_Inifile

contains

  !> Read and validate `filename`. On return `message` is empty, or it holds
  !> every problem found: the caller reports it and stops. Nothing stops here.
  subroutine Read_Inifile ( filename, message )
    use Finer,                only: file_ini
    use Morph_Register_m,     only: Register_Morph
    use Morph_Backend_INI_m,  only: Load_Ini, Check_Unknown_Options, Check_Value_Lengths, Check_Value_Reads
    use Morph_Input_Registry, only: Validate_Registry
    implicit none
    character(len=*), intent(in)  :: filename
    character(len=*), intent(out) :: message
    ! Local
    type(file_ini)      :: fini
    character(len=1024) :: unknown, long, unread, invalid
    logical             :: present

    message = ''

    inquire( file=trim(filename), exist=present )
    if ( .not. present ) then
      message = '[ERROR] '//trim(filename)//' not found in the working directory'
      return
    endif

    call fini%load( filename=trim(filename) )

    ! Build the registry, then load the values from the file
    call Register_Morph()
    call Load_Ini( fini )

    call Check_Unknown_Options( fini, unknown )
    call Check_Value_Lengths( fini, long )
    call Check_Value_Reads( fini, unread )
    invalid = Validate_Registry()

    call add( unknown )
    call add( long )
    call add( unread )
    call add( invalid )

  contains

    subroutine add ( piece )
      character(len=*), intent(in) :: piece
      if ( len_trim(piece) == 0 ) return
      if ( len_trim(message) == 0 ) then
        message = trim(piece)
      else
        message = trim(message)//new_line('a')//trim(piece)
      endif
    end subroutine add

  end subroutine Read_Inifile

end module Morph_Read_Ini_m
