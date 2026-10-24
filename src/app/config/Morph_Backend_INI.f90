!>@brief Reading input.ini into the registry, and the checks that catch what
!>       the registry cannot see: FUSS's src/lib/config/Backend_INI.f90 at
!>       2943354 (Load_Ini, Check_Unknown_Options, Check_Value_Lengths), with
!>       FUSS's mid-run re-read and its multigrid and probe scans left out, and
!>       one check FUSS does not have (Check_Value_Reads).
module Morph_Backend_INI_m
  use iso_fortran_env, only: R8 => real64
  implicit none
  private
  public :: Load_Ini, Check_Unknown_Options, Check_Value_Lengths, Check_Value_Reads

contains

  !> Assign every registered parameter that input.ini sets, and mark it set.
  subroutine Load_Ini ( fini )
    use Morph_Input_Registry, only: reg
    use Finer,                only: file_ini
    implicit none
    type(file_ini), intent(inout) :: fini
    ! Local
    integer :: i, error

    do i = 1, reg%size
      error = 1

      ! Each registry entry has exactly one associated typed pointer.
      if (associated(reg%params(i)%value%i)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%i, error=error)
      else if (associated(reg%params(i)%value%r)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%r, error=error)
      else if (associated(reg%params(i)%value%l)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%l, error=error)
      else if (associated(reg%params(i)%value%s)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%s, error=error)
      else if (associated(reg%params(i)%value%iarr)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%iarr, error=error)
      else if (associated(reg%params(i)%value%rarr)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, val=reg%params(i)%value%rarr, error=error)
      end if

      if (error == 0) reg%params(i)%is_set = .true.
    end do

  end subroutine Load_Ini


  !> Report input.ini content that MORPH would silently not apply: a key in a
  !> [MORPH-*] section that nothing registered (a misspelt key otherwise falls
  !> back to its default without a word), and a repeated [MORPH-*] section
  !> (FiNeR keeps only the first occurrence). Sections of other tools are not
  !> MORPH's business and are not looked at.
  subroutine Check_Unknown_Options ( fini, out )
    use Finer,                only: file_ini
    use Morph_Input_Registry, only: reg
    use Morph_Config_Types_m, only: codename
    implicit none
    type(file_ini),   intent(inout) :: fini
    character(len=*), intent(out)   :: out
    ! Local
    character(len=:), allocatable :: sections(:), pairs(:)
    character(len=256) :: sec, opt
    integer :: s, i, nbad
    logical :: known

    out   = ''
    nbad  = 0

    call fini%get_sections_list ( sections )
    if ( .not. allocated(sections) ) return

    do s = 1, size(sections)
      sec = adjustl(sections(s))
      if ( .not. Section_Is_Ours(sec) ) cycle
      do i = 1, s - 1
        if ( Same(sections(i), sec) ) then
          nbad = nbad + 1
          if ( nbad == 1 ) then
            out = '[ERROR] repeated section in input.ini: ['//trim(sec)//']'
          else
            out = trim(out)//' ; repeated section ['//trim(sec)//']'
          endif
          exit
        endif
      enddo
    enddo

    do s = 1, size(sections)
      sec = adjustl(sections(s))
      if ( .not. Section_Is_Ours(sec) ) cycle

      ! Each section's loop must be run to completion: FiNeR keeps the loop cursor
      ! in the section object (v2.1.0+, hence intent(inout) on fini), so
      ! abandoning one part-way would corrupt the next.
      do while ( fini%loop( section_name=trim(sec), option_pairs=pairs ) )
        opt = adjustl(pairs(1))

        known = .false.
        do i = 1, reg%size
          if ( Same(reg%params(i)%section, sec) .and. Same(reg%params(i)%name, opt) ) then
            known = .true.
            exit
          endif
        enddo

        if ( .not. known ) then
          nbad = nbad + 1
          if ( nbad == 1 ) then
            out = '[ERROR] unknown option in input.ini: ['//trim(sec)//'] '//trim(opt)
          else
            out = trim(out)//' ; ['//trim(sec)//'] '//trim(opt)
          endif
        endif
      enddo
    enddo

  contains

    logical function Section_Is_Ours ( name ) result ( ours )
      character(len=*), intent(in) :: name
      ours = ( index(name, codename//'-') == 1 )
    end function Section_Is_Ours

    !> Trimmed, case-insensitive comparison.
    logical function Same ( a, b ) result ( eq )
      character(len=*), intent(in) :: a, b
      character(len=256) :: x, y
      integer :: k, c
      x = adjustl(a); y = adjustl(b)
      do k = 1, len_trim(x)
        c = iachar(x(k:k)); if (c >= 65 .and. c <= 90) x(k:k) = achar(c+32)
      enddo
      do k = 1, len_trim(y)
        c = iachar(y(k:k)); if (c >= 65 .and. c <= 90) y(k:k) = achar(c+32)
      enddo
      eq = ( trim(x) == trim(y) )
    end function Same

  end subroutine Check_Unknown_Options


  !> A string value longer than the variable it is read into is cut by FiNeR
  !> without a message (FUSS, measured 2026-10-15). Every string key is checked
  !> against the length of its target.
  subroutine Check_Value_Lengths ( fini, out )
    use Finer,                only: file_ini
    use Morph_Input_Registry, only: reg
    implicit none
    type(file_ini),   intent(inout) :: fini
    character(len=*), intent(out)   :: out
    ! Local
    character(len=4096) :: buf
    character(len=16)   :: sn, sl
    integer :: i, error, n, limit, nbad

    out  = ''
    nbad = 0
    do i = 1, reg%size
      if ( .not. associated(reg%params(i)%value%s) ) cycle
      buf = ''
      call fini%get(reg%params(i)%section, reg%params(i)%name, val=buf, error=error)
      if ( error /= 0 ) cycle
      n     = len_trim(adjustl(buf))
      limit = len(reg%params(i)%value%s)
      if ( n <= limit ) cycle
      nbad = nbad + 1
      write(sn,'(I0)') n
      write(sl,'(I0)') limit
      if ( nbad == 1 ) then
        out = '[ERROR] value too long in input.ini: ['//trim(reg%params(i)%section)//'] '// &
              trim(reg%params(i)%name)//' has '//trim(sn)//' characters, the limit is '//trim(sl)
      else
        out = trim(out)//' ; ['//trim(reg%params(i)%section)//'] '//trim(reg%params(i)%name)// &
              ' has '//trim(sn)//', limit '//trim(sl)
      endif
    enddo

  end subroutine Check_Value_Lengths


  !> A value FiNeR cannot convert leaves its key at the default without a word,
  !> and so does a vector with too many values; with too few, only the leading
  !> elements are assigned (FiNeR v2.2.0, finer_option_t.F90 get_a_option). A
  !> scalar given two values keeps the first. Every key input.ini sets is
  !> therefore checked: Load_Ini could read it, and a numeric or logical key has
  !> exactly as many values as its target (1 for a scalar). Strings are not
  !> counted: Check_Value_Lengths covers them. Not in FUSS's Backend_INI.
  subroutine Check_Value_Reads ( fini, out )
    use Finer,                only: file_ini
    use Morph_Input_Registry, only: reg
    implicit none
    type(file_ini),   intent(inout) :: fini
    character(len=*), intent(out)   :: out
    ! Local
    character(len=4096) :: buf
    character(len=512)  :: what
    character(len=16)   :: sn, sm
    integer :: i, error, n, m, nbad

    out  = ''
    nbad = 0
    do i = 1, reg%size
      if ( associated(reg%params(i)%value%s) ) cycle
      buf = ''
      call fini%get(reg%params(i)%section, reg%params(i)%name, val=buf, error=error)
      if ( error /= 0 ) cycle                          ! not in input.ini: the default stands

      m = 1
      if ( associated(reg%params(i)%value%rarr) ) m = size(reg%params(i)%value%rarr)
      if ( associated(reg%params(i)%value%iarr) ) m = size(reg%params(i)%value%iarr)
      n = Count_Tokens ( buf )

      what = ''
      if ( n /= m ) then
        write(sn,'(I0)') n
        write(sm,'(I0)') m
        what = 'has '//trim(sn)//' value(s), '//trim(sm)//' expected'
      else if ( .not. reg%params(i)%is_set ) then
        what = 'cannot be read: "'//trim(adjustl(buf))//'"'
      endif
      if ( len_trim(what) == 0 ) cycle

      nbad = nbad + 1
      if ( nbad == 1 ) then
        out = '[ERROR] bad value in input.ini: ['//trim(reg%params(i)%section)//'] '// &
              trim(reg%params(i)%name)//' '//trim(what)
      else
        out = trim(out)//' ; ['//trim(reg%params(i)%section)//'] '//trim(reg%params(i)%name)//' '//trim(what)
      endif
    enddo

  contains

    !> Number of blank-separated tokens.
    integer function Count_Tokens ( str ) result ( n )
      character(len=*), intent(in) :: str
      integer :: k
      logical :: inside, blank
      n = 0
      inside = .false.
      do k = 1, len_trim(str)
        blank = ( str(k:k) == ' ' .or. str(k:k) == achar(9) )
        if ( .not. blank .and. .not. inside ) n = n + 1
        inside = .not. blank
      enddo
    end function Count_Tokens

  end subroutine Check_Value_Reads

end module Morph_Backend_INI_m
