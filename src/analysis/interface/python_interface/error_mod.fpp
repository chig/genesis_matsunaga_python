!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> @brief   utilities for error handling
!! @authors CK
!
!--------1---------2---------3---------4---------5---------6---------7---------8
module error_mod
  use iso_c_binding, only: c_int, c_char, c_null_char
  implicit none
  private
  public :: s_error, error_init, error_clear, error_set, error_has, &
            fi_msg_len, error_to_c

  integer,      public, parameter :: ERROR_CODE = 101
  integer(c_int),  parameter :: MSG_LEN_ = 2048

  type :: s_error
     integer :: code = 0
     character(len=:), allocatable :: msg
  end type s_error

contains
  integer(c_int) function fi_msg_len() bind(C, name="fi_msg_len")
    fi_msg_len = MSG_LEN_
  end function

  subroutine error_init(e)
    type(s_error), intent(out) :: e
    e%code = 0
    if (allocated(e%msg)) deallocate(e%msg)
  end subroutine

  subroutine error_clear(e)
    type(s_error), intent(inout) :: e
    e%code = 0
    if (allocated(e%msg)) e%msg = ''
  end subroutine

  subroutine error_set(e, code, text)
    type(s_error), intent(inout) :: e
    integer, intent(in) :: code
    character(*), intent(in) :: text
    e%code = code
    e%msg  = text
  end subroutine

  logical function error_has(e)
    type(s_error), intent(in) :: e
    error_has = (e%code /= 0)
  end function

  subroutine error_to_c(err, status, msg, msglen)
    use iso_fortran_env, only: error_unit
    implicit none
    type(s_error),         intent(in)   :: err
    integer(c_int),       intent(inout) :: status
    character(kind=c_char), intent(inout) :: msg(*)
    integer(c_int),       value         :: msglen

    integer :: i, n
    character(len=:), allocatable       :: s

    status = err%code
    if (msglen > 0) msg(1) = c_null_char

    if (allocated(err%msg)) then
       s = trim(err%msg)
    else
       s = ''
    end if
    n = min(len_trim(s),msglen-1)
    do i = 1, n
      msg(i) = s(i:i)
    end do
    if (msglen > 0) msg(n+1) = c_null_char
    ! write(error_unit, '(A,I0,2A)') 'fi_error :', status, ' ',s
    ! flush(error_unit)
  end subroutine 

end module error_mod

