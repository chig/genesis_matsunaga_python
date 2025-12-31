!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> Program  dr_main
!! @brief   analysis of Drms
!! @authors Chigusa Kobayashi (CK), Daisuke Matsuoka (DM), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module drms_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use drms_impl_mod

  use trajectory_str_mod
  use molecules_str_mod
  use trj_source_mod
  use string_mod
  use error_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  public :: drms_analysis_c
  public :: drms_analysis_lazy_c

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    drms_analysis_c
  !> @brief        DRMS analysis (zerocopy, pre-allocated result array)
  !! @authors      Claude Code
  !! @param[in]    contact_list_ptr : pointer to contact atom pairs (2, n_contact)
  !! @param[in]    contact_dist_ptr : pointer to reference distances
  !! @param[in]    n_contact        : number of contacts
  !! @param[in]    s_trajes_c       : trajectories C structure
  !! @param[in]    ana_period       : analysis period
  !! @param[in]    pbc_correct      : apply PBC correction (0 or 1)
  !! @param[in]    result_ptr       : pointer to pre-allocated result array
  !! @param[in]    result_size      : size of pre-allocated result array
  !! @param[out]   nstru_out        : actual number of structures analyzed
  !! @param[out]   status           : error status
  !! @param[out]   msg              : error message
  !! @param[in]    msglen           : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine drms_analysis_c(contact_list_ptr, contact_dist_ptr, &
                             n_contact, s_trajes_c, ana_period, &
                             pbc_correct, result_ptr, result_size, &
                             nstru_out, status, msg, msglen) &
        bind(C, name="drms_analysis_c")
    implicit none

    ! Arguments
    type(c_ptr), value :: contact_list_ptr
    type(c_ptr), value :: contact_dist_ptr
    integer(c_int), value :: n_contact
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer(c_int), value :: ana_period
    integer(c_int), value :: pbc_correct
    type(c_ptr), value :: result_ptr
    integer(c_int), value :: result_size
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(s_error) :: err
    integer, pointer :: contact_list_f(:,:)
    real(wp), pointer :: contact_dist_f(:)
    real(wp), pointer :: result_f(:)
    logical :: pbc_flag
    integer :: nstru_local

    ! Initialize
    call error_init(err)
    status = 0
    nstru_out = 0

    ! Validate inputs
    if (n_contact <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_c: n_contact must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(contact_list_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_c: contact_list_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(contact_dist_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_c: contact_dist_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(result_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_c: result_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (result_size <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_c: result_size must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    ! Create zero-copy views of arrays from Python
    call C_F_POINTER(contact_list_ptr, contact_list_f, [2, n_contact])
    call C_F_POINTER(contact_dist_ptr, contact_dist_f, [n_contact])
    call C_F_POINTER(result_ptr, result_f, [result_size])

    ! Convert pbc_correct to logical
    pbc_flag = (pbc_correct /= 0)

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    ! Run analysis
    write(MsgOut,'(A)') '[STEP1] DRMS Analysis'
    write(MsgOut,'(A)') ' '

    call analyze(contact_list_f, contact_dist_f, n_contact, &
                 s_trajes_c, ana_period, pbc_flag, &
                 result_f, nstru_local)

    nstru_out = nstru_local

  end subroutine drms_analysis_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    drms_analysis_lazy_c
  !> @brief        DRMS analysis with lazy DCD loading (memory efficient)
  !! @authors      Claude Code
  !! @param[in]    dcd_filename     : DCD file path (C string)
  !! @param[in]    filename_len     : length of filename
  !! @param[in]    trj_type         : trajectory type (1=COOR, 2=COOR+BOX)
  !! @param[in]    contact_list_ptr : pointer to contact atom pairs (2, n_contact)
  !! @param[in]    contact_dist_ptr : pointer to reference distances
  !! @param[in]    n_contact        : number of contacts
  !! @param[in]    n_atoms          : number of atoms (for DCD validation)
  !! @param[in]    ana_period       : analysis period
  !! @param[in]    pbc_correct      : apply PBC correction (0 or 1)
  !! @param[in]    result_ptr       : pointer to pre-allocated result array
  !! @param[in]    result_size      : size of result array
  !! @param[out]   nstru_out        : actual number of structures analyzed
  !! @param[out]   dcd_nframe_out   : total frames in DCD
  !! @param[out]   dcd_natom_out    : atoms per frame in DCD
  !! @param[out]   status           : error status
  !! @param[out]   msg              : error message
  !! @param[in]    msglen           : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine drms_analysis_lazy_c(dcd_filename, filename_len, trj_type, &
                                  contact_list_ptr, contact_dist_ptr, &
                                  n_contact, n_atoms, ana_period, &
                                  pbc_correct, result_ptr, result_size, &
                                  nstru_out, dcd_nframe_out, dcd_natom_out, &
                                  status, msg, msglen) &
        bind(C, name="drms_analysis_lazy_c")
    implicit none

    ! Arguments
    character(kind=c_char), intent(in) :: dcd_filename(*)
    integer(c_int), value :: filename_len
    integer(c_int), value :: trj_type
    type(c_ptr), value :: contact_list_ptr
    type(c_ptr), value :: contact_dist_ptr
    integer(c_int), value :: n_contact
    integer(c_int), value :: n_atoms
    integer(c_int), value :: ana_period
    integer(c_int), value :: pbc_correct
    type(c_ptr), value :: result_ptr
    integer(c_int), value :: result_size
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: dcd_nframe_out
    integer(c_int), intent(out) :: dcd_natom_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(s_error) :: err
    type(s_trj_source) :: source
    type(s_trajectory) :: trajectory
    character(MaxFilename) :: filename_f
    integer, pointer :: contact_list_f(:,:)
    real(wp), pointer :: contact_dist_f(:)
    real(wp), pointer :: result_f(:)
    logical :: pbc_flag
    integer :: nstru, frame_status
    integer :: i, i1, i2
    real(wp) :: d12(3), r12, tmp, dn, drms
    real(wp) :: box(3)

    ! Initialize
    call error_init(err)
    status = 0
    nstru_out = 0
    dcd_nframe_out = 0
    dcd_natom_out = 0

    ! Convert C string to Fortran string
    filename_f = ''
    do i = 1, min(filename_len, MaxFilename)
      if (dcd_filename(i) == c_null_char) exit
      filename_f(i:i) = dcd_filename(i)
    end do

    ! Validate inputs
    if (n_contact <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_lazy_c: n_contact must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(contact_list_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_lazy_c: contact_list_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(contact_dist_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_lazy_c: contact_dist_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(result_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_lazy_c: result_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (result_size <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "drms_analysis_lazy_c: result_size must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    ! Create zero-copy views of arrays from Python
    call C_F_POINTER(contact_list_ptr, contact_list_f, [2, n_contact])
    call C_F_POINTER(contact_dist_ptr, contact_dist_f, [n_contact])
    call C_F_POINTER(result_ptr, result_f, [result_size])

    ! Convert pbc_correct to logical
    pbc_flag = (pbc_correct /= 0)

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    ! Initialize lazy DCD source
    write(MsgOut,'(A)') '[STEP1] Initialize Lazy DCD Source for DRMS'
    write(MsgOut,'(A)') ' '

    call init_source_lazy_dcd(source, trim(filename_f), trj_type, ana_period)

    ! Return DCD info
    dcd_nframe_out = source%dcd_nframe
    dcd_natom_out = source%dcd_natom

    ! Check atom count
    if (source%dcd_natom /= n_atoms) then
      call error_set(err, ERROR_ATOM_COUNT, &
                     "drms_analysis_lazy_c: atom count mismatch")
      call error_to_c(err, status, msg, msglen)
      call finalize_source(source)
      return
    end if

    ! Main analysis loop
    write(MsgOut,'(A)') '[STEP2] DRMS Analysis (lazy loading)'
    write(MsgOut,'(A)') ' '

    nstru = 0
    dn = real(n_contact, wp)

    do while (has_more_frames(source))

      ! Get next frame via lazy loading
      call get_next_frame(source, trajectory, frame_status)
      if (frame_status /= 0) exit

      nstru = nstru + 1
      write(MsgOut,*) '      number of structures = ', nstru

      ! Get box size for PBC
      box(1) = trajectory%pbc_box(1,1)
      box(2) = trajectory%pbc_box(2,2)
      box(3) = trajectory%pbc_box(3,3)

      ! Compute DRMS
      drms = 0.0_wp

      !$omp parallel do default(none)                          &
      !$omp private(i, d12, r12, tmp, i1, i2)                  &
      !$omp shared(trajectory, contact_list_f, contact_dist_f, n_contact, box, pbc_flag) &
      !$omp reduction(+:drms)
      do i = 1, n_contact
        i1 = contact_list_f(1, i)
        i2 = contact_list_f(2, i)
        d12(1:3) = trajectory%coord(1:3, i1) - trajectory%coord(1:3, i2)
        if (pbc_flag) &
          d12(1:3) = d12(1:3) - anint(d12(1:3) / box(1:3)) * box(1:3)
        r12 = sqrt(d12(1)*d12(1) + d12(2)*d12(2) + d12(3)*d12(3))
        tmp = r12 - contact_dist_f(i)
        drms = drms + tmp * tmp
      end do
      !$omp end parallel do

      if (dn > EPS) drms = sqrt(drms / dn)

      ! Write result to zerocopy array
      if (nstru <= result_size) then
        result_f(nstru) = drms
      end if

      write(MsgOut,'(a,f10.5)') '              DRMS = ', drms
      write(MsgOut,*) ''

    end do

    nstru_out = nstru

    ! Output summary
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> DRMS lazy analysis completed'
    write(MsgOut,'(A)') ''

    ! Cleanup
    call finalize_source(source)
    if (allocated(trajectory%coord)) deallocate(trajectory%coord)

  end subroutine drms_analysis_lazy_c

end module drms_c_mod
