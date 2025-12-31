!--------1---------2---------3---------4---------5---------6---------7---------8
!
!> Program  rg_main
!! @brief   RG analysis
!! @authors Motoshi Kamiya (MK), Takaharu Mori (TM), Yuji Sugita (YS)
!
!  (c) Copyright 2016 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rg_c_mod
  use, intrinsic :: iso_c_binding
  use s_molecule_c_mod
  use s_trajectories_c_mod
  use rg_impl_mod

  use trajectory_str_mod
  use molecules_str_mod
  use trj_source_mod
  use string_mod
  use error_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
  implicit none

  public :: rg_analysis_c
  public :: rg_analysis_lazy_c

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    rg_analysis_c
  !> @brief        RG analysis (zerocopy, pre-allocated result array)
  !! @authors      Claude Code
  !! @param[in]    mass_ptr        : pointer to mass array (from Python NumPy)
  !! @param[in]    n_atoms         : number of atoms
  !! @param[in]    s_trajes_c      : trajectories C structure
  !! @param[in]    ana_period      : analysis period
  !! @param[in]    analysis_idx    : analysis atom indices (1-indexed)
  !! @param[in]    n_analysis      : number of analysis atoms
  !! @param[in]    mass_weighted   : use mass weighting (0 or 1)
  !! @param[in]    result_ptr      : pointer to pre-allocated result array
  !! @param[in]    result_size     : size of result array
  !! @param[out]   nstru_out       : actual number of structures analyzed
  !! @param[out]   status          : error status
  !! @param[out]   msg             : error message
  !! @param[in]    msglen          : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine rg_analysis_c(mass_ptr, n_atoms, s_trajes_c, ana_period, &
                           analysis_idx, n_analysis, mass_weighted, &
                           result_ptr, result_size, nstru_out, &
                           status, msg, msglen) &
        bind(C, name="rg_analysis_c")
    implicit none

    ! Arguments
    type(c_ptr), value :: mass_ptr
    integer(c_int), value :: n_atoms
    type(s_trajectories_c), intent(in) :: s_trajes_c
    integer(c_int), value :: ana_period
    type(c_ptr), value :: analysis_idx
    integer(c_int), value :: n_analysis
    integer(c_int), value :: mass_weighted
    type(c_ptr), value :: result_ptr
    integer(c_int), value :: result_size
    integer(c_int), intent(out) :: nstru_out
    integer(c_int), intent(out) :: status
    character(kind=c_char), intent(out) :: msg(*)
    integer(c_int), value :: msglen

    ! Local variables
    type(s_error) :: err
    real(wp), pointer :: mass_f(:)
    real(wp), pointer :: result_f(:)
    integer, pointer :: idx_f(:)
    integer, allocatable :: idx_copy(:)
    logical :: use_mass
    integer :: nstru

    ! Initialize
    call error_init(err)
    status = 0
    nstru_out = 0

    ! Validate inputs
    if (n_analysis <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rg_analysis_c: n_analysis must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(mass_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rg_analysis_c: mass_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(result_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rg_analysis_c: result_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (result_size <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rg_analysis_c: result_size must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    ! Create zero-copy view of mass array from Python
    call C_F_POINTER(mass_ptr, mass_f, [n_atoms])

    ! Create zero-copy view of result array from Python
    call C_F_POINTER(result_ptr, result_f, [result_size])

    ! Convert analysis indices from C pointer to Fortran array
    call C_F_POINTER(analysis_idx, idx_f, [n_analysis])
    allocate(idx_copy(n_analysis))
    idx_copy(:) = idx_f(:)

    ! Convert mass_weighted to logical
    use_mass = (mass_weighted /= 0)

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    ! Run analysis
    write(MsgOut,'(A)') '[STEP1] RG Analysis'
    write(MsgOut,'(A)') ' '

    call analyze(mass_f, s_trajes_c, ana_period, &
                 idx_copy, n_analysis, use_mass, result_f, nstru)

    nstru_out = nstru

    ! Cleanup local arrays (mass_f, result_f are views, don't deallocate)
    deallocate(idx_copy)

  end subroutine rg_analysis_c

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    rg_analysis_lazy_c
  !> @brief        RG analysis with lazy DCD loading (memory efficient)
  !! @authors      Claude Code
  !! @param[in]    dcd_filename    : DCD file path (C string)
  !! @param[in]    filename_len    : length of filename
  !! @param[in]    trj_type        : trajectory type (1=COOR, 2=COOR+BOX)
  !! @param[in]    mass_ptr        : pointer to mass array
  !! @param[in]    n_atoms         : number of atoms
  !! @param[in]    ana_period      : analysis period
  !! @param[in]    analysis_idx    : analysis atom indices (1-indexed)
  !! @param[in]    n_analysis      : number of analysis atoms
  !! @param[in]    mass_weighted   : use mass weighting (0 or 1)
  !! @param[in]    result_ptr      : pointer to pre-allocated result array
  !! @param[in]    result_size     : size of result array
  !! @param[out]   nstru_out       : actual number of structures analyzed
  !! @param[out]   dcd_nframe_out  : total frames in DCD
  !! @param[out]   dcd_natom_out   : atoms per frame in DCD
  !! @param[out]   status          : error status
  !! @param[out]   msg             : error message
  !! @param[in]    msglen          : max length of error message
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine rg_analysis_lazy_c(dcd_filename, filename_len, trj_type, &
                                mass_ptr, n_atoms, ana_period, &
                                analysis_idx, n_analysis, mass_weighted, &
                                result_ptr, result_size, &
                                nstru_out, dcd_nframe_out, dcd_natom_out, &
                                status, msg, msglen) &
        bind(C, name="rg_analysis_lazy_c")
    implicit none

    ! Arguments
    character(kind=c_char), intent(in) :: dcd_filename(*)
    integer(c_int), value :: filename_len
    integer(c_int), value :: trj_type
    type(c_ptr), value :: mass_ptr
    integer(c_int), value :: n_atoms
    integer(c_int), value :: ana_period
    type(c_ptr), value :: analysis_idx
    integer(c_int), value :: n_analysis
    integer(c_int), value :: mass_weighted
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
    real(wp), pointer :: mass_f(:)
    real(wp), pointer :: result_f(:)
    integer, pointer :: idx_f(:)
    integer, allocatable :: idx_copy(:)
    logical :: use_mass
    integer :: nstru, frame_status
    integer :: i, iatom, idx
    real(wp) :: com(3), weight, tot_weight, rg

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
    if (n_analysis <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rg_analysis_lazy_c: n_analysis must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(mass_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rg_analysis_lazy_c: mass_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (.not. c_associated(result_ptr)) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rg_analysis_lazy_c: result_ptr is null")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    if (result_size <= 0) then
      call error_set(err, ERROR_INVALID_PARAM, &
                     "rg_analysis_lazy_c: result_size must be positive")
      call error_to_c(err, status, msg, msglen)
      return
    end if

    ! Create zero-copy views of arrays from Python
    call C_F_POINTER(mass_ptr, mass_f, [n_atoms])
    call C_F_POINTER(result_ptr, result_f, [result_size])

    ! Convert analysis indices from C pointer to Fortran array
    call C_F_POINTER(analysis_idx, idx_f, [n_analysis])
    allocate(idx_copy(n_analysis))
    idx_copy(:) = idx_f(:)

    ! Convert mass_weighted to logical
    use_mass = (mass_weighted /= 0)

    ! Set MPI variables for analysis
    my_city_rank = 0
    nproc_city   = 1
    main_rank    = .true.

    ! Initialize lazy DCD source
    write(MsgOut,'(A)') '[STEP1] Initialize Lazy DCD Source for RG'
    write(MsgOut,'(A)') ' '

    call init_source_lazy_dcd(source, trim(filename_f), trj_type, ana_period)

    ! Return DCD info
    dcd_nframe_out = source%dcd_nframe
    dcd_natom_out = source%dcd_natom

    ! Check atom count
    if (source%dcd_natom /= n_atoms) then
      call error_set(err, ERROR_ATOM_COUNT, &
                     "rg_analysis_lazy_c: atom count mismatch")
      call error_to_c(err, status, msg, msglen)
      call finalize_source(source)
      deallocate(idx_copy)
      return
    end if

    ! Main analysis loop
    write(MsgOut,'(A)') '[STEP2] RG Analysis (lazy loading)'
    write(MsgOut,'(A)') ' '

    nstru = 0

    do while (has_more_frames(source))

      ! Get next frame via lazy loading
      call get_next_frame(source, trajectory, frame_status)
      if (frame_status /= 0) exit

      nstru = nstru + 1
      write(MsgOut,*) '      number of structures = ', nstru

      ! Compute center of mass
      com(1:3) = 0.0_wp
      tot_weight = 0.0_wp
      do iatom = 1, n_analysis
        idx = idx_copy(iatom)
        weight = 1.0_wp
        if (use_mass) weight = mass_f(idx)
        com(:) = com(:) + weight * trajectory%coord(:,idx)
        tot_weight = tot_weight + weight
      end do
      com(:) = com(:) / tot_weight

      ! Compute radius of gyration
      rg = 0.0_wp
      do iatom = 1, n_analysis
        idx = idx_copy(iatom)
        weight = 1.0_wp
        if (use_mass) weight = mass_f(idx)
        do i = 1, 3
          rg = rg + weight * ((trajectory%coord(i,idx) - com(i)) ** 2)
        end do
      end do
      rg = sqrt(rg / tot_weight)

      ! Write result to zerocopy array
      if (nstru <= result_size) then
        result_f(nstru) = rg
      end if

      write(MsgOut,'(a,f10.5)') '              RG of analysis atoms = ', rg
      write(MsgOut,*) ''

    end do

    nstru_out = nstru

    ! Output summary
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> RG lazy analysis completed'
    write(MsgOut,'(A)') ''

    ! Cleanup
    call finalize_source(source)
    deallocate(idx_copy)
    if (allocated(trajectory%coord)) deallocate(trajectory%coord)

  end subroutine rg_analysis_lazy_c

end module rg_c_mod
