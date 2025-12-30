!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   crd_convert_convert
!> @brief   convert trajectory files
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module crd_convert_impl_mod

  use cc_option_mod
  use cc_option_str_mod
  use pbc_correct_mod
  use fitting_mod
  use fitting_str_mod
  use trajectory_str_mod
  use output_str_mod
  use select_atoms_mod
  use molecules_mod
  use molecules_str_mod
  use constants_mod
  use measure_mod
  use fileio_trj_mod
  use fileio_mod
  use messages_mod
  use error_mod
  use string_mod
  use fileio_pdb_mod
  use s_trajectories_c_mod

  implicit none
  private

  ! subroutines
  public  :: convert
  public  :: get_info
  public  :: convert_zerocopy
  private :: centering
  private :: output_split_trjpdb
  private :: get_filename

  real(wp), allocatable, save   :: tmp_coord(:,:)

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    convert
  !> @brief        convert trajectory files
  !! @authors      NT
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trj_list   : trajectory list information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] fitting    : fitting information
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !! @param[inout] s_trajs_c_array : output multi flame trajectories
  !! @param[out]   num_trajs       : length of s_trajs_c_array
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

 subroutine convert(molecule,   &
                    trj_list,   &
                    trajectory, &
                    fitting,    &
                    option,     &
                    output,     &
                    s_trajs_c_array, &
                    num_trajs,  &
                    err)
    use, intrinsic :: iso_c_binding

    ! formal arguments
    type(s_molecule),   intent(inout) :: molecule
    type(s_trj_list),   intent(inout) :: trj_list
    type(s_trajectory), intent(inout) :: trajectory
    type(s_fitting),    intent(inout) :: fitting
    type(s_option),     intent(inout) :: option
    type(s_output),     intent(inout) :: output
    type(c_ptr), intent(inout) :: s_trajs_c_array
    integer(c_int), intent(out) :: num_trajs
    type(s_error),                   intent(inout) :: err

    ! local variables
    type(s_trj_file)         :: trj_in, trj_out
    integer                  :: nstru, irun, itrj
    integer                  :: rms_out, trr_out
    type(s_trajectories_c), pointer :: s_trajs_c_buf(:)


    ! check-only
    if (option%check_only) &
      return

    ! open output files
    ! if (output%trjfile /= '' .and. .not. option%split_trjpdb) &
    !   call open_trj (trj_out,              &
    !                  output%trjfile,       &
    !                  option%trjout_format, &
    !                  option%trjout_type,   &
    !                  IOFileOutputNew)

    if (output%rmsfile /= '') &
      call open_file(rms_out, output%rmsfile, IOFileOutputNew)

    if (output%trrfile /= '') &
      call open_file(trr_out, output%trrfile, IOFileOutputNew)

    nstru = 0

    ! init s_trajs_c_array
    num_trajs = size(trj_list%md_steps)
    allocate(s_trajs_c_buf(num_trajs))

    do irun = 1, size(trj_list%md_steps)

      call open_trj(trj_in,                   &
                    trj_list%filenames(irun), &
                    trj_list%trj_format,      &
                    trj_list%trj_type,        &
                    IOFileInput)

      do itrj = 1, trj_list%md_steps(irun)

        ! input trj
        !
        call read_trj(trj_in, trajectory)

        if (itrj == 1) then
            ! Initialize trajectory with selected atom count
            if (allocated(option%trjout_atom%idx)) then
              call init_empty_s_trajectories_c( &
                  s_trajs_c_buf(irun), size(option%trjout_atom%idx), sum(trj_list%md_steps))
            else
              call init_empty_s_trajectories_c( &
                  s_trajs_c_buf(irun), size(trajectory%coord, dim=2), sum(trj_list%md_steps))
            end if
        end if


        if (mod(itrj, trj_list%ana_periods(irun)) == 0) then

          nstru = nstru + 1
          write(MsgOut,*) '      number of structures = ', nstru

          ! selection
          !
          call reselect_atom(molecule, &
                             option%trjout_atom_exp, &
                             trajectory%coord, &
                             option%trjout_atom, &
                             option%trjout_atom_trj)

          ! centering
          !
          call centering(molecule, trajectory%coord, option)

          ! pbc-correct
          !
          call run_pbc_correct(option%pbcc_mode, &
                               molecule, &
                               trajectory)

          ! fitting
          !
          if (fitting%mass_weight) then
            call run_fitting(fitting, &
                             molecule%atom_coord, &
                             trajectory%coord, &
                             trajectory%coord, &
                             molecule%mass)

          else
            call run_fitting(fitting, &
                             molecule%atom_coord, &
                             trajectory%coord, &
                             trajectory%coord)

          end if

          ! write data
          !
          if (output%rmsfile /= '') &
            call out_rmsd (rms_out, nstru, fitting)

          if (output%trrfile /= '') &
            call out_trrot(trr_out, nstru, fitting)

          if (output%trjfile /= '') then
            if (option%split_trjpdb) then
              call output_split_trjpdb(nstru, molecule, trajectory, option, output, err)
            else
              ! call write_trj(trj_out, trajectory, option%trjout_atom_trj,molecule)
            end if
          end if

          ! Filter trajectory coordinates based on selected atoms
          if (allocated(option%trjout_atom%idx)) then
            call set_frame_filtered(s_trajs_c_buf(irun), trajectory, itrj, option%trjout_atom)
          else
            call set_frame(s_trajs_c_buf(irun), trajectory, itrj)
          end if
        end if

      end do

      call close_trj(trj_in)

    end do

    if (output%trrfile /= '') call close_file(trr_out)
    if (output%rmsfile /= '') call close_file(rms_out)
    ! if (output%trjfile /= '' .and. .not. option%split_trjpdb) &
    !                           call close_trj (trj_out)

    s_trajs_c_array = c_loc(s_trajs_c_buf)
    return

  end subroutine convert

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    centering
  !> @brief        move the COM of the fitting target to the origin
  !! @authors      DM
  !! @param[in]    molecule   : molecule information
  !! @param[inout] coord      : atom coordinates
  !! @param[in]    fitting    : fitting information
  !! @param[in]    option     : option information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine centering(molecule, coord, option)

  ! formal arguments
  type(s_molecule), intent(in)    :: molecule
  real(wp),         intent(inout) :: coord(:,:)
  type(s_option),   intent(in)    :: option

  ! local variables
  integer  :: iatm, natm
  real(wp) :: com(3)

  if (.not. option%centering) return

  natm = size(molecule%atom_no)
  com  = compute_com(coord, molecule%mass, option%centering_atom%idx)

  do iatm = 1, natm
    coord(:,iatm) = coord(:,iatm) - com(:) + option%center_coord(:)
  end do

  return

  end subroutine centering

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_split_trjpdb
  !> @brief        output split PDB files as trajectory
  !! @authors      TM
  !! @param[inout] nstru      : structure index
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_split_trjpdb(nstru, molecule, trajectory, option, output, err)

    ! formal arguments
    integer,                 intent(inout) :: nstru
    type(s_molecule),        intent(inout) :: molecule
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_option),          intent(inout) :: option
    type(s_output),          intent(inout) :: output
    type(s_error),           intent(inout) :: err

    ! local variables
    type(s_pdb)              :: tmp_pdb


    ! allocate tmp_coord to save molecule%atom_coord
    !
    if(.not. allocated(tmp_coord)) then
      allocate(tmp_coord(3,molecule%num_atoms))
    end if

    tmp_coord(:,:) = molecule%atom_coord(:,:)
    molecule%atom_coord(:,:) = trajectory%coord(:,:)

    call export_molecules(molecule, option%trjout_atom, tmp_pdb)

    if (option%trjout_type == TrjTypeCoorBox) then
      tmp_pdb%cryst_rec = .true.
      tmp_pdb%pbc_box(1:3,1:3) = trajectory%pbc_box(1:3,1:3)
    else
      tmp_pdb%cryst_rec = .false.
    end if

    tmp_pdb%model_rec = .false.

    call output_pdb(get_filename(output%trjfile,nstru, err), tmp_pdb)

    molecule%atom_coord(:,:) = tmp_coord(:,:)

    return

  end subroutine output_split_trjpdb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_filename
  !> @brief        insert snapshot index into {} in the filename
  !! @authors      TM
  !! @param[in]    filename      : filename
  !! @param[in]    no            : index
  !! @note         this subroutine was originally made by NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_filename(filename, no, err)

    ! return
    character(Maxfilename)   :: get_filename

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: no
    type(s_error),           intent(inout) :: err

    ! local variables
    integer                  :: bl, br
    character(100)           :: fid


    bl = index(filename, '{', back=.true.)
    br = index(filename, '}', back=.true.)

    if (bl == 0 .or. br == 0 .or. bl > br) then
      call error_set(err, ERROR_INVALID_PARAM, &
                    'Get_Filename> {} is not found in the output trjfile name')
      return
    endif

    write(fid,'(i0)') no
    get_filename = filename(:bl-1) // trim(fid) // filename(br+1:)

    return

  end function get_filename

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_info
  !> @brief        Get trajectory info without loading data
  !! @authors      YS
  !! @param[in]    molecule       : molecule information
  !! @param[in]    trj_filenames  : trajectory file paths (packed)
  !! @param[in]    n_trj_files    : number of trajectory files
  !! @param[in]    filename_len   : max length per filename
  !! @param[in]    trj_format     : trajectory format (TrjFormatDCD etc)
  !! @param[in]    trj_type       : trajectory type (TrjTypeCoor etc)
  !! @param[out]   frame_counts   : frame counts per trajectory
  !! @param[out]   n_trajs        : number of trajectories
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_info(molecule, trj_filenames, n_trj_files, filename_len, &
                      trj_format, trj_type, frame_counts, n_trajs, err)
    use, intrinsic :: iso_c_binding

    ! formal arguments
    type(s_molecule),   intent(in)    :: molecule
    character(kind=c_char), intent(in) :: trj_filenames(*)
    integer(c_int), value, intent(in) :: n_trj_files
    integer(c_int), value, intent(in) :: filename_len
    integer(c_int), value, intent(in) :: trj_format
    integer(c_int), value, intent(in) :: trj_type
    integer(c_int), intent(out)       :: frame_counts(n_trj_files)
    integer(c_int), intent(out)       :: n_trajs
    type(s_error),        intent(inout) :: err

    ! local variables
    integer :: i, offset
    character(len=MaxFilename) :: filename_f

    n_trajs = n_trj_files

    ! Loop through trajectory files to get frame counts
    do i = 1, n_trj_files
      ! Extract filename from packed string
      offset = (i - 1) * filename_len
      call c2f_string_at_offset(trj_filenames, offset, filename_len, filename_f)

      ! Get number of frames from trajectory header
      frame_counts(i) = get_num_steps_trj(trim(filename_f), trj_format, trj_type)

      if (frame_counts(i) <= 0) then
        call error_set(err, ERROR_FILE_FORMAT, &
                       'get_info> Cannot read frame count from: ' // trim(filename_f))
        return
      end if
    end do

  end subroutine get_info

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    convert_zerocopy
  !> @brief        Convert trajectory files with zerocopy pattern
  !! @authors      YS
  !! @param[in]    molecule           : molecule information
  !! @param[in]    trj_filenames      : trajectory file paths (packed)
  !! @param[in]    n_trj_files        : number of trajectory files
  !! @param[in]    filename_len       : max length per filename
  !! @param[in]    trj_format         : trajectory format
  !! @param[in]    trj_type           : trajectory type
  !! @param[in]    selected_indices   : selected atom indices (1-indexed)
  !! @param[in]    n_selected         : number of selected atoms
  !! @param[in]    fitting_method     : fitting method
  !! @param[in]    fitting_indices    : fitting atom indices (1-indexed)
  !! @param[in]    n_fitting          : number of fitting atoms
  !! @param[in]    mass_weighted      : use mass weighting for fitting
  !! @param[in]    do_centering       : enable centering
  !! @param[in]    centering_indices  : centering atom indices
  !! @param[in]    n_centering        : number of centering atoms
  !! @param[in]    center_coord       : target center coordinates
  !! @param[in]    pbcc_mode          : PBC correction mode
  !! @param[in]    ana_period         : analysis period
  !! @param[in]    frame_counts       : frame counts per trajectory
  !! @param[inout] coords_ptrs        : pre-allocated coords arrays
  !! @param[inout] pbc_box_ptrs       : pre-allocated pbc_box arrays
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine convert_zerocopy(molecule, &
                              trj_filenames, n_trj_files, filename_len, &
                              trj_format, trj_type, &
                              selected_indices, n_selected, &
                              fitting_method, fitting_indices, n_fitting, mass_weighted, &
                              do_centering, centering_indices, n_centering, center_coord, &
                              pbcc_mode, ana_period, &
                              frame_counts, &
                              coords_ptrs, pbc_box_ptrs, &
                              err)
    use, intrinsic :: iso_c_binding
    use fitting_mod
    use fitting_str_mod
    use select_atoms_str_mod

    ! formal arguments
    type(s_molecule),   intent(inout) :: molecule
    character(kind=c_char), intent(in) :: trj_filenames(*)
    integer(c_int), value, intent(in) :: n_trj_files
    integer(c_int), value, intent(in) :: filename_len
    integer(c_int), value, intent(in) :: trj_format
    integer(c_int), value, intent(in) :: trj_type
    integer(c_int), intent(in)        :: selected_indices(*)
    integer(c_int), value, intent(in) :: n_selected
    integer(c_int), value, intent(in) :: fitting_method
    integer(c_int), intent(in)        :: fitting_indices(*)
    integer(c_int), value, intent(in) :: n_fitting
    integer(c_int), value, intent(in) :: mass_weighted
    integer(c_int), value, intent(in) :: do_centering
    integer(c_int), intent(in)        :: centering_indices(*)
    integer(c_int), value, intent(in) :: n_centering
    real(c_double), intent(in)        :: center_coord(3)
    integer(c_int), value, intent(in) :: pbcc_mode
    integer(c_int), value, intent(in) :: ana_period
    integer(c_int), intent(in)        :: frame_counts(n_trj_files)
    type(c_ptr), intent(in)           :: coords_ptrs(n_trj_files)
    type(c_ptr), intent(in)           :: pbc_box_ptrs(n_trj_files)
    type(s_error),        intent(inout) :: err

    ! local variables
    type(s_trj_file)    :: trj_in
    type(s_trajectory)  :: trajectory
    type(s_fitting)     :: fitting
    integer :: i, j, k, irun, itrj, frame_idx, offset
    integer :: natom_full
    character(len=MaxFilename) :: filename_f
    real(wp), pointer :: coords_f(:,:,:)   ! (3, n_selected, nframes)
    real(wp), pointer :: pbc_box_f(:,:,:)  ! (3, 3, nframes)
    real(wp) :: com(3)
    integer, allocatable :: sel_idx(:), fit_idx(:), cen_idx(:)

    ! Allocate temporary arrays for indices
    allocate(sel_idx(n_selected))
    do i = 1, n_selected
      sel_idx(i) = selected_indices(i)
    end do

    if (n_fitting > 0) then
      allocate(fit_idx(n_fitting))
      do i = 1, n_fitting
        fit_idx(i) = fitting_indices(i)
      end do
    end if

    if (n_centering > 0) then
      allocate(cen_idx(n_centering))
      do i = 1, n_centering
        cen_idx(i) = centering_indices(i)
      end do
    end if

    ! Allocate trajectory structure for reading (uses full molecule atom count)
    natom_full = molecule%num_atoms
    call alloc_trajectory(trajectory, natom_full)

    ! Setup fitting structure if needed
    if (fitting_method /= FittingMethodNO .and. n_fitting > 0) then
      call alloc_selatoms(fitting%fitting_atom, n_fitting)
      fitting%fitting_method = fitting_method
      fitting%fitting_atom%idx(1:n_fitting) = fit_idx(1:n_fitting)
      fitting%mass_weight = (mass_weighted /= 0)
    end if

    ! Process each trajectory file
    do irun = 1, n_trj_files
      ! Extract filename
      offset = (irun - 1) * filename_len
      call c2f_string_at_offset(trj_filenames, offset, filename_len, filename_f)

      ! Open trajectory file
      call open_trj(trj_in, trim(filename_f), trj_format, trj_type, IOFileInput)

      ! Get Fortran pointers to pre-allocated Python arrays
      call c_f_pointer(coords_ptrs(irun), coords_f, [3, n_selected, frame_counts(irun)])
      call c_f_pointer(pbc_box_ptrs(irun), pbc_box_f, [3, 3, frame_counts(irun)])

      frame_idx = 0

      ! Read frames
      do itrj = 1, frame_counts(irun)
        call read_trj(trj_in, trajectory)

        ! Apply ana_period filter
        if (mod(itrj, ana_period) /= 0) cycle

        frame_idx = frame_idx + 1

        ! Centering
        if (do_centering /= 0 .and. n_centering > 0) then
          natom_full = size(trajectory%coord, dim=2)
          com = compute_com(trajectory%coord, molecule%mass, cen_idx)
          do i = 1, natom_full
            trajectory%coord(:,i) = trajectory%coord(:,i) - com(:) + center_coord(:)
          end do
        end if

        ! PBC correction
        if (pbcc_mode /= PBCCModeNO) then
          call run_pbc_correct(pbcc_mode, molecule, trajectory)
        end if

        ! Fitting
        if (fitting_method /= FittingMethodNO .and. n_fitting > 0) then
          if (fitting%mass_weight) then
            call run_fitting(fitting, molecule%atom_coord, trajectory%coord, &
                           trajectory%coord, molecule%mass)
          else
            call run_fitting(fitting, molecule%atom_coord, trajectory%coord, &
                           trajectory%coord)
          end if
        end if

        ! Copy selected atoms to output buffer
        do j = 1, n_selected
          k = sel_idx(j)
          coords_f(:, j, frame_idx) = trajectory%coord(:, k)
        end do

        ! Copy PBC box
        pbc_box_f(:, :, frame_idx) = trajectory%pbc_box(:, :)
      end do

      call close_trj(trj_in)
    end do

    ! Cleanup
    if (allocated(sel_idx)) deallocate(sel_idx)
    if (allocated(fit_idx)) deallocate(fit_idx)
    if (allocated(cen_idx)) deallocate(cen_idx)
    if (fitting_method /= FittingMethodNO .and. n_fitting > 0) then
      call dealloc_fitting(fitting)
    end if
    call dealloc_trajectory(trajectory)

  end subroutine convert_zerocopy

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    c2f_string_at_offset
  !> @brief        Convert C string at offset in packed array to Fortran string
  !! @param[in]    c_str_array  : packed C string array
  !! @param[in]    offset       : byte offset into array
  !! @param[in]    max_len      : max length of each string
  !! @param[out]   f_str        : output Fortran string
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine c2f_string_at_offset(c_str_array, offset, max_len, f_str)
    use, intrinsic :: iso_c_binding

    character(kind=c_char), intent(in) :: c_str_array(*)
    integer, intent(in) :: offset, max_len
    character(len=*), intent(out) :: f_str

    integer :: i, pos

    f_str = ''
    do i = 1, min(max_len, len(f_str))
      pos = offset + i
      if (c_str_array(pos) == c_null_char) exit
      f_str(i:i) = c_str_array(pos)
    end do

  end subroutine c2f_string_at_offset

end module crd_convert_impl_mod
