module define_molecule
  use iso_c_binding
  use conv_f_c_util
  use input_mod
  use molecules_mod
  use fileio_pdb_mod
  use fileio_psf_mod
  use constants_mod      ! Likely contains MaxFilename
  use molecules_str_mod  ! Contains s_molecule definition
  use s_molecule_c_mod
  implicit none

  private
  public :: define_molecule_from_file

  ! Define MaxFilename if it"s not available from constants_mod
  integer, parameter :: MaxFilename = 256  ! Adjust this value as needed

contains
  subroutine define_molecule_from_file( &
          pdb_filename, &
          top_filename, &
          gpr_filename, &
          psf_filename, &
          ref_filename, &
          fit_filename, &
          prmtop_filename, &
          ambcrd_filename, &
          ambref_filename, &
          grotop_filename, &
          grocrd_filename, &
          groref_filename, &
          out_mol_ptr) &
      bind(C, name="define_molecule_from_file")
    use fileio_top_mod
    use fileio_par_mod
    use fileio_gpr_mod
    use fileio_pdb_mod
    use fileio_psf_mod
    use fileio_crd_mod
    use fileio_prmtop_mod
    use fileio_ambcrd_mod
    use fileio_grotop_mod
    use fileio_grocrd_mod
    use fileio_mode_mod
    use iso_c_binding
    use conv_f_c_util
    implicit none
    ! Input parameters
    type(c_ptr), value :: pdb_filename, top_filename, gpr_filename
    type(c_ptr), value :: psf_filename, ref_filename, fit_filename
    type(c_ptr), value :: prmtop_filename, ambcrd_filename, ambref_filename
    type(c_ptr), value :: grotop_filename, grocrd_filename, groref_filename
    character(kind=c_char), pointer :: pdb_ptr(:)
    character(kind=c_char), pointer :: top_ptr(:)
    character(kind=c_char), pointer :: gpr_ptr(:)
    character(kind=c_char), pointer :: psf_ptr(:)
    character(kind=c_char), pointer :: ref_ptr(:)
    character(kind=c_char), pointer :: fit_ptr(:)
    character(kind=c_char), pointer :: prmtop_ptr(:)
    character(kind=c_char), pointer :: ambcrd_ptr(:)
    character(kind=c_char), pointer :: ambref_ptr(:)
    character(kind=c_char), pointer :: grotop_ptr(:)
    character(kind=c_char), pointer :: grocrd_ptr(:)
    character(kind=c_char), pointer :: groref_ptr(:)

    type(c_ptr), value :: out_mol_ptr
    type(s_molecule_c), pointer :: out_mol

    ! Local variables
    type(s_inp_info) :: inp_info
    type(s_pdb) :: pdb
    type(s_top) :: top
    type(s_gpr) :: gpr
    type(s_psf) :: psf
    type(s_pdb) :: ref
    type(s_pdb) :: fit
    type(s_prmtop) :: prmtop
    type(s_ambcrd) :: ambcrd
    type(s_ambcrd) :: ambref
    type(s_grotop) :: grotop
    type(s_grocrd) :: grocrd
    type(s_grocrd) :: groref
    type(s_molecule) :: molecule
    character(len=MaxFilename) :: filename
    logical :: ok


    if (c_associated(pdb_filename)) then
      call c_f_pointer(pdb_filename, pdb_ptr, [MaxFilename])
      call c2f_string(pdb_ptr, filename)

      if (filename(1:1) /= ' ') then
        inp_info%pdbfile = trim(filename)
        call input_files(inp_info, pdb=pdb)
      end if
    endif

    if (c_associated(top_filename)) then
      call c_f_pointer(top_filename, top_ptr, [MaxFilename])
      call c2f_string(top_ptr, filename)
      if (filename(1:1) /= ' ') then
        inp_info%topfile = trim(filename)
        call input_files(inp_info, top=top)
      end if
    endif
    if (c_associated(gpr_filename)) then
      call c_f_pointer(gpr_filename, gpr_ptr, [MaxFilename])
      call c2f_string(gpr_ptr, filename)
      if (filename(1:1) /= ' ') then
        inp_info%gprfile = trim(filename)
        call input_files(inp_info, gpr=gpr)
      end if
    endif

    if (c_associated(psf_filename)) then
      call c_f_pointer(psf_filename, psf_ptr, [MaxFilename])
      call c2f_string(psf_ptr, filename)
      if (filename(1:1) /= ' ') then
        inp_info%psffile = trim(filename)
        call input_files(inp_info, psf=psf)
      end if
    endif
    if (c_associated(ref_filename)) then
      call c_f_pointer(ref_filename, ref_ptr, [MaxFilename])
      call c2f_string(ref_ptr, filename)
      if (filename(1:1) /= ' ') then
        inp_info%reffile = trim(filename)
        call input_files(inp_info, ref=ref)
      end if
    endif
    if (c_associated(fit_filename)) then
      call c_f_pointer(fit_filename, fit_ptr, [MaxFilename])
      call c2f_string(fit_ptr, filename)
      if (filename(1:1) /= ' ') then
        inp_info%fitfile = trim(filename)
        call input_files(inp_info, fit=fit)
      end if
    end if
    if (c_associated(prmtop_filename)) then
      call c_f_pointer(prmtop_filename, prmtop_ptr, [MaxFilename])
      call c2f_string(prmtop_ptr, filename)
      if (filename(1:1) /= ' ') then
        inp_info%prmtopfile = trim(filename)
        call input_files(inp_info, prmtop=prmtop)
      end if
    end if
    if (c_associated(ambcrd_filename)) then
      call c_f_pointer(ambcrd_filename, ambcrd_ptr, [MaxFilename])
      call c2f_string(ambcrd_ptr, filename)
      if (filename(1:1) /= ' ') then
        inp_info%ambcrdfile = trim(filename)
        call input_files(inp_info, ambcrd=ambcrd)
      end if
    end if
    if (c_associated(ambref_filename)) then
      call c_f_pointer(ambref_filename, ambref_ptr, [MaxFilename])
      call c2f_string(ambref_ptr, filename)
      if (filename(1:1) /= ' ') then
        inp_info%ambreffile = trim(filename)
        call input_files(inp_info, ambref=ambref)
      end if
    end if
    if (c_associated(grotop_filename)) then
      call c_f_pointer(grotop_filename, grotop_ptr, [MaxFilename])
      call c2f_string(grotop_ptr, filename)
      if (filename(1:1) /= ' ') then
        inp_info%grotopfile = trim(filename)
        call input_files(inp_info, grotop=grotop)
      end if
    end if
    if (c_associated(grocrd_filename)) then
      call c_f_pointer(grocrd_filename, grocrd_ptr, [MaxFilename])
      call c2f_string(grocrd_ptr, filename)
      if (filename(1:1) /= ' ') then
        inp_info%grocrdfile = trim(filename)
        call input_files(inp_info, grocrd=grocrd)
      end if
    end if
    if (c_associated(groref_filename)) then
      call c_f_pointer(groref_filename, groref_ptr, [MaxFilename])
      call c2f_string(groref_ptr, filename)
      if (filename(1:1) /= ' ') then
        inp_info%groreffile = trim(filename)
        call input_files(inp_info, groref=groref)
      end if
    end if
    call define_molecules(molecule, &
        pdb=pdb, top=top, gpr=gpr, psf=psf, ref=ref, fit=fit, &
        prmtop=prmtop, ambcrd=ambcrd, ambref=ambref, &
        grotop=grotop, grocrd=grocrd, groref=groref)
    call c_f_pointer(out_mol_ptr, out_mol)
    call f2c_s_molecule(molecule, out_mol)

  end subroutine define_molecule_from_file

end module define_molecule
