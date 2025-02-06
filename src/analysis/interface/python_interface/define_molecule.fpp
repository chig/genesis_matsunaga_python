module define_molecule
  use, intrinsic :: iso_c_binding
  use input_mod
  use molecules_mod
  use fileio_pdb_mod
  use fileio_psf_mod
  use constants_mod      ! Likely contains MaxFilename
  use molecules_str_mod  ! Contains s_molecule definition
  use s_molecule_c_mod
  implicit none

  private
  public :: define_molecule_from_pdb
  public :: define_molecule_from_pdb_psf

  ! Define MaxFilename if it"s not available from constants_mod
  integer, parameter :: MaxFilename = 256  ! Adjust this value as needed

contains
  subroutine define_molecule_from_pdb(pdb_filename, out_mol) &
      bind(C, name="define_molecule_from_pdb")
    use conv_f_c_util
    implicit none
    ! Input parameters
    character(kind=c_char), intent(in) :: pdb_filename(*)
    ! Output parameters
    type(s_molecule_c), intent(out) :: out_mol
    ! Local variables
    type(s_inp_info) :: inp_info
    type(s_pdb) :: pdb
    type(s_molecule) :: molecule
    character(MaxFilename) :: filename

    call c2f_string(pdb_filename, filename)
    inp_info%pdbfile = trim(filename)
    call input_files(inp_info, pdb=pdb)
    call define_molecules(molecule, pdb=pdb)
    call f2c_s_molecule(molecule, out_mol)
  end subroutine define_molecule_from_pdb

  subroutine define_molecule_from_pdb_psf(pdb_path, psf_path, out_mol) &
      bind(C, name="define_molecule_from_pdb_psf")
    use conv_f_c_util
    implicit none
    ! Input parameters
    character(kind=c_char), intent(in) :: pdb_path(*)
    character(kind=c_char), intent(in) :: psf_path(*)
    ! Output parameters
    type(s_molecule_c), intent(out) :: out_mol
    ! Local variables
    type(s_inp_info) :: inp_info
    type(s_pdb) :: pdb
    type(s_psf) :: psf
    type(s_molecule) :: molecule
    character(:), allocatable :: fort_pdb_path
    character(:), allocatable :: fort_psf_path
    call c2f_string_allocate(pdb_path, fort_pdb_path)
    call c2f_string_allocate(psf_path, fort_psf_path)
    inp_info%pdbfile = trim(fort_pdb_path)
    inp_info%psffile = trim(fort_psf_path)
    call input_files(inp_info, pdb=pdb, psf=psf)
    call define_molecules(molecule, pdb=pdb, psf=psf)
    call f2c_s_molecule(molecule, out_mol)
  end subroutine define_molecule_from_pdb_psf
end module define_molecule
