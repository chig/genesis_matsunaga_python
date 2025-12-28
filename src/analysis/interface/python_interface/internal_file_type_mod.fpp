module internal_file_type_mod
  use fileio_pdb_mod
  use messages_mod
  use dynamic_string_mod
  use error_mod
  implicit none

contains
  subroutine write_pdb_to_string(dst_str, pdb, err)
    implicit none
    character(len=:), allocatable, intent(out) :: dst_str
    type(s_pdb),             intent(in)    :: pdb
    type(s_error),            intent(inout) :: err
    allocate(character(len=1024) :: dst_str)
    dst_str = ' '
    call append_pdb_to_string(dst_str, pdb, err)
  end subroutine write_pdb_to_string

  subroutine append_pdb_to_string(dst_str, pdb, err)
    implicit none
    ! formal arguments
    character(len=:), allocatable, intent(inout) :: dst_str
    type(s_pdb),             intent(in)    :: pdb
    type(s_error),            intent(inout) :: err

    ! local variables
    integer                  :: i, j, len
    integer                  :: iter, num_atoms, num_ters
    character(80)            :: fmt_a, fmt_b, fmt_t
    character(6)             :: crec
    character(4)             :: catm, cres, cstr, cseg
    character(7)             :: nchar_atom
    logical                  :: use_cid

    character(len=80) :: line_buf

    if (.not.allocated(pdb%hetatm)) then
      call error_set(err, ERROR_NOT_ALLOCATED, &
      'Out_Pdb> not allocated: pdb%hetatm')
      return
    end if

    if (.not.allocated(pdb%atom_name)) then
      call error_set(err, ERROR_NOT_ALLOCATED, &
      'Out_Pdb> not allocated: pdb%atom_name')
      return
    end if
    if (.not.allocated(pdb%residue_name)) then
      call error_set(err, ERROR_NOT_ALLOCATED, &
      'Out_Pdb> not allocated: pdb%residue_name')
      return
    end if
    if (.not.allocated(pdb%segment_name)) then
      call error_set(err, ERROR_NOT_ALLOCATED, &
      'Out_Pdb> not allocated: pdb%segment_name')
      return
    end if
    if (.not.allocated(pdb%chain_id)) then
      call error_set(err, ERROR_NOT_ALLOCATED, &
      'Out_Pdb> not allocated: pdb%chain_id')
      return
    end if
    if (.not.allocated(pdb%atom_no)) then
      call error_set(err, ERROR_NOT_ALLOCATED, &
      'Out_Pdb> not allocated: pdb%atom_no')
      return
    end if
    if (.not.allocated(pdb%residue_no)) then
      call error_set(err, ERROR_NOT_ALLOCATED, &
      'Out_Pdb> not allocated: pdb%residue_no')
      return
    end if
    if (.not.allocated(pdb%atom_coord)) then
      call error_set(err, ERROR_NOT_ALLOCATED, &
      'Out_Pdb> not allocated: pdb%atom_coord')
      return
    end if
    if (.not.allocated(pdb%atom_occupancy)) then
      call error_set(err, ERROR_NOT_ALLOCATED, &
      'Out_Pdb> not allocated: pdb%atom_occupancy')
      return
    end if
    if (.not.allocated(pdb%atom_temp_factor)) then
      call error_set(err, ERROR_NOT_ALLOCATED, &
      'Out_Pdb> not allocated: pdb%atom_temp_factor')
      return
    end if

    num_atoms = size(pdb%hetatm)
    num_ters  = size(pdb%ter_line_no)

    if (pdb%atom_col7) then
      if (pdb%res_col6) then
        use_cid = .false.
        fmt_a   = '(a4,a7,1x,a4,1x,a4,i6,3x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a4,a7,6x,a4,i6,48x,a1)'
      else if (pdb%res_col5) then
        use_cid = .false.
        fmt_a   = '(a4,a7,1x,a4,1x,a4,i5,4x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a4,a7,6x,a4,i5,49x,a1)'
      else
        use_cid = .true.
        fmt_a   = '(a4,a7,1x,a4,1x,a4,a1,i4,4x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_b   = '(a6,a5,1x,a4,1x,a4,a1,i5,3x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a4,a7,6x,a4,a1,i4,49x,a1)'
      end if
    else
      if (pdb%res_col6) then
        use_cid = .false.
        fmt_a   = '(a6,a5,1x,a4,1x,a4,i6,3x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a6,a5,6x,a4,i6,48x,a1)'
      else if (pdb%res_col5) then
        use_cid = .false.
        fmt_a   = '(a6,a5,1x,a4,1x,a4,i5,4x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a6,a5,6x,a4,i5,49x,a1)'
      else
        use_cid = .true.
        fmt_a   = '(a6,a5,1x,a4,1x,a4,a1,i4,4x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_b   = '(a6,a5,1x,a4,1x,a4,a1,i5,3x,3f8.3,f6.2,f6.2,6x,a4)'
        fmt_t   = '(a6,a5,6x,a4,a1,i4,49x,a1)'
      end if
    end if

    ! write CRYST record
    if (pdb%cryst_rec) then
      write(line_buf, fmt='(a6,3(f9.3),a42)') 'CRYST1', &
            (pdb%pbc_box(j,j), j=1,3), &
            '  90.00  90.00  90.00 P 1           1       '
      call append_to_dynamic_string_ln(dst_str, line_buf)
    end if

    ! write MODEL record
    if (pdb%model_rec) then
      write(line_buf, fmt='(a5,5x,i4,61x,a1)') 'MODEL', &
            pdb%model_no, ' '
      call append_to_dynamic_string_ln(dst_str, line_buf)
    end if


    iter = 1

    do i = 1, num_atoms

      ! write ATOM/HETATM record
      if (pdb%hetatm_rec) then
         if (.not. pdb%hetatm(i) .or. num_atoms > 99999) then
            crec = 'ATOM  '
         else
            crec = 'HETATM'
         end if
      else
         crec = 'ATOM  '
      end if

      if (pdb%atom_col7) then
        write(nchar_atom,'(i7)') i
      else
        if (i <= 99999) then
          write(nchar_atom,'(i5)') i
        else
          write(nchar_atom,'(a5)') '*****'
        end if
      end if

      read(pdb%atom_name(i), *) cstr
      len = len_trim(cstr)
      if (len < 4) then
        write(catm, fmt='(1x,a3)') cstr
      else
        catm = pdb%atom_name(i)
      end if

      read(pdb%residue_name(i),*) cstr
      len = len_trim(cstr)
      if (len == 2) then
        write(cres, fmt='(1x,a2)') cstr
      else if (len == 1) then
        write(cres, fmt='(2x,a1)') cstr
      else
        cres = cstr
      end if

      if (pdb%segment) then
        cseg = pdb%segment_name(i)
      else
        cseg = '    '
      end if

      if (use_cid) then
        if (pdb%residue_no(i) < 10000) then
          write(line_buf, fmt=fmt_a)          &
                crec,                         &
                nchar_atom,                   &
                catm,                         &
                cres,                         &
                pdb%chain_id(i),              &
                pdb%residue_no(i),            &
                (pdb%atom_coord(j,i), j=1,3), &
                pdb%atom_occupancy(i),        &
                pdb%atom_temp_factor(i),      &
                cseg
          call append_to_dynamic_string_ln(dst_str, line_buf)
        else
          write(line_buf, fmt=fmt_b)              &
                crec,                         &
                nchar_atom,                   &
                catm,                         &
                cres,                         &
                pdb%chain_id(i),              &
                pdb%residue_no(i),            &
                (pdb%atom_coord(j,i), j=1,3), &
                pdb%atom_occupancy(i),        &
                pdb%atom_temp_factor(i),      &
                cseg
          call append_to_dynamic_string_ln(dst_str, line_buf)
        end if
      else
        write(line_buf, fmt=fmt_a)              &
              crec,                         &
              nchar_atom,                   &
              catm,                         &
              cres,                         &
              pdb%residue_no(i),            &
              (pdb%atom_coord(j,i), j=1,3), &
              pdb%atom_occupancy(i),        &
              pdb%atom_temp_factor(i),      &
              cseg
        call append_to_dynamic_string_ln(dst_str, line_buf)
      end if

      ! write TER record
      if (.not. pdb%ter_rec .or. iter > num_ters) then
        cycle
      end if

      if (pdb%ter_line_no(iter) <= pdb%atom_no(i)) then
        iter = iter + 1
        cycle
      end if

      if (i < num_atoms) then
        if (pdb%ter_line_no(iter) > pdb%atom_no(i+1)) &
          cycle
      end if

      if (pdb%atom_col7) then
        write(nchar_atom,'(i7)') pdb%atom_no(i)+1
      else
        if (i <= 99999) then
          write(nchar_atom,'(i5)') pdb%atom_no(i)+1
        else
          write(nchar_atom,'(a5)') '*****'
        end if
      end if
      if (use_cid) then

        write(line_buf, fmt=fmt_t)   &
             'TER   '          , &
             nchar_atom        , &
             cres              , &
             pdb%chain_id(i)   , &
             pdb%residue_no(i) , &
             ' '
        call append_to_dynamic_string_ln(dst_str, line_buf)

      else

        write(line_buf, fmt=fmt_t)   &
             'TER   '          , &
             nchar_atom        , &
             cres              , &
             pdb%residue_no(i) , &
             ' '
        call append_to_dynamic_string_ln(dst_str, line_buf)

      end if

      iter = iter + 1

    end do

    ! write ENDMDL record
    if (pdb%model_rec) then
      write(line_buf, fmt='(a6,69x,a1)') 'ENDMDL', ' '
      call append_to_dynamic_string_ln(dst_str, line_buf)
    end if

    ! write END record
    if (pdb%end_rec) then
      write(line_buf, fmt='(a6,69x,a1)') 'END   ', ' '
      call append_to_dynamic_string_ln(dst_str, line_buf)
    end if

    return

  end subroutine append_pdb_to_string
end module internal_file_type_mod
