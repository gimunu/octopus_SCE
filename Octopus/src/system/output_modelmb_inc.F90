!! Copyright (C) 2009 N. Helbig and M. Verstraete
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id: states.F90 5022 2009-03-03 17:47:58Z nitsche $

  ! ---------------------------------------------------------
  !
  !> routine for output of model many-body quantities.
  !
  subroutine X(output_modelmb) (dir, gr, st, geo, outp)
    type(states_t),         intent(inout) :: st
    type(grid_t),           intent(inout) :: gr
    character(len=*),       intent(in)    :: dir
    type(geometry_t),       intent(in)    :: geo
    type(output_t),         intent(in)    :: outp

    integer :: mm, iunit, itype
    integer :: ierr
    integer :: tdrun
    integer :: ncombo
    integer, allocatable :: ndiagrams(:)
    integer, allocatable :: young_used(:)
    logical :: symmetries_satisfied, impose_exch_symmetry
    R_TYPE, allocatable :: wf(:)
    character(len=80) :: dirname
    character(len=80) :: filename
    type(modelmb_denmat_t) :: denmat
    type(unit_t)  :: fn_unit


    PUSH_SUB(X(output_modelmb))

    impose_exch_symmetry = .true.
    ! this is ugly, but no other way to tell if we are in td run.
    ! other option is to extract present routine and call explicitly from outside output_all. Don`t wanna.
    tdrun = index(dir, 'td.')
    if (tdrun > 0) impose_exch_symmetry = .false.

    ! make sure directory exists
    call io_mkdir(trim(dir))
    ! all model mb stuff should be in this directory
    dirname = trim(dir)//'modelmb'
    call io_mkdir(trim(dirname))

    SAFE_ALLOCATE(wf(1:gr%mesh%np))

    call modelmb_density_matrix_nullify(denmat)
    if(iand(outp%what, C_OUTPUT_MMB_DEN) /= 0) then
      call modelmb_density_matrix_init(dirname, st, denmat)
    end if
 
    ! open file for Young diagrams and projection info
    write (filename,'(a,a)') trim(dirname), '/youngprojections'
    iunit = io_open(trim(filename), action='write')

    ! treat all particle types
    SAFE_ALLOCATE(ndiagrams(1:st%modelmbparticles%ntype_of_particle))
    ndiagrams = 1
    do itype = 1, st%modelmbparticles%ntype_of_particle
      write (iunit, '(a, i6)') '  Young diagrams for particle type ', itype
      call young_write_allspins (iunit, st%modelmbparticles%nparticles_per_type(itype))
      call young_ndiagrams (st%modelmbparticles%nparticles_per_type(itype), ndiagrams(itype))
    end do
 
    ncombo = product(ndiagrams)
    write (iunit, '(a, i6)') ' # of possible combinations of Young diagrams for all types = ', &
      ncombo

    SAFE_ALLOCATE(young_used(1:ncombo))
    young_used = 0


    ! write header
    write (iunit, '(a)') '  state        eigenvalue         projection   nspindown Young# for each type'

    do mm = 1, st%nst
      call states_get_state(st, gr%mesh, 1, mm, 1, wf)

      symmetries_satisfied = .true.
      if (impose_exch_symmetry) then
        if (mm > 1) then
          ! if eigenval is not degenerate reset young_used
          if (abs(st%eigenval(mm,1) - st%eigenval(mm-1,1)) > 1.e-5) then
            young_used = 0
          end if
        end if

        call X(modelmb_sym_state)(st%eigenval(mm,1), iunit, gr, mm, &
          st%modelmbparticles, ncombo, young_used, wf, symmetries_satisfied, .true.)
        write (iunit,'(a,l1)') "symmetries_satisfied1 ", symmetries_satisfied
        young_used = 0

        call X(modelmb_sym_state)(st%eigenval(mm,1), iunit, gr, mm, &
          st%modelmbparticles, ncombo, young_used, wf, symmetries_satisfied, .true.)
        write (iunit,'(a,l1)') "symmetries_satisfied2 ", symmetries_satisfied
      end if

      if(iand(outp%what, C_OUTPUT_MMB_DEN) /= 0 .and. symmetries_satisfied) then
        call X(modelmb_density_matrix_write)(gr, st, wf, mm, denmat)
      end if

      if(gr%mesh%parallel_in_domains) then
      end if

      if(iand(outp%what, C_OUTPUT_MMB_WFS) /= 0 .and. symmetries_satisfied) then
        fn_unit = units_out%length**(-gr%mesh%sb%dim)
        write(filename, '(a,i4.4)') 'wf-st', mm
          call X(io_function_output)(outp%how, trim(dirname), trim(filename), gr%mesh, wf, &
            fn_unit, ierr, is_tmp = .false., geo = geo)
      end if

    end do

    call io_close(iunit)

    SAFE_DEALLOCATE_A(ndiagrams)
    SAFE_DEALLOCATE_A(young_used)

    SAFE_DEALLOCATE_A(wf)

    if(iand(outp%what, C_OUTPUT_MMB_DEN) /= 0) then
      call modelmb_density_matrix_end (denmat)
    end if
 
    POP_SUB(X(output_modelmb))

  end subroutine X(output_modelmb)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
