!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id: vibrations.F90 11057 2013-08-02 01:48:17Z xavier $

#include "global.h"

module vibrations_m
  use geometry_m
  use global_m
  use io_m
  use lalg_adv_m
  use messages_m
  use mpi_m
  use profiling_m
  use simul_box_m
  use species_m
  use unit_m
  use unit_system_m

  implicit none

  private
  public :: &
       vibrations_t,                    &
       vibrations_init,                 &
       vibrations_end,                  &
       vibrations_symmetrize_dyn_matrix, &
       vibrations_normalize_dyn_matrix, &
       vibrations_out_dyn_matrix,       &
       vibrations_norm_factor,          &
       vibrations_diag_dyn_matrix,      &
       vibrations_get_index,            &
       vibrations_get_atom,             &
       vibrations_get_dir,              &
       vibrations_output,               &
       vibrations_get_suffix
  
  type vibrations_t
    integer :: num_modes
    integer :: ndim
    integer :: natoms
    FLOAT, pointer :: dyn_matrix(:,:), normal_mode(:,:), freq(:)
    FLOAT :: disp
    FLOAT :: total_mass
    integer :: dyn_mat_unit
    character (len=2) :: suffix
  end type vibrations_t

contains

  ! ---------------------------------------------------------
  subroutine vibrations_init(this, geo, sb, suffix)
    type(vibrations_t), intent(out) :: this
    type(geometry_t),   intent(in)  :: geo
    type(simul_box_t),  intent(in)  :: sb
    character (len=2),  intent(in)  :: suffix

    integer :: iatom

    PUSH_SUB(vibrations_init)

    this%ndim = sb%dim
    this%natoms = geo%natoms
    this%num_modes = geo%natoms*sb%dim
    SAFE_ALLOCATE(this%dyn_matrix(1:this%num_modes, 1:this%num_modes))
    SAFE_ALLOCATE(this%normal_mode(1:this%num_modes, 1:this%num_modes))
    SAFE_ALLOCATE(this%freq(1:this%num_modes))

    this%total_mass = M_ZERO
    do iatom = 1, geo%natoms
      this%total_mass = this%total_mass + species_weight(geo%atom(iatom)%spec)
    end do

    this%suffix = suffix
    if(mpi_grp_is_root(mpi_world)) then
      call io_mkdir(VIB_MODES_DIR)
      this%dyn_mat_unit = io_open(VIB_MODES_DIR//'dynamical_matrix_'//trim(this%suffix), action='write')
    endif

    POP_SUB(vibrations_init)
  end subroutine vibrations_init


  ! ---------------------------------------------------------
  subroutine vibrations_end(this)
    type(vibrations_t), intent(inout) :: this

    PUSH_SUB(vibrations_end)

    SAFE_DEALLOCATE_P(this%dyn_matrix)
    SAFE_DEALLOCATE_P(this%freq)
    SAFE_DEALLOCATE_P(this%normal_mode)

    if(mpi_grp_is_root(mpi_world)) call io_close(this%dyn_mat_unit)

    POP_SUB(vibrations_end)
  end subroutine vibrations_end


  ! ---------------------------------------------------------
  character(len=2) function vibrations_get_suffix(this)
    type(vibrations_t), intent(in) :: this

    PUSH_SUB(vibrations_get_suffix)
    vibrations_get_suffix = this%suffix

    POP_SUB(vibrations_get_suffix)
  end function vibrations_get_suffix


  ! ---------------------------------------------------------
  subroutine vibrations_symmetrize_dyn_matrix(this)
    type(vibrations_t), intent(inout) :: this

    integer :: imat, jmat
    FLOAT :: average

    PUSH_SUB(vibrations_symmetrize_dyn_matrix)

    do imat = 1, this%num_modes
      do jmat = 1, this%num_modes
        average = M_HALF * (this%dyn_matrix(imat, jmat) + this%dyn_matrix(jmat, imat))
        this%dyn_matrix(imat, jmat) = average
        this%dyn_matrix(jmat, imat) = average
      end do
    end do

    POP_SUB(vibrations_symmetrize_dyn_matrix)
  end subroutine vibrations_symmetrize_dyn_matrix

  ! ---------------------------------------------------------
  subroutine vibrations_normalize_dyn_matrix(this, geo)
    type(vibrations_t), intent(inout) :: this
    type(geometry_t),   intent(in)    :: geo

    FLOAT :: factor
    integer :: iatom, idir, jatom, jdir, imat, jmat

    PUSH_SUB(vibrations_normalize_dyn_matrix)

    do iatom = 1, this%natoms
      do idir = 1, this%ndim
        
        imat = vibrations_get_index(this, iatom, idir)

        do jatom = 1, this%natoms
          do jdir = 1, this%ndim
            
            jmat = vibrations_get_index(this, jatom, jdir)

            factor = vibrations_norm_factor(this, geo, iatom, jatom)
            
            this%dyn_matrix(imat, jmat) = this%dyn_matrix(imat, jmat) * factor

          end do
        end do

      end do
    end do

    POP_SUB(vibrations_normalize_dyn_matrix)
  end subroutine vibrations_normalize_dyn_matrix

  ! ---------------------------------------------------------
  FLOAT pure function vibrations_norm_factor(this, geo, iatom, jatom)
    type(vibrations_t), intent(in) :: this
    type(geometry_t),   intent(in) :: geo
    integer,            intent(in) :: iatom
    integer,            intent(in) :: jatom

    vibrations_norm_factor = this%total_mass / &
      (sqrt(species_weight(geo%atom(iatom)%spec)) * sqrt(species_weight(geo%atom(jatom)%spec)))

  end function vibrations_norm_factor

  ! ---------------------------------------------------------
  subroutine vibrations_out_dyn_matrix(this, imat, jmat)
    type(vibrations_t), intent(in) :: this
    integer,            intent(in) :: imat
    integer,            intent(in) :: jmat

    integer :: iatom, idir, jatom, jdir

    if(.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(vibrations_out_dyn_matrix)

    iatom = vibrations_get_atom(this, imat)
    idir  = vibrations_get_dir (this, imat)
    jatom = vibrations_get_atom(this, jmat)
    jdir  = vibrations_get_dir (this, jmat)
    
    write(this%dyn_mat_unit, '(i6, i3, i6, i3, e20.12)') iatom, idir, jatom, jdir, &
      units_from_atomic(unit_invcm**2, this%dyn_matrix(imat, jmat))

    POP_SUB(vibrations_out_dyn_matrix)
  end subroutine vibrations_out_dyn_matrix

  ! ---------------------------------------------------------
  subroutine vibrations_diag_dyn_matrix(this)
    type(vibrations_t), intent(inout) :: this
    
    integer :: imode

    PUSH_SUB(vibrations_diag_dyn_matrix)

    this%normal_mode = M_ZERO
    
    this%normal_mode = this%dyn_matrix
    call lalg_eigensolve(this%num_modes, this%normal_mode, this%freq)

    this%freq(1:this%num_modes) = -this%freq(1:this%num_modes) / this%total_mass

    if(any(this%freq(1:this%num_modes) < -M_EPSILON)) then
      message(1) = "There are imaginary vibrational frequencies (represented as negative)."
      call messages_warning(1)
    endif

    do imode = 1, this%num_modes
      if(this%freq(imode) > M_EPSILON) then
        this%freq(imode) =  sqrt(abs(this%freq(imode)))
      else
        this%freq(imode) = -sqrt(abs(this%freq(imode)))
      endif

      ! make the largest component positive, to specify the phase
      if( maxval(this%normal_mode(:, imode)) - abs(minval(this%normal_mode(:, imode))) < -M_EPSILON) then
        this%normal_mode(:, imode) = -this%normal_mode(:, imode)
      endif
    enddo

    POP_SUB(vibrations_diag_dyn_matrix)
  end subroutine vibrations_diag_dyn_matrix


  ! ---------------------------------------------------------
  integer pure function vibrations_get_index(this, iatom, idim)
    type(vibrations_t), intent(in) :: this
    integer,            intent(in) :: iatom
    integer,            intent(in) :: idim

    vibrations_get_index = (iatom - 1)*this%ndim + idim
  end function vibrations_get_index


  ! ---------------------------------------------------------
  integer pure function vibrations_get_atom(this, index)
    type(vibrations_t), intent(in) :: this
    integer,            intent(in) :: index

    vibrations_get_atom = 1 + (index - 1)/ this%ndim 
  end function vibrations_get_atom


  ! ---------------------------------------------------------
  integer pure function vibrations_get_dir(this, index)
    type(vibrations_t), intent(in) :: this
    integer,            intent(in) :: index

    vibrations_get_dir =  1 + mod(index - 1, this%ndim)
  end function vibrations_get_dir


  ! ---------------------------------------------------------
  subroutine vibrations_output(this)
    type(vibrations_t), intent(in) :: this
    
    integer :: iunit, i, j

    if(.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(vibrations_output)

    ! output frequencies and eigenvectors
    iunit = io_open(VIB_MODES_DIR//'normal_frequencies_'//trim(this%suffix), action='write')
    do i = 1, this%num_modes
      write(iunit, '(i6,f14.5)') i, units_from_atomic(unit_invcm, this%freq(i))
    end do
    call io_close(iunit)

    ! output eigenvectors
    iunit = io_open(VIB_MODES_DIR//'normal_modes_'//trim(this%suffix), action='write')
    do i = 1, this%num_modes
      write(iunit, '(i6)', advance='no') i
      do j = 1, this%num_modes
        write(iunit, '(es14.5)', advance='no') this%normal_mode(j, i)
      end do
      write(iunit, '(1x)')
    end do
    call io_close(iunit)

    POP_SUB(vibrations_output)
  end subroutine vibrations_output

end module vibrations_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
