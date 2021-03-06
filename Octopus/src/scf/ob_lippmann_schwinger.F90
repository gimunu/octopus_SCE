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
!! $Id: scf.F90 4182 2008-05-14 14:02:30Z acastro $

!> This module solves the Schroedinger equation for a system with open
!! boundaries for a prescribed energy.

#include "global.h"

module ob_lippmann_schwinger_m
  use blas_m
  use eigensolver_m
  use global_m
  use grid_m
  use hamiltonian_m
  use hamiltonian_base_m
  use io_m
  use lalg_basic_m
  use loct_m
  use math_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use mpi_lib_m
  use ob_interface_m
  use ob_lead_m
  use profiling_m
  use simul_box_m
  use solvers_m
  use states_m

  implicit none

  private
  public :: &
    lippmann_schwinger

  type p_se_t
    CMPLX, pointer           :: self_energy(:, :, :)
  end type p_se_t

  !> Pointers to communicate with iterative linear solver.
  integer, pointer             :: ist_p, ik_p
  FLOAT, pointer               :: energy_p
  type(p_se_t), pointer        :: lead_p(:)
  type(grid_t), pointer        :: gr_p
  type(hamiltonian_t), pointer :: hm_p
  type(states_t), pointer      :: st_p

contains

  ! ---------------------------------------------------------
  !> Solve the Lippmann-Schwinger equation for the open boundary
  !! system. Use convergence criteria in eigens.
  subroutine lippmann_schwinger(eigens, hm, gr, st)
    type(eigensolver_t),         intent(inout) :: eigens
    type(hamiltonian_t), target, intent(inout) :: hm
    type(grid_t), target,        intent(inout) :: gr
    type(states_t), target,      intent(inout) :: st

    integer                    :: il, iter, np, np_part, lead_np, idim, dim
    integer, target            :: ist, ik
    FLOAT, target              :: energy
    FLOAT                      :: res
    CMPLX, allocatable         :: rhs2(:, :), rhs(:), psi(:)
    type(p_se_t), target       :: lead(2*MAX_DIM)
    logical                    :: conv
#ifdef HAVE_MPI
    integer :: outcount
    FLOAT, allocatable :: ldiff(:), leigenval(:)
#endif

    PUSH_SUB(lippmann_schwinger)

    ! set up pointer for QMR solver
    call mesh_init_mesh_aux(gr%mesh)

    np = gr%mesh%np
    np_part = gr%mesh%np_part
    dim = st%d%dim
    SAFE_ALLOCATE(psi(1:np_part*dim))
    SAFE_ALLOCATE(rhs2(1:np, 1:dim))
    SAFE_ALLOCATE(rhs(1:np*dim))
    do il = 1, NLEADS
      if(gr%intf(il)%reducible) then
        lead_np = gr%intf(il)%np_intf
      else
        lead_np = gr%intf(il)%np_uc
      end if
      SAFE_ALLOCATE(lead(il)%self_energy(1:lead_np, 1:lead_np, 1:dim))
    end do

    eigens%converged = 0
    eigens%matvec    = 0

    ist_p    => ist
    ik_p     => ik
    lead_p   => lead
    gr_p     => gr
    hm_p     => hm
    st_p     => st
    energy_p => energy

    ASSERT(ubound(st%zphi, dim = 1) == np_part)

    ! We have many k-points, so show the progress, but only if not in debug mode since
    ! for every k-point and state the convergence process is shown
    if(.not. in_debug_mode .and. mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, st%d%kpt%nlocal)

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        ! Solve Lippmann-Schwinger equation for this energy.
        energy = st%ob_eigenval(ist, ik)
        st%eigenval(ist, ik) = energy
        do il = 1, NLEADS
          lead(il)%self_energy(:, :, :) = st%ob_lead(il)%self_energy(:, :, :, ist, ik)
        end do

        ! Calculate right hand side e-T-V0-sum(a)[H_ca*g_a*H_ac].
        rhs2(1:np, :) = st%zphi(1:np, :, ist, ik)
        do idim = 1, dim
          psi((idim-1)*np+1:idim*np) = rhs2(1:np, idim)
        end do

        call calc_rhs(rhs2)

        if (associated(hm%ep%A_static)) call calc_rhs(rhs2, transposed = .true.) ! multiply transposed version
        
        ! put in continuous array
        do idim = 1, dim
          rhs((idim-1)*np+1:idim*np) = rhs2(1:np, idim)
        end do

        ! Solve linear system lhs psi = rhs.
        iter = eigens%es_maxiter

        conv = .false.
        if (associated(hm%ep%A_static)) then ! magnetic gs
          call zqmr_sym(dim*np, psi, rhs, lhs_symmetrized, zmf_dotu_aux, zmf_nrm2_aux, &
                        ls_qmr_prec, iter, residue = res, threshold = eigens%tolerance, &
                        converged = conv, showprogress = in_debug_mode)
        else
          call zqmr_sym(dim*np, psi, rhs, lhs, zmf_dotu_aux, zmf_nrm2_aux, ls_qmr_prec, &
                        iter, residue=res, threshold = eigens%tolerance, &
                        converged = conv, showprogress = in_debug_mode)
        end if
        do idim = 1, dim
          call states_set_state(st, gr%mesh, idim, ist, ik, psi((idim-1)*np+1:idim*np))
        end do

        if(in_debug_mode) then ! write info
          write(message(1), '(a,i8,e10.3)') 'Iterations, Residual: ', iter, res
          call messages_info(1)
        end if

        eigens%matvec = eigens%matvec + iter + 1 + 2
        if(conv) eigens%converged = eigens%converged + 1
        eigens%diff(ist, ik) = res
      end do
      if(.not. in_debug_mode .and. mpi_grp_is_root(mpi_world)) then
        call loct_progress_bar(ik-st%d%kpt%start+1, st%d%kpt%nlocal)
      end if
    end do

#ifdef HAVE_MPI
    if(st%d%kpt%parallel) then
      ! every node needs to know all eigenvalues (and diff)
      SAFE_ALLOCATE(ldiff(1:st%d%kpt%nlocal))
      SAFE_ALLOCATE(leigenval(1:st%d%kpt%nlocal))
      do ist = st%st_start, st%st_end
        ldiff(1:st%d%kpt%nlocal) = eigens%diff(ist, st%d%kpt%start:st%d%kpt%end)
        leigenval(1:st%d%kpt%nlocal) = st%eigenval(ist, st%d%kpt%start:st%d%kpt%end)
        call lmpi_gen_allgatherv(st%d%kpt%nlocal, ldiff, outcount, &
                                 eigens%diff(ist, :), st%d%kpt%mpi_grp)
        ASSERT(outcount == st%d%nik)
        call lmpi_gen_allgatherv(st%d%kpt%nlocal, leigenval, outcount, &
                                 st%eigenval(ist, :), st%d%kpt%mpi_grp)
        ASSERT(outcount == st%d%nik)
      end do
      SAFE_DEALLOCATE_A(ldiff)
      SAFE_DEALLOCATE_A(leigenval)
    end if
#endif

    SAFE_DEALLOCATE_A(rhs)
    SAFE_DEALLOCATE_A(rhs2)
    SAFE_DEALLOCATE_A(psi)
    do il = 1, NLEADS
      SAFE_DEALLOCATE_P(lead(il)%self_energy)
    end do

    POP_SUB(lippmann_schwinger)
  end subroutine lippmann_schwinger


  ! ---------------------------------------------------------
  !> The right hand side of the Lippmann-Schwinger equation
  !! e-T-V0-sum(a)[H_ca*g_a*H_ac].
  subroutine calc_rhs(rhs, transposed)
    CMPLX, intent(inout)          :: rhs(:, :)
    logical, optional, intent(in) :: transposed !< needed only for the non-Hermitian part

    integer :: ip, idim, il, np
    integer :: start(1:3), finish(1:3), start_lead(1:3), finish_lead(1:3)
    logical :: transposed_
    CMPLX, allocatable :: tmp(:, :)
    FLOAT, allocatable :: tmp_pot(:)

    PUSH_SUB(calc_rhs)

    np = gr_p%mesh%np

    SAFE_ALLOCATE(tmp(1:gr_p%mesh%np_part, 1:st_p%d%dim))
    SAFE_ALLOCATE(tmp_pot(1:np))

    transposed_ = optional_default(transposed, .false.)

    if(transposed_) then ! the usual conjugate trick for the hermitian part
      tmp(1:np, :) = conjg(rhs(1:np, :))
    else
      tmp(1:np, :) = rhs(1:np, :)
    end if
    ! Calculate right hand side e-T-V0-sum(a)[H_ca*g_a*H_ac].
    rhs(:, :) = M_z0

    call zhamiltonian_apply(hm_p, gr_p%der, tmp, rhs, ist_p, ik_p, terms = TERM_KINETIC)

    ! Apply lead potential. Left and right lead potential are assumed to be equal.
    start(1:3) = gr_p%mesh%idx%nr(1, 1:3) + gr_p%mesh%idx%enlarge(1:3)
    finish(1:3) = gr_p%mesh%idx%nr(2, 1:3) - gr_p%mesh%idx%enlarge(1:3)

    start_lead(1:3) = gr_p%ob_grid%lead(LEFT)%mesh%idx%nr(1, 1:3) + gr_p%ob_grid%lead(LEFT)%mesh%idx%enlarge(1:3)
    finish_lead(1:3) = gr_p%ob_grid%lead(LEFT)%mesh%idx%nr(2, 1:3) - gr_p%ob_grid%lead(LEFT)%mesh%idx%enlarge(1:3)

    do idim = 1, st_p%d%dim
      tmp_pot(:) = M_ZERO
      call dmf_add(gr_p%ob_grid%lead(LEFT)%mesh, start_lead, finish_lead, hm_p%lead(LEFT)%vks(:, idim), &
                    gr_p%mesh, start, finish, tmp_pot(:), TRANS_DIR)
      forall(ip = 1:gr_p%mesh%np)
        rhs(ip, idim) = rhs(ip, idim) + tmp_pot(ip) * tmp(ip, idim)
      end forall
    end do
    ! Add energy.
    forall(ip = 1:np ) rhs(ip, :) = energy_p * tmp(ip, :) - rhs(ip, :)

    if(transposed_) then
      rhs = conjg(rhs)
      tmp = conjg(tmp) ! restore original
    end if

    do il = 1, NLEADS
      do idim = 1, st_p%d%dim
        if(transposed_) then
          call interface_apply_op(gr_p%intf(il), -M_z1, transpose(lead_p(il)%self_energy(:, :, idim)), &
                                tmp(:, idim), rhs(:, idim))
        else
          call interface_apply_op(gr_p%intf(il), -M_z1, lead_p(il)%self_energy(:, :, idim), &
                                tmp(:, idim), rhs(:, idim))
        end if
      end do
    end do

    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(tmp_pot)
    POP_SUB(calc_rhs)
  end subroutine calc_rhs


  ! ---------------------------------------------------------
  !> The left hand side of the Lippmann-Schwinger equation
  !! e-H-sum(a)[H_ca*g_a*H_ac].
  !! Used by the iterative linear solver.
  subroutine lhs(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    CMPLX, allocatable :: tmp_x(:, :)
    CMPLX, allocatable :: tmp_y(:, :)
    integer            :: np, np_part, idim, il, dim

    ! no push_sub, called too frequently

    np_part = gr_p%mesh%np_part
    np      = gr_p%mesh%np
    dim     = st_p%d%dim

    SAFE_ALLOCATE(tmp_x(1:np_part, 1:dim))
    SAFE_ALLOCATE(tmp_y(1:np, 1:dim))

    do idim = 1, dim
      tmp_x(1:np, idim) = x((idim - 1)*np + 1:idim*np)
    end do
    call zhamiltonian_apply(hm_p, gr_p%der, tmp_x, tmp_y, ist_p, ik_p)

    ! y <- e x - tmp_y
    do idim = 1, dim
      tmp_y(1:np, idim) = energy_p * x((idim - 1)*np + 1:idim*np) - tmp_y(1:np, idim)
    end do

    do il = 1, NLEADS
      do idim = 1, dim
        call interface_apply_op(gr_p%intf(il), -M_z1, &
          lead_p(il)%self_energy(:, :, idim), tmp_x(:, idim), tmp_y(:, idim))
      end do
    end do

    do idim = 1, dim
      y((idim - 1)*np + 1:idim*np) = tmp_y(1:np, idim)
    end do

    SAFE_DEALLOCATE_A(tmp_x)
    SAFE_DEALLOCATE_A(tmp_y)

  end subroutine lhs


  ! ---------------------------------------------------------
  !> The left hand side of the Lippmann-Schwinger equation
  !! (e-H-sum(a)[H_ca*g_a*H_ac])^T.
  !! Used by the iterative linear solver.
  subroutine lhs_t(y)
    CMPLX, intent(inout) :: y(:)

    CMPLX, allocatable :: tmp_x(:, :)
    CMPLX, allocatable :: tmp_y(:, :)
    integer            :: np, np_part, idim, il, dim

    ! no push_sub, called too frequently

    np      = gr_p%mesh%np
    np_part = gr_p%mesh%np_part
    dim     = st_p%d%dim

    SAFE_ALLOCATE(tmp_x(1:np_part, 1:dim))
    SAFE_ALLOCATE(tmp_y(1:np, 1:dim))

    do idim = 1, dim
      tmp_x(1:np, idim) = conjg(y((idim-1)*np+1:idim*np))
    end do
    call zhamiltonian_apply(hm_p, gr_p%der, tmp_x, tmp_y, ist_p, ik_p)

    ! y <- e x - tmp_y
    do idim = 1, dim
      tmp_y(1:np, idim) = energy_p * tmp_x(1:np, idim) - tmp_y(1:np, idim)
    end do
    tmp_y = conjg(tmp_y)
    tmp_x = conjg(tmp_x) ! restore for the non-Hermitian part

    do il = 1, NLEADS
      do idim = 1, dim
        call interface_apply_op(gr_p%intf(il), -M_z1, transpose(lead_p(il)%self_energy(:, :, idim)), &
                                tmp_x(:, idim), tmp_y(:, idim))
      end do
    end do

    do idim = 1, dim
      y((idim-1)*np+1:idim*np) = tmp_y(1:np, idim)
    end do

    SAFE_DEALLOCATE_A(tmp_x)
    SAFE_DEALLOCATE_A(tmp_y)

  end subroutine lhs_t


  ! ---------------------------------------------------------
  !> The left hand side of the Lippmann-Schwinger equation
  !! (e-H-sum(a)[H_ca*g_a*H_ac])^T*(e-H-sum(a)[H_ca*g_a*H_ac]).
  !! Used by the iterative linear solver.
  subroutine lhs_symmetrized(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    ! no push_sub, called too frequently

    call lhs(x, y)
    call lhs_t(y)

  end subroutine lhs_symmetrized


  ! ---------------------------------------------------------
  !> Identity preconditioner. Since preconditioning with the inverse of
  !! the diagonal did not improve the convergence we put identity here
  !! until we have something better.
  subroutine ls_qmr_prec(x, y)
    CMPLX, intent(in)  :: x(:)
    CMPLX, intent(out) :: y(:)

    ! no push_sub, called too frequently

    y(:) = x(:)

  end subroutine ls_qmr_prec

end module ob_lippmann_schwinger_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
