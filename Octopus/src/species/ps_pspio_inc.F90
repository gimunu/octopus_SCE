!! Copyright (C) 2012 M. Oliveira
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
!! $Id: ps_pspio_inc.F90 10978 2013-07-11 15:28:46Z micael $

  ! ---------------------------------------------------------
  subroutine ps_pspio_init(ps, z, lmax, lloc, ispin, filename)
    type(ps_t),        intent(out)   :: ps
    integer,           intent(inout) :: lmax
    integer,           intent(in)    :: lloc, ispin
    FLOAT,             intent(in)    :: z
    character(len=*),  intent(in)    :: filename

#if HAVE_PSPIO

    logical :: found, has_kb
    integer :: idir
    character(len=3) :: psp_dir(3) = (/"PSF", "FHI", "UPF"/)
    character(len=256) :: filename2
    type(pspio_f90_pspdata_t)   :: pspdata

    PUSH_SUB(ps_pspio_init)

    call messages_experimental("Reading pseudopotential file using PSPIO library")

    ! Find out the file
    filename2 = trim(filename)
    inquire(file=filename2, exist=found)

    if(.not. found) then
      do idir = 1, size(psp_dir)
        filename2 = trim(conf%share) // "/PP/" // trim(psp_dir(idir)) // "/" // trim(filename)
        inquire(file=filename2, exist=found)

        if(found) exit

        if (idir == size(psp_dir)) then
          message(1) = "Pseudopotential file '" // trim(filename) // " not found"
          call messages_fatal(1)
        end if
      end do
    end if

    message(1) = "Reading pseudopotential from file:"
    write(message(2), '(6x,3a)') "'", trim(filename2), "'"
    call messages_info(2)

    ! Init pspio data structure and parse file
    call check_error(pspio_f90_pspdata_init(pspdata))
    call check_error(pspio_f90_pspdata_read(pspdata, PSPIO_UNKNOWN, filename2))

    ! General info
    ps%ispin = ispin
    ps%z = z
    ps%conf%z = nint(z)
    call ps_pspio_read_info(ps, pspdata)
    lmax = ps%l_max
    write(message(1), '(a,i2,a)') "Info: l = ", ps%l_max, " is maximum angular momentum considered."
    call messages_info(1)

    ! Mesh
    call ps_pspio_read_mesh(ps, pspdata)

    ! XC
    call ps_pspio_read_xc(ps, pspdata)

    ! We will first try to read the KB projectors
    call ps_pspio_read_kb_projectors(ps, pspdata, has_kb)

    ! States
    call ps_pspio_read_states(ps, pspdata)

    ! If we do not have KB projectors, then we read the pseudopotentials
    if (.not. has_kb) then
      call ps_pspio_read_potentials(ps, pspdata)
    end if

    !No variable description, as it is already in ps.F90
    call parse_float(datasets_check('SpeciesProjectorSphereThreshold'), &
      CNST(0.001), ps%projectors_sphere_threshold)
    if(ps%projectors_sphere_threshold <= M_ZERO) call input_error('SpeciesProjectorSphereThreshold')
    ps%has_long_range = .true.
    ps%is_separated = .false.

    !Free memory
    call pspio_f90_pspdata_free(pspdata)

#else
    message(1) = 'PSPIO selected for pseudopotential parsing, but the code was compiled witout PSPIO support.'
    call messages_fatal(1)
#endif

    POP_SUB(ps_pspio_init)
  end subroutine ps_pspio_init


#if HAVE_PSPIO

  ! ---------------------------------------------------------
  subroutine ps_pspio_read_info(ps, pspdata)
    type(ps_t),                intent(inout) :: ps
    type(pspio_f90_pspdata_t), intent(in)    :: pspdata

    PUSH_SUB(ps_pspio_read_info)

    call pspio_f90_pspdata_get_symbol(pspdata, ps%conf%symbol)
    call pspio_f90_pspdata_get_l_max(pspdata, ps%l_max)
    call pspio_f90_pspdata_get_zvalence(pspdata, ps%z_val)

    POP_SUB(ps_pspio_read_info)
  end subroutine ps_pspio_read_info

  ! ---------------------------------------------------------
  subroutine ps_pspio_read_mesh(ps, pspdata)
    type(ps_t),                intent(inout) :: ps
    type(pspio_f90_pspdata_t), intent(in)    :: pspdata

    integer :: ip
    type(pspio_f90_mesh_t) :: mesh
    FLOAT, allocatable :: r_tmp(:)

    PUSH_SUB(ps_pspio_read_mesh)

    call pspio_f90_pspdata_get_mesh(pspdata, mesh)
    call pspio_f90_mesh_get_np(mesh, ps%g%nrval)

    SAFE_ALLOCATE(r_tmp(1:ps%g%nrval))
    call pspio_f90_mesh_get_r(mesh, r_tmp(1))
    if(any(abs(r_tmp(2:ps%g%nrval)) < M_EPSILON)) then
      ! only the first point is allowed to be zero
      message(1) = "Illegal zero values in PSPIO radial grid"
      call messages_fatal(1)
    endif
    if (r_tmp(1) == M_ZERO) then
      ip = 1
    else
      ip = 2
      ps%g%nrval = ps%g%nrval + 1
    end if

    nullify(ps%g%drdi, ps%g%s)
    SAFE_ALLOCATE(ps%g%rofi(1:ps%g%nrval))
    SAFE_ALLOCATE(ps%g%r2ofi(1:ps%g%nrval))
    ps%g%rofi(1) = M_ZERO
    ps%g%rofi(ip:ps%g%nrval) = r_tmp(1:ps%g%nrval-ip+1)
    ps%g%r2ofi = ps%g%rofi**2
    SAFE_DEALLOCATE_A(r_tmp)

    POP_SUB(ps_pspio_read_mesh)
  end subroutine ps_pspio_read_mesh

  ! ---------------------------------------------------------
  subroutine ps_pspio_read_states(ps, pspdata)
    type(ps_t),                intent(inout) :: ps
    type(pspio_f90_pspdata_t), intent(in)    :: pspdata

    integer :: is, ist, l
    FLOAT :: x, j
    type(pspio_f90_state_t) :: state
    FLOAT, allocatable :: wfs(:)

    PUSH_SUB(ps_pspio_read_states)

    call pspio_f90_pspdata_get_n_states(pspdata, ps%conf%p)
    SAFE_ALLOCATE(ps%ur   (1:ps%conf%p, 1:ps%ispin))
    SAFE_ALLOCATE(ps%ur_sq(1:ps%conf%p, 1:ps%ispin))
    do ist = 1, ps%conf%p
      call pspio_f90_pspdata_get_state(pspdata, ist, state)

      !Quantum numbers
      call pspio_f90_state_get_qn(state, ps%conf%n(ist), ps%conf%l(ist), j)

      !Occupations
      call pspio_f90_state_get_occ(state, ps%conf%occ(ist, 1))
      if(ps%ispin == 2) then
        ! Spin-dependent pseudopotentials are not supported, so we need to fix the occupations
        ! if we want to have a spin-dependent atomic density.      
        x = ps%conf%occ(l, 1)
        ps%conf%occ(ist, 1) = min(x, real(2*ps%conf%l(ist)+1, REAL_PRECISION))
        ps%conf%occ(ist, 2) = x - ps%conf%occ(l, 1)
      end if

      !Wavefunctions
      SAFE_ALLOCATE(wfs(ps%g%nrval))
      call pspio_f90_state_wf_eval(state, ps%g%nrval, ps%g%rofi, wfs)
      do is = 1, ps%ispin
        call spline_fit(ps%g%nrval, ps%g%rofi, wfs, ps%ur(ist, is))
      end do
      do is = 1, ps%ispin
        call spline_fit(ps%g%nrval, ps%g%r2ofi, wfs, ps%ur_sq(ist, is))
      end do
      SAFE_DEALLOCATE_A(wfs)
    end do

    POP_SUB(ps_pspio_read_states)
  end subroutine ps_pspio_read_states

  ! ---------------------------------------------------------
  subroutine ps_pspio_read_kb_projectors(ps, pspdata, has_kb)
    type(ps_t),                intent(inout) :: ps
    type(pspio_f90_pspdata_t), intent(in)    :: pspdata
    logical,                   intent(out)   :: has_kb

    integer :: ir, n_kbproj, wave_eq, ikb, ikbc, l
    FLOAT :: j
    FLOAT, parameter :: threshold = CNST(0.5e-7)
    type(pspio_f90_potential_t) :: vlocal
    type(pspio_f90_projector_t) :: kb_projector
    logical, allocatable :: was_init(:,:)
    FLOAT, allocatable :: v_local(:), proj(:)

    PUSH_SUB(ps_pspio_read_kb_projectors)

    call pspio_f90_pspdata_get_n_kbproj(pspdata, n_kbproj)

    has_kb = .true.
    if (n_kbproj == 0) then
      has_kb = .false.
      POP_SUB(ps_pspio_read_kb_projectors)
      return
    end if

    ! Local potential
    call pspio_f90_pspdata_get_l_local(pspdata, ps%l_loc)
    call pspio_f90_pspdata_get_vlocal(pspdata, vlocal)
    SAFE_ALLOCATE(v_local(ps%g%nrval))
    call pspio_f90_potential_eval(vlocal, ps%g%nrval, ps%g%rofi, v_local)
    do ir = ps%g%nrval-1, 2, -1
      if(abs(v_local(ir)*ps%g%rofi(ir) + ps%z_val) > threshold) exit
    end do
    ps%rc_max = ps%g%rofi(ir + 1)
    call spline_init(ps%vl)
    call spline_fit(ps%g%nrval, ps%g%rofi, v_local, ps%vl)
    SAFE_DEALLOCATE_A(v_local)

    ! KB projectors 
    call pspio_f90_pspdata_get_wave_eq(pspdata, wave_eq)
    if (wave_eq == PSPIO_DIRAC) then
      ps%kbc = 2
    else
      ps%kbc = 1
    end if

    SAFE_ALLOCATE(ps%kb (0:ps%l_max, 1:ps%kbc))
    SAFE_ALLOCATE(ps%dkb(0:ps%l_max, 1:ps%kbc))
    SAFE_ALLOCATE(ps%h  (0:ps%l_max, 1:ps%kbc, 1:ps%kbc))
    SAFE_ALLOCATE(ps%k  (0:ps%l_max, 1:ps%kbc, 1:ps%kbc))
    SAFE_ALLOCATE(was_init(0:ps%l_max, 1:ps%kbc))
    SAFE_ALLOCATE(proj(ps%g%nrval))
    call spline_init(ps%kb)
    call spline_init(ps%dkb)

    was_init = .false.
    ps%h = M_ZERO
    ps%k = M_ZERO
    do ikb = 1, n_kbproj
      call pspio_f90_pspdata_get_kb_projector(pspdata, ikb, kb_projector)

      call pspio_f90_projector_get_l(kb_projector, l)
      call pspio_f90_projector_get_j(kb_projector, j)
      ikbc = 1
      if (j == real(l, REAL_PRECISION) - M_HALF) ikbc = 2

      call pspio_f90_projector_get_energy(kb_projector, ps%h(l, ikbc, ikbc))

      call pspio_f90_projector_eval(kb_projector, ps%g%nrval, ps%g%rofi, proj)
      do ir = ps%g%nrval-1, 2, -1
        if(abs(proj(ir)) > threshold) exit
      end do
      ps%rc_max = max(ps%g%rofi(ir + 1), ps%rc_max)
      proj(ir+1:ps%g%nrval) = M_ZERO
      call spline_fit(ps%g%nrval, ps%g%rofi, proj, ps%kb(l, ikbc))
      was_init(l, ikbc) = .true.
    end do

    !Make sure all projector splines were initialized
    proj = M_ZERO
    do l = 0, ps%l_max
      do ikbc = 1, ps%kbc
        if (.not. was_init(l, ikbc)) then
          call spline_fit(ps%g%nrval, ps%g%rofi, proj, ps%kb(l, ikbc))
        end if
      end do
    end do
    SAFE_DEALLOCATE_A(proj)
    SAFE_DEALLOCATE_A(was_init)

    ! Increasing radius a little, just in case.
    ! I have hard-coded a larger increase of the cutoff for the filtering.
    ps%rc_max = ps%rc_max*CNST(1.5)

    POP_SUB(ps_pspio_read_kb_projectors)
  end subroutine ps_pspio_read_kb_projectors

  ! ---------------------------------------------------------
  subroutine ps_pspio_read_potentials(ps, pspdata)
    type(ps_t),                intent(inout) :: ps
    type(pspio_f90_pspdata_t), intent(in)    :: pspdata

    PUSH_SUB(ps_pspio_read_potentials)

    message(1) = "Not yet implemented"
    call messages_fatal(1)

    POP_SUB(ps_pspio_read_potentials)
  end subroutine ps_pspio_read_potentials

  ! ---------------------------------------------------------
  subroutine ps_pspio_read_xc(ps, pspdata)
    type(ps_t),                intent(inout) :: ps
    type(pspio_f90_pspdata_t), intent(in)    :: pspdata

    logical :: has_nlcc
    integer :: ir, nrc
    FLOAT, parameter :: threshold = CNST(0.5e-7)
    type(pspio_f90_xc_t) :: xc
    FLOAT, allocatable :: rho(:)

    PUSH_SUB(ps_pspio_read_xc)

    !Non-linear core-corrections
    call pspio_f90_pspdata_get_xc(pspdata, xc)

    call spline_init(ps%core)
    call pspio_f90_xc_has_nlcc(xc, has_nlcc)
    if (has_nlcc) then
      ps%icore=''

      ! get core density
      SAFE_ALLOCATE(rho(ps%g%nrval))
      call pspio_f90_xc_nlcc_eval(xc, ps%g%nrval, ps%g%rofi, rho)

      ! find cutoff radius
      do ir = ps%g%nrval-1, 1, -1
        if(rho(ir) > eps) then
          nrc = ir + 1
          exit
        end if
      end do
      rho(nrc:ps%g%nrval) = M_ZERO

      call spline_fit(ps%g%nrval, ps%g%rofi, rho, ps%core)

      SAFE_DEALLOCATE_A(rho)
    else
      ps%icore = 'nc'
    end if

    POP_SUB(ps_pspio_read_xc)
  end subroutine ps_pspio_read_xc

  ! ---------------------------------------------------------
  subroutine check_error(ierr)
    integer, intent(in) :: ierr

    integer :: ierr2

    if (ierr /= PSPIO_SUCCESS) then
      ierr2 = pspio_f90_error_flush()
      message(1) = "PSPIO error"
      call messages_fatal(1)
    end if

  end subroutine check_error

#endif

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
