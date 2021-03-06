!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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
!! $Id: poisson.F90 2660 2007-01-23 15:11:54Z lorenzen $

#include "global.h"

module poisson_fft_m
  use cube_function_m
  use cube_m
  use datasets_m
  use fft_m
  use fourier_space_m
  use geometry_m
  use global_m
  use lalg_basic_m
  use loct_m
  use loct_math_m
  use math_m
  use mesh_cube_parallel_map_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mpi_m
  use par_vec_m
  use parser_m
  use poisson_cutoff_m
  use profiling_m
  use simul_box_m
  use splines_m
  use unit_m
  use unit_system_m

  implicit none

  private
  public ::                  &
    poisson_fft_t,           &
    poisson_fft_init,        &
    poisson_fft_end,         &
    poisson_fft_solve

  integer, public, parameter ::                &
       POISSON_FFT_KERNEL_NONE      = -1,      &
       POISSON_FFT_KERNEL_SPH       =  0,      &
       POISSON_FFT_KERNEL_CYL       =  1,      &
       POISSON_FFT_KERNEL_PLA       =  2,      &
       POISSON_FFT_KERNEL_NOCUT     =  3,      &
       POISSON_FFT_KERNEL_CORRECTED =  4

  type poisson_fft_t
    type(fourier_space_op_t) :: coulb  !< object for Fourier space operations
    integer                  :: kernel !< choice of kernel, one of options above
    FLOAT                    :: qq(MAX_DIM) !< q-point for exchange in periodic system
  end type poisson_fft_t
contains

  subroutine poisson_fft_init(this, mesh, cube, kernel, soft_coulb_param, qq)
    type(poisson_fft_t), intent(out)   :: this
    type(mesh_t),        intent(in)    :: mesh
    type(cube_t),        intent(inout) :: cube
    integer,             intent(in)    :: kernel
    FLOAT, optional,     intent(in)    :: soft_coulb_param
    FLOAT, optional,     intent(in)    :: qq(:) !< (1:mesh%sb%periodic_dim)

    PUSH_SUB(poisson_fft_init)

    this%kernel = kernel
    this%qq = M_ZERO

    if(present(qq) .and. simul_box_is_periodic(mesh%sb)) then
      ASSERT(ubound(qq, 1) >= mesh%sb%periodic_dim)
      this%qq(1:mesh%sb%periodic_dim) = qq(1:mesh%sb%periodic_dim)
    endif
    
    select case(mesh%sb%dim)
    case(1)
      ASSERT(present(soft_coulb_param))
      select case(kernel)
      case(POISSON_FFT_KERNEL_SPH)
        call poisson_fft_build_1d_0d(this, mesh, cube, soft_coulb_param)
      case(POISSON_FFT_KERNEL_NOCUT)
        call poisson_fft_build_1d_1d(this, mesh, cube, soft_coulb_param)
      case default
        message(1) = "Invalid Poisson FFT kernel for 1D."
        call messages_fatal(1)
      end select

    case(2)
      select case(kernel)
      case(POISSON_FFT_KERNEL_SPH)
        call poisson_fft_build_2d_0d(this, mesh, cube)
      case(POISSON_FFT_KERNEL_CYL)
        call poisson_fft_build_2d_1d(this, mesh, cube)
      case(POISSON_FFT_KERNEL_NOCUT)
        call poisson_fft_build_2d_2d(this, mesh, cube)
      case default
        message(1) = "Invalid Poisson FFT kernel for 2D."
        call messages_fatal(1)
      end select
      
    case(3)
      select case(kernel)
      case(POISSON_FFT_KERNEL_SPH, POISSON_FFT_KERNEL_CORRECTED)
        call poisson_fft_build_3d_0d(this, mesh, cube, kernel)
        
      case(POISSON_FFT_KERNEL_CYL)
        call poisson_fft_build_3d_1d(this, mesh, cube)
        
      case(POISSON_FFT_KERNEL_PLA)
        call poisson_fft_build_3d_2d(this, mesh, cube)
        
      case(POISSON_FFT_KERNEL_NOCUT)
        call poisson_fft_build_3d_3d(this, mesh, cube)

      case default
        message(1) = "Invalid Poisson FFT kernel for 3D."
        call messages_fatal(1)        
      end select
    end select

    POP_SUB(poisson_fft_init)
  end subroutine poisson_fft_init

  !-----------------------------------------------------------------

  subroutine get_cutoff(default_r_c, r_c)
    FLOAT, intent(in)  :: default_r_c
    FLOAT, intent(out) :: r_c

    PUSH_SUB(get_cutoff)

    call parse_float(datasets_check('PoissonCutoffRadius'), default_r_c, r_c, units_inp%length)
    
    call messages_write('Info: Poisson Cutoff Radius     =')
    call messages_write(r_c, units = units_out%length, fmt = '(f6.1)')
    call messages_info()
    
    if ( r_c > default_r_c + M_EPSILON) then
      call messages_write('Poisson cutoff radius is larger than cell size.', new_line = .true.)
      call messages_write('You can see electrons in neighboring cell(s).')
      call messages_warning()
    end if

    POP_SUB(get_cutoff)
  end subroutine get_cutoff

  !-----------------------------------------------------------------
  subroutine poisson_fft_gg_transform(gg_in, sb, qq, gg, modg2)
    FLOAT,             intent(in)    :: gg_in(:)
    type(simul_box_t), intent(in)    :: sb
    FLOAT,             intent(in)    :: qq(:)
    FLOAT,             intent(inout) :: gg(:)
    FLOAT,             intent(out)   :: modg2

    integer :: idir

    ! no PUSH_SUB, called too frequently

    gg(1:3) = matmul(gg_in(1:3), sb%klattice_primitive(1:3,1:3))
    do idir = 1, 3
      gg(idir) = gg(idir) / lalg_nrm2(3, sb%klattice_primitive(1:3, idir))
    end do
    gg(1:sb%periodic_dim) = gg(1:sb%periodic_dim) + qq(1:sb%periodic_dim)

    modg2 = sum(gg(1:3)**2)

  end subroutine poisson_fft_gg_transform

  !-----------------------------------------------------------------
  subroutine poisson_fft_build_3d_3d(this, mesh, cube)
    type(poisson_fft_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(cube_t),        intent(inout) :: cube

    integer :: ix, iy, iz, ixx(3), db(3)
    FLOAT :: temp(3), modg2
    FLOAT :: gg(3)
    FLOAT, allocatable :: fft_Coulb_FS(:,:,:)
    
    PUSH_SUB(poisson_fft_build_3d_3d)

    db(1:3) = cube%rs_n_global(1:3)

    ! store the Fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:cube%fs_n_global(1), 1:cube%fs_n_global(2), 1:cube%fs_n_global(3)))
    fft_Coulb_FS = M_ZERO
    temp(1:3) = M_TWO*M_PI/(db(1:3)*mesh%spacing(1:3))

    do ix = 1, cube%fs_n_global(1)
      ixx(1) = pad_feq(ix, db(1), .true.)
      do iy = 1, cube%fs_n_global(2)
        ixx(2) = pad_feq(iy, db(2), .true.)
        do iz = 1, cube%fs_n_global(3)
          ixx(3) = pad_feq(iz, db(3), .true.)

          call poisson_fft_gg_transform(temp * ixx, mesh%sb, this%qq, gg, modg2)

          if(abs(modg2) > M_EPSILON) then
            fft_Coulb_FS(ix, iy, iz) = M_ONE/modg2
          else
            fft_Coulb_FS(ix, iy, iz) = M_ZERO
          end if
        end do
      end do

    end do

    forall(iz=1:cube%fs_n_global(3), iy=1:cube%fs_n_global(2), ix=1:cube%fs_n_global(1))
      fft_Coulb_FS(ix, iy, iz) = M_FOUR*M_PI*fft_Coulb_FS(ix, iy, iz)
    end forall

    call dfourier_space_op_init(this%coulb, cube, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    POP_SUB(poisson_fft_build_3d_3d)
  end subroutine poisson_fft_build_3d_3d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  !> C. A. Rozzi et al., Phys. Rev. B 73, 205119 (2006), Table I
  subroutine poisson_fft_build_3d_2d(this, mesh, cube)
    type(poisson_fft_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(cube_t),        intent(inout) :: cube

    integer :: ix, iy, iz, ixx(3), db(3)
    FLOAT :: temp(3), modg2
    FLOAT :: gpar, gz, r_c, gg(3), default_r_c
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_3d_2d)

    db(1:3) = cube%rs_n_global(1:3)

    !%Variable PoissonCutoffRadius
    !%Type float
    !%Section Hamiltonian::Poisson
    !%Description
    !% When <tt>PoissonSolver = fft</tt> and <tt>PoissonFFTKernel</tt> is neither <tt>multipole_corrections</tt>
    !% nor <tt>fft_nocut</tt>,
    !% this variable controls the distance after which the electron-electron interaction goes to zero.
    !% A warning will be written if the value is too large and will cause spurious interactions between images.
    !% The default is half of the FFT box max dimension in a finite direction.
    !%End

    default_r_c = db(1)*mesh%spacing(1)/M_TWO
    call get_cutoff(default_r_c, r_c)

    ! store the fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:cube%fs_n_global(1), 1:cube%fs_n_global(2), 1:cube%fs_n_global(3)))
    fft_Coulb_FS = M_ZERO

    temp(1:3) = M_TWO*M_PI/(db(1:3)*mesh%spacing(1:3))

    do ix = 1, cube%fs_n_global(1)
      ixx(1) = pad_feq(ix, db(1), .true.)
      do iy = 1, cube%fs_n_global(2)
        ixx(2) = pad_feq(iy, db(2), .true.)
        do iz = 1, cube%fs_n_global(3)
          ixx(3) = pad_feq(iz, db(3), .true.)

          call poisson_fft_gg_transform(temp * ixx, mesh%sb, this%qq, gg, modg2)

          if(abs(modg2) > M_EPSILON) then
            gz = abs(gg(3))
            gpar = hypot(gg(1), gg(2))
            ! note: if gpar = 0, then modg2 = gz**2
            fft_Coulb_FS(ix, iy, iz) = poisson_cutoff_3D_2D(gpar,gz,r_c)/modg2
          else
            fft_Coulb_FS(ix, iy, iz) = -M_HALF*r_c**2
          end if
        end do
      end do

    end do

    forall(iz=1:cube%fs_n_global(3), iy=1:cube%fs_n_global(2), ix=1:cube%fs_n_global(1))
      fft_Coulb_FS(ix, iy, iz) = M_FOUR*M_PI*fft_Coulb_FS(ix, iy, iz)
    end forall

    call dfourier_space_op_init(this%coulb, cube, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    POP_SUB(poisson_fft_build_3d_2d)
  end subroutine poisson_fft_build_3d_2d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  !> C. A. Rozzi et al., Phys. Rev. B 73, 205119 (2006), Table I
  subroutine poisson_fft_build_3d_1d(this, mesh, cube)
    type(poisson_fft_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(cube_t),        intent(inout) :: cube

    type(spline_t)     :: cylinder_cutoff_f
    FLOAT, allocatable :: x(:), y(:)
    integer :: ix, iy, iz, ixx(3), db(3), k, ngp
    FLOAT :: temp(3), modg2, xmax
    FLOAT :: gperp, gx, gy, gz, r_c, gg(3), default_r_c
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_3d_1d)

    db(1:3) = cube%rs_n_global(1:3)

    default_r_c = maxval(db(2:3)*mesh%spacing(2:3)/M_TWO)
    call get_cutoff(default_r_c, r_c)

    ! store the fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:cube%fs_n_global(1), 1:cube%fs_n_global(2), 1:cube%fs_n_global(3)))
    fft_Coulb_FS = M_ZERO

    temp(1:3) = M_TWO*M_PI/(db(1:3)*mesh%spacing(1:3))

    if( mesh%sb%periodic_dim == 0 ) then
      ngp = 8*db(2)
      SAFE_ALLOCATE(x(1:ngp))
      SAFE_ALLOCATE(y(1:ngp))
    end if


    do ix = 1, cube%fs_n_global(1)
      ixx(1) = pad_feq(ix, db(1), .true.)
      gx = temp(1)*ixx(1)

      if( mesh%sb%periodic_dim == 0 ) then
        call spline_init(cylinder_cutoff_f)
        xmax = sqrt((temp(2)*db(2)/2)**2 + (temp(3)*db(3)/2)**2)
        do k = 1, ngp
          x(k) = (k-1)*(xmax/(ngp-1))
          y(k) = poisson_cutoff_3D_1D_finite(gx, x(k), M_TWO*mesh%sb%xsize, M_TWO*mesh%sb%rsize)
        end do
        call spline_fit(ngp, x, y, cylinder_cutoff_f)
      end if

      do iy = 1, cube%fs_n_global(2)
        ixx(2) = pad_feq(iy, db(2), .true.)
        do iz = 1, db(3)
          ixx(3) = pad_feq(iz, db(3), .true.)

          call poisson_fft_gg_transform(temp * ixx, mesh%sb, this%qq, gg, modg2)

          if(abs(modg2) > M_EPSILON) then
            gperp = hypot(gg(2), gg(3))
            if (mesh%sb%periodic_dim==1) then
              fft_Coulb_FS(ix, iy, iz) = poisson_cutoff_3D_1D(abs(gx), gperp, r_c)/modg2
            else if (mesh%sb%periodic_dim==0) then
              gy = gg(2)
              gz = gg(3)
              if ((gz >= M_ZERO) .and. (gy >= M_ZERO)) then
                fft_Coulb_FS(ix, iy, iz) = spline_eval(cylinder_cutoff_f, gperp)
              end if
              if ((gz >= M_ZERO) .and. (gy < M_ZERO)) then
                fft_Coulb_FS(ix, iy, iz) = fft_Coulb_FS(ix, -ixx(2) + 1, iz)
              end if
              if ((gz < M_ZERO) .and. (gy >= M_ZERO)) then
                fft_Coulb_FS(ix, iy, iz) = fft_Coulb_FS(ix, iy, -ixx(3) + 1)
              end if
              if ((gz < M_ZERO) .and. (gy < M_ZERO) ) then
                fft_Coulb_FS(ix, iy, iz) = fft_Coulb_FS(ix, -ixx(2) + 1, -ixx(3) + 1)
              end if
            end if

          else
            if (mesh%sb%periodic_dim == 1) then
              fft_Coulb_FS(ix, iy, iz) = -(M_HALF*log(r_c) - M_FOURTH)*r_c**2
            else if (mesh%sb%periodic_dim == 0) then
              fft_Coulb_FS(ix, iy, iz) = poisson_cutoff_3D_1D_finite(M_ZERO, M_ZERO, &
                M_TWO*mesh%sb%xsize, M_TWO*mesh%sb%rsize)
            end if

          end if
        end do
      end do
 
      if( mesh%sb%periodic_dim == 0 ) call spline_end(cylinder_cutoff_f)
    end do

    forall(iz=1:cube%fs_n_global(3), iy=1:cube%fs_n_global(2), ix=1:cube%fs_n_global(1))
      fft_Coulb_FS(ix, iy, iz) = M_FOUR*M_PI*fft_Coulb_FS(ix, iy, iz)
    end forall

    call dfourier_space_op_init(this%coulb, cube, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)
    POP_SUB(poisson_fft_build_3d_1d)
  end subroutine poisson_fft_build_3d_1d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  !> C. A. Rozzi et al., Phys. Rev. B 73, 205119 (2006), Table I
  subroutine poisson_fft_build_3d_0d(this, mesh, cube, kernel)
    type(poisson_fft_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(cube_t),        intent(inout) :: cube
    integer,             intent(in)    :: kernel

    integer :: ix, iy, iz, ixx(3), db(3), lx, ly, lz, n1, n2, n3
    FLOAT :: temp(3), modg2
    FLOAT :: r_c, gg(3), default_r_c
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_3d_0d)

    db(1:3) = cube%rs_n_global(1:3)

    if (kernel /= POISSON_FFT_KERNEL_CORRECTED) then
      default_r_c = maxval(db(1:3)*mesh%spacing(1:3)/M_TWO)
      call get_cutoff(default_r_c, r_c)
    end if

    n1 = max(1, cube%fs_n(1))
    n2 = max(1, cube%fs_n(2))
    n3 = max(1, cube%fs_n(3))

    ! store the fourier transform of the Coulomb interaction
    ! store only the relevant part if PFFT is used
    SAFE_ALLOCATE(fft_Coulb_FS(1:n1,1:n2,1:n3))
    fft_Coulb_FS = M_ZERO

    temp(1:3) = M_TWO*M_PI/(db(1:3)*mesh%spacing(1:3))
    do lx = 1, n1
      ix = cube%fs_istart(1) + lx - 1
      ixx(1) = pad_feq(ix, db(1), .true.)
      do ly = 1, n2
        iy = cube%fs_istart(2) + ly - 1
        ixx(2) = pad_feq(iy, db(2), .true.)
        do lz = 1, n3
          iz = cube%fs_istart(3) + lz - 1
          ixx(3) = pad_feq(iz, db(3), .true.)
            
          call poisson_fft_gg_transform(temp * ixx, mesh%sb, this%qq, gg, modg2)

          if(abs(modg2) > M_EPSILON) then
            select case(kernel)
            case(POISSON_FFT_KERNEL_SPH)
              fft_Coulb_FS(lx, ly, lz) = poisson_cutoff_3D_0D(sqrt(modg2),r_c)/modg2
            case(POISSON_FFT_KERNEL_CORRECTED)
              fft_Coulb_FS(lx, ly, lz) = M_ONE/modg2
            end select
          else
            select case(kernel)
            case(POISSON_FFT_KERNEL_SPH)
              fft_Coulb_FS(lx, ly, lz) = r_c**2/M_TWO
            case (POISSON_FFT_KERNEL_CORRECTED)
              fft_Coulb_FS(lx, ly, lz) = M_ZERO
            end select
          end if
        end do
      end do
    end do

    forall(iz=1:cube%fs_n(3), iy=1:cube%fs_n(2), ix=1:cube%fs_n(1))
      fft_Coulb_FS(ix, iy, iz) = M_FOUR*M_PI*fft_Coulb_FS(ix, iy, iz)
    end forall

    call dfourier_space_op_init(this%coulb, cube, fft_coulb_fs, in_device = (kernel /= POISSON_FFT_KERNEL_CORRECTED))

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    POP_SUB(poisson_fft_build_3d_0d)
  end subroutine poisson_fft_build_3d_0d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  !> A. Castro et al., Phys. Rev. B 80, 033102 (2009)
  subroutine poisson_fft_build_2d_0d(this, mesh, cube)
    type(poisson_fft_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(cube_t),        intent(inout) :: cube

    type(spline_t) :: besselintf
    integer :: i, ix, iy, ixx(2), db(2), npoints
    FLOAT :: temp(2), vec, r_c, maxf, dk, default_r_c
    FLOAT, allocatable :: x(:), y(:)
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_2d_0d)

    db(1:2) = cube%rs_n_global(1:2)

    default_r_c = maxval(db(1:2)*mesh%spacing(1:2)/M_TWO)
    call get_cutoff(default_r_c, r_c)

    call spline_init(besselintf)

    ! store the fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:cube%fs_n_global(1), 1:cube%fs_n_global(2), 1:cube%fs_n_global(3)))
    fft_Coulb_FS = M_ZERO
    temp(1:2) = M_TWO*M_PI/(db(1:2)*mesh%spacing(1:2))

    maxf = r_c * sqrt((temp(1)*db(1)/2)**2 + (temp(2)*db(2)/2)**2)
    dk = CNST(0.25) ! This seems to be reasonable.
    npoints = nint(maxf/dk)
    SAFE_ALLOCATE(x(1:npoints))
    SAFE_ALLOCATE(y(1:npoints))
    x(1) = M_ZERO
    y(1) = M_ZERO
    do i = 2, npoints
      x(i) = (i-1) * maxf / (npoints-1)
      y(i) = y(i-1) + poisson_cutoff_2D_0D(x(i-1), x(i))
    end do
    call spline_fit(npoints, x, y, besselintf)

    do iy = 1, cube%fs_n_global(2)
      ixx(2) = pad_feq(iy, db(2), .true.)
      do ix = 1, cube%fs_n_global(1)
        ixx(1) = pad_feq(ix, db(1), .true.)
        vec = sqrt( (temp(1)*ixx(1))**2 + (temp(2)*ixx(2))**2)
        if (vec > M_ZERO) then
           fft_coulb_fs(ix, iy, 1) = (M_TWO * M_PI / vec) * spline_eval(besselintf, vec*r_c)
        else
           fft_coulb_fs(ix, iy, 1) = M_TWO * M_PI * r_c
        end if
      end do
    end do

    call dfourier_space_op_init(this%coulb, cube, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(y)
    call spline_end(besselintf)
    POP_SUB(poisson_fft_build_2d_0d)
  end subroutine poisson_fft_build_2d_0d
  !-----------------------------------------------------------------
    

  !-----------------------------------------------------------------
  !> A. Castro et al., Phys. Rev. B 80, 033102 (2009)
  subroutine poisson_fft_build_2d_1d(this, mesh, cube)
    type(poisson_fft_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(cube_t),        intent(inout) :: cube

    integer :: ix, iy, ixx(2), db(2)
    FLOAT :: temp(2), r_c, gx, gy, default_r_c
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_2d_1d)

    db(1:2) = cube%rs_n_global(1:2)

    default_r_c = db(2)*mesh%spacing(2)/M_TWO
    call get_cutoff(default_r_c, r_c)

    ! store the fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:cube%fs_n_global(1), 1:cube%fs_n_global(2), 1:cube%fs_n_global(3)))
    fft_Coulb_FS = M_ZERO
    temp(1:2) = M_TWO*M_PI/(db(1:2)*mesh%spacing(1:2))

    ! First, the term ix = 0 => gx = 0.
    fft_coulb_fs(1, 1, 1) = -M_FOUR * r_c * (log(r_c)-M_ONE)
    do iy = 2, cube%fs_n_global(2)
      ixx(2) = pad_feq(iy, db(2), .true.)
      gy = temp(2)*ixx(2)
      fft_coulb_fs(1, iy, 1) = -M_FOUR * poisson_cutoff_intcoslog(r_c, gy, M_ONE )
    end do

    do ix = 2, cube%fs_n_global(1)
      ixx(1) = pad_feq(ix, db(1), .true.)
      gx = temp(1)*ixx(1)
      do iy = 1, db(2)
        ixx(2) = pad_feq(iy, db(2), .true.)
        gy = temp(2)*ixx(2)
        fft_coulb_fs(ix, iy, 1) = poisson_cutoff_2d_1d(gy, gx, r_c)
      end do
    end do

    call dfourier_space_op_init(this%coulb, cube, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)

    POP_SUB(poisson_fft_build_2d_1d)
  end subroutine poisson_fft_build_2d_1d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  !> A. Castro et al., Phys. Rev. B 80, 033102 (2009)
  subroutine poisson_fft_build_2d_2d(this, mesh, cube)
    type(poisson_fft_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(cube_t),        intent(inout) :: cube

    integer :: ix, iy, ixx(2), db(2)
    FLOAT :: temp(2), vec
    FLOAT, allocatable :: fft_coulb_FS(:,:,:)

    PUSH_SUB(poisson_fft_build_2d_2d)

    db(1:2) = cube%rs_n_global(1:2)

    ! store the fourier transform of the Coulomb interaction
    SAFE_ALLOCATE(fft_Coulb_FS(1:cube%fs_n_global(1), 1:cube%fs_n_global(2), 1:cube%fs_n_global(3)))
    fft_Coulb_FS = M_ZERO
    temp(1:2) = M_TWO*M_PI/(db(1:2)*mesh%spacing(1:2))

    do iy = 1, cube%fs_n_global(2)
      ixx(2) = pad_feq(iy, db(2), .true.)
      do ix = 1, cube%fs_n_global(1)
        ixx(1) = pad_feq(ix, db(1), .true.)
        vec = sqrt( (temp(1)*ixx(1))**2 + (temp(2)*ixx(2))**2)
        if (vec > M_ZERO) fft_coulb_fs(ix, iy, 1) = M_TWO * M_PI / vec
      end do
    end do

    call dfourier_space_op_init(this%coulb, cube, fft_Coulb_FS)

    SAFE_DEALLOCATE_A(fft_Coulb_FS)
    POP_SUB(poisson_fft_build_2d_2d)
  end subroutine poisson_fft_build_2d_2d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  subroutine poisson_fft_build_1d_1d(this, mesh, cube, poisson_soft_coulomb_param)
    type(poisson_fft_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(cube_t),        intent(inout) :: cube
    FLOAT,               intent(in)    :: poisson_soft_coulomb_param

    integer            :: ix
    FLOAT              :: g
    FLOAT, allocatable :: fft_coulb_fs(:, :, :)

    PUSH_SUB(poisson_fft_build_1d_1d)

    SAFE_ALLOCATE(fft_coulb_fs(1:cube%fs_n_global(1), 1:cube%fs_n_global(2), 1:cube%fs_n_global(3)))
    fft_coulb_fs = M_ZERO

    ! Fourier transform of Soft Coulomb interaction.
    do ix = 1, cube%fs_n_global(1)
      g = ix*M_PI/mesh%sb%lsize(1) ! note that g is always positive with this definition
      fft_coulb_fs(ix, 1, 1) = M_TWO * loct_bessel_k0(poisson_soft_coulomb_param*g)
    end do

    call dfourier_space_op_init(this%coulb, cube, fft_coulb_fs)
    SAFE_DEALLOCATE_A(fft_coulb_fs)
    
    POP_SUB(poisson_fft_build_1d_1d)
  end subroutine poisson_fft_build_1d_1d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  subroutine poisson_fft_build_1d_0d(this, mesh, cube, poisson_soft_coulomb_param)
    type(poisson_fft_t), intent(inout) :: this
    type(mesh_t),        intent(in)    :: mesh
    type(cube_t),        intent(inout) :: cube
    FLOAT,               intent(in)    :: poisson_soft_coulomb_param

    integer            :: box(1), ixx(1), ix
    FLOAT              :: temp(1), g, r_c, default_r_c
    FLOAT, allocatable :: fft_coulb_fs(:, :, :)

    PUSH_SUB(poisson_fft_build_1d_0d)

    box(1:1) = cube%rs_n_global(1:1)

    default_r_c = box(1)*mesh%spacing(1)/M_TWO
    call get_cutoff(default_r_c, r_c)

    SAFE_ALLOCATE(fft_coulb_fs(1:cube%fs_n_global(1), 1:cube%fs_n_global(2), 1:cube%fs_n_global(3)))
    fft_coulb_fs = M_ZERO
    temp(1:1) = M_TWO*M_PI/(box(1:1)*mesh%spacing(1:1))

    ! Fourier transform of Soft Coulomb interaction.
    do ix = 1, cube%fs_n_global(1)
      ixx(1) = pad_feq(ix, box(1), .true.)
      g = temp(1)*ixx(1)
      fft_coulb_fs(ix, 1, 1) = poisson_cutoff_1D_0D(g, poisson_soft_coulomb_param, r_c)
    end do

    call dfourier_space_op_init(this%coulb, cube, fft_coulb_fs)
    SAFE_DEALLOCATE_A(fft_coulb_fs)
    
    POP_SUB(poisson_fft_build_1d_0d)
  end subroutine poisson_fft_build_1d_0d
  !-----------------------------------------------------------------


  !-----------------------------------------------------------------
  subroutine poisson_fft_end(this)
    type(poisson_fft_t), intent(inout) :: this

    PUSH_SUB(poisson_fft.end)

    call fourier_space_op_end(this%coulb)

    POP_SUB(poisson_fft.end)
  end subroutine poisson_fft_end

  !-----------------------------------------------------------------

  subroutine poisson_fft_solve(this, mesh, cube, pot, rho, mesh_cube_map, average_to_zero)
    type(poisson_fft_t),            intent(inout) :: this
    type(mesh_t),                   intent(in)    :: mesh
    type(cube_t),                   intent(inout) :: cube
    FLOAT,                          intent(out)   :: pot(:)
    FLOAT,                          intent(in)    :: rho(:)
    type(mesh_cube_parallel_map_t), intent(in)    :: mesh_cube_map
    logical,              optional, intent(in)    :: average_to_zero !< default is false

    logical :: average_to_zero_
    FLOAT :: average
    type(cube_function_t) :: cf

    PUSH_SUB(poisson_fft)
    
    average_to_zero_ = .false.
    if (present(average_to_zero)) average_to_zero_ = average_to_zero
    average = M_ZERO !this avoids a non-initialized warning

    call cube_function_null(cf)    
    call dcube_function_alloc_RS(cube, cf, in_device = (this%kernel /= POISSON_FFT_KERNEL_CORRECTED)) 

    ! put the density in the cube
    if (cube%parallel_in_domains) then
      call dmesh_to_cube_parallel(mesh, rho, cube, cf, mesh_cube_map)
    else
      if(mesh%parallel_in_domains) then
        call dmesh_to_cube(mesh, rho, cube, cf, local = .true.)
      else 
        call dmesh_to_cube(mesh, rho, cube, cf)
      end if
    end if

    ! apply the Couloumb term in Fourier space
    call dfourier_space_op_apply(this%coulb, cube, cf)

    !now the cube has the potential
    if(average_to_zero_) average = cube_function_surface_average(cube, cf)
    
    ! move the potential back to the mesh
    if (cube%parallel_in_domains) then
      call dcube_to_mesh_parallel(cube, cf, mesh, pot, mesh_cube_map)
    else
      if(mesh%parallel_in_domains) then
        call dcube_to_mesh(cube, cf, mesh, pot, local=.true.)
      else
        call dcube_to_mesh(cube, cf, mesh, pot)
      end if
    end if

    if(average_to_zero_) pot(1:mesh%np) = pot(1:mesh%np) - average
    
    call dcube_function_free_RS(cube, cf) ! memory is no longer needed

    POP_SUB(poisson_fft)
  end subroutine poisson_fft_solve
  
end module poisson_fft_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
