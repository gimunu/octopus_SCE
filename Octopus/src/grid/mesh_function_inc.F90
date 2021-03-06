!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Verstraete
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
!! $Id: mesh_function_inc.F90 11008 2013-07-14 19:01:51Z acastro $


! ---------------------------------------------------------
!> integrates a function
R_TYPE function X(mf_integrate) (mesh, ff) result(dd)
  type(mesh_t), intent(in) :: mesh
  R_TYPE,       intent(in) :: ff(:)  !< (mesh%np)

  integer :: ip

  call profiling_in(C_PROFILING_MF_INTEGRATE, 'MF_INTEGRATE')
  PUSH_SUB(X(mf_integrate))

  ASSERT(ubound(ff, dim = 1) == mesh%np .or. ubound(ff, dim = 1) == mesh%np_part)

  dd = M_ZERO
  if (mesh%use_curvilinear) then
    do ip = 1, mesh%np
      dd = dd + ff(ip)*mesh%vol_pp(ip)
    end do
  else
    do ip = 1, mesh%np
      dd = dd + ff(ip)
    end do
  end if

  dd = dd*mesh%volume_element

  if(mesh%parallel_in_domains) then
    call profiling_in(C_PROFILING_MF_REDUCE, "MF_REDUCE")
    call comm_allreduce(mesh%mpi_grp%comm, dd)
    call profiling_out(C_PROFILING_MF_REDUCE)
  end if
  
  POP_SUB(X(mf_integrate))
  call profiling_out(C_PROFILING_MF_INTEGRATE)

end function X(mf_integrate)


!> ---------------------------------------------------------
!! This function returns the dot product between two vectors,
!! but using the mesh_aux defined as a global object in this
!! module. This way it can be called by external libraries,
!! passing only the two vectors. First, one has to 
!! make sure that mesh_aux is pointing to some defined
!! mesh data structure, by calling mesh_init_mesh_aux.
!! ---------------------------------------------------------
R_TYPE function X(mf_dotp_aux)(f1, f2) result(dotp)
  R_TYPE,            intent(in) :: f1(:), f2(:)

  PUSH_SUB(X(mf_dotp_aux))

  ASSERT(associated(mesh_aux))
  dotp = X(mf_dotp)(mesh_aux, f1, f2)

  POP_SUB(X(mf_dotp_aux))
end function X(mf_dotp_aux)

!> Same as above, but no conjugation.
!! ---------------------------------------------------------
R_TYPE function X(mf_dotu_aux)(f1, f2) result(dotu)
  R_TYPE,            intent(in) :: f1(:), f2(:)

  PUSH_SUB(X(mf_dotu_aux))

  ASSERT(associated(mesh_aux))
  dotu = X(mf_dotp)(mesh_aux, f1, f2, dotu = .true.)

  POP_SUB(X(mf_dotu_aux))
end function X(mf_dotu_aux)

!> Same as above, but for norm.
!! ---------------------------------------------------------
FLOAT function X(mf_nrm2_aux)(ff) result(norm)
  R_TYPE,            intent(in) :: ff(:)

  PUSH_SUB(X(mf_nrm2_aux))

  ASSERT(associated(mesh_aux))
  norm = X(mf_nrm2)(mesh_aux, ff)

  POP_SUB(X(mf_nrm2_aux))
end function X(mf_nrm2_aux)


! ---------------------------------------------------------
!> this function returns the dot product between two vectors
R_TYPE function X(mf_dotp_1)(mesh, f1, f2, reduce, dotu) result(dotp)
  type(mesh_t),      intent(in) :: mesh
  R_TYPE,            intent(in) :: f1(:), f2(:)
  logical, optional, intent(in) :: reduce
  logical, optional, intent(in) :: dotu
     !< if true, use blas_dotu instead of blas_dot;
     !! no complex conjugation.  Default is false.
     !! has no effect if working with real version

#ifdef R_TCOMPLEX
  logical             :: dotu_
#endif
  integer             :: ip

  call profiling_in(C_PROFILING_MF_DOTP, "MF_DOTP")
  PUSH_SUB(X(mf_dotp_1))

  ASSERT(ubound(f1, dim = 1) == mesh%np .or. ubound(f1, dim = 1) == mesh%np_part)
  ASSERT(ubound(f2, dim = 1) == mesh%np .or. ubound(f2, dim = 1) == mesh%np_part)

#ifdef R_TCOMPLEX
  dotu_ = optional_default(dotu, .false.)
#endif

  if(mesh%use_curvilinear) then
    dotp = M_ZERO
    ! preprocessor conditionals necessary since blas_dotu only exists for complex input
#ifdef R_TCOMPLEX
    if (.not. dotu_) then
#endif
      do ip = 1, mesh%np
        dotp = dotp + mesh%vol_pp(ip)*R_CONJ(f1(ip))*f2(ip)
      end do
#ifdef R_TCOMPLEX
    else
      do ip = 1, mesh%np
        dotp = dotp + mesh%vol_pp(ip)*f1(ip)*f2(ip)
      end do
    endif
#endif
    call profiling_count_operations(mesh%np*(2*R_ADD + R_MUL))
  else
#ifdef R_TCOMPLEX
    if (.not. dotu_) then
#endif
      dotp = blas_dot(mesh%np, f1(1), 1, f2(1), 1)
#ifdef R_TCOMPLEX
    else
      dotp = blas_dotu(mesh%np, f1(1), 1, f2(1), 1)
    endif
#endif
    call profiling_count_operations(mesh%np*(R_ADD + R_MUL))

  end if

  dotp = dotp*mesh%volume_element

  if(mesh%parallel_in_domains .and. optional_default(reduce, .true.)) then
    call profiling_in(C_PROFILING_MF_REDUCE, "MF_REDUCE")
    call comm_allreduce(mesh%vp%comm, dotp)
    call profiling_out(C_PROFILING_MF_REDUCE)
  end if

  POP_SUB(X(mf_dotp_1))
  call profiling_out(C_PROFILING_MF_DOTP)

end function X(mf_dotp_1)


! ---------------------------------------------------------
R_TYPE function X(mf_dotp_2)(mesh, dim, f1, f2, reduce, dotu) result(dotp)
  type(mesh_t),      intent(in) :: mesh
  integer,           intent(in) :: dim
  R_TYPE,            intent(in) :: f1(:,:), f2(:,:)
  logical, optional, intent(in) :: reduce
  logical, optional, intent(in) :: dotu
     !< if true, use lalg_dotu instead of lalg_dot;
     !! no complex conjugation.  Default is false.

  integer :: idim

  PUSH_SUB(X(mf_dotp_2))

  dotp = R_TOTYPE(M_ZERO)
  do idim = 1, dim
    dotp = dotp + X(mf_dotp_1)(mesh, f1(:, idim), f2(:, idim), reduce = .false., dotu = dotu)
  end do

  if(mesh%parallel_in_domains .and. optional_default(reduce, .true.)) then
    call profiling_in(C_PROFILING_MF_REDUCE, "MF_REDUCE")
    call comm_allreduce(mesh%vp%comm, dotp)
    call profiling_out(C_PROFILING_MF_REDUCE)
  end if

  POP_SUB(X(mf_dotp_2))

end function X(mf_dotp_2)


! ---------------------------------------------------------
!> this function returns the the norm of a vector
FLOAT function X(mf_nrm2_1)(mesh, ff, reduce) result(nrm2)
  type(mesh_t),      intent(in) :: mesh
  R_TYPE,            intent(in) :: ff(:)
  logical, optional, intent(in) :: reduce
 
  R_TYPE, allocatable :: ll(:)

  call profiling_in(C_PROFILING_MF_NRM2, "MF_NRM2")
  PUSH_SUB(X(mf_nrm2_1))

  if(mesh%use_curvilinear) then
    SAFE_ALLOCATE(ll(1:mesh%np))
    ll(1:mesh%np) = ff(1:mesh%np)*sqrt(mesh%vol_pp(1:mesh%np))
    nrm2 = lalg_nrm2(mesh%np, ll)
    SAFE_DEALLOCATE_A(ll)
  else
    nrm2 = lalg_nrm2(mesh%np, ff)
  end if

  nrm2 = nrm2*sqrt(mesh%volume_element)

  if(mesh%parallel_in_domains .and. optional_default(reduce, .true.)) then
    call profiling_in(C_PROFILING_MF_REDUCE, "MF_REDUCE")
    nrm2 = nrm2**2
    call comm_allreduce(mesh%vp%comm, nrm2)
    nrm2 = sqrt(nrm2)
    call profiling_out(C_PROFILING_MF_REDUCE)
  end if

  POP_SUB(X(mf_nrm2_1))
  call profiling_out(C_PROFILING_MF_NRM2)

end function X(mf_nrm2_1)

! ---------------------------------------------------------
FLOAT function X(mf_nrm2_2)(mesh, dim, ff, reduce) result(nrm2)
  type(mesh_t),      intent(in) :: mesh
  integer,           intent(in) :: dim
  R_TYPE,            intent(in) :: ff(:,:)
  logical, optional, intent(in) :: reduce

  integer :: idim

  PUSH_SUB(X(mf_nrm2_2))

  nrm2 = M_ZERO

  do idim = 1, dim
    nrm2 = hypot(nrm2, X(mf_nrm2)(mesh, ff(:, idim), reduce = reduce))
  end do

  POP_SUB(X(mf_nrm2_2))

end function X(mf_nrm2_2)


! ---------------------------------------------------------
!> This function calculates the "order" moment of the function ff
R_TYPE function X(mf_moment) (mesh, ff, idir, order) result(rr)
  type(mesh_t), intent(in) :: mesh
  R_TYPE,       intent(in) :: ff(:)
  integer,      intent(in) :: idir
  integer,      intent(in) :: order

  R_TYPE, allocatable :: fxn(:)

  PUSH_SUB(X(mf_moment))

  SAFE_ALLOCATE(fxn(1:mesh%np))

  fxn(1:mesh%np) = ff(1:mesh%np)*mesh%x(1:mesh%np, idir)**order
  rr = X(mf_integrate)(mesh, fxn)

  SAFE_DEALLOCATE_A(fxn)
  POP_SUB(X(mf_moment))

end function X(mf_moment)


! ---------------------------------------------------------
!> This subroutine generates a Gaussian wavefunction at a
!! random position in space.
subroutine X(mf_random)(mesh, ff, seed)
  type(mesh_t),      intent(in)  :: mesh
  R_TYPE,            intent(out) :: ff(:)
  integer, optional, intent(in)  :: seed

  integer, save :: iseed = 123
  integer :: idim, ip
  FLOAT :: aa(MAX_DIM), rnd, rr

  PUSH_SUB(X(mf_random))

  if(present(seed)) then
    iseed = iseed + seed
  end if

  aa = M_ZERO
  do idim = 1, mesh%sb%dim
    call quickrnd(iseed, rnd)
    aa(idim) = M_TWO*(2*rnd - 1) * M_FOUR * mesh%spacing(idim)
  end do

  !$omp parallel do private(rr)
  do ip = 1, mesh%np
    rr = sum( ((mesh%x(ip, 1:mesh%sb%dim) - aa(1:mesh%sb%dim)) / (M_FOUR * mesh%spacing(1:mesh%sb%dim)) ) **2)
    if ( rr < CNST(100.0) ) then 
      ff(ip) = exp(-M_HALF*rr)
    else
      ff(ip) = M_ZERO
    end if
  end do
  !$omp end parallel do

  rr = X(mf_nrm2)(mesh, ff)
  call lalg_scal(mesh%np, M_ONE/rr, ff)

  POP_SUB(X(mf_random))

end subroutine X(mf_random)


! --------------------------------------------------------- 
!> This function receives a function f_in defined in a mesh, and returns
!! the interpolated values of the function over the npoints_in defined
!! by x_in.

subroutine X(mf_interpolate_points) (ndim, npoints_in, x_in, f_in, npoints_out, x_out, f_out)
  integer,              intent(in)  :: ndim, npoints_in, npoints_out
  R_TYPE, target,       intent(in)  :: f_in(:)    !< (npoints_in)
  FLOAT,  target,       intent(in)  :: x_in(:, :)
  FLOAT,                intent(in)  :: x_out(:,:)
  R_TYPE,               intent(out) :: f_out(:)   !< (npoints_out)

  real(8) :: pp(MAX_DIM)
  R_DOUBLE, pointer :: rf_in(:)
  real(8),  pointer :: rx_in(:, :)
  integer :: ip
  type(qshep_t) :: interp
#ifndef R_TCOMPLEX
  type(spline_t) :: interp1d
#endif

  PUSH_SUB(X(mf_interpolate_points))

#ifdef SINGLE_PRECISION
  SAFE_ALLOCATE(rx_in(1:npoints_in, 1:ndim))
  rx_in = x_in
  SAFE_ALLOCATE(rf_in(1:npoints_in))
  rf_in = f_in
#else
  rx_in => x_in
  rf_in => f_in
#endif

  select case(ndim)
  case(2)
    call init_qshep(interp, npoints_in, rf_in, rx_in(:, 1), rx_in(:, 2))
    do ip = 1, npoints_out
      pp(1:2)   = x_out(ip, 1:2)
      f_out(ip) = qshep_interpolate(interp, rf_in, pp(1:2))
    end do
    call kill_qshep(interp)

  case(3)
    call init_qshep(interp, npoints_in, rf_in, rx_in(:, 1), rx_in(:, 2), rx_in(:, 3))
    do ip = 1, npoints_out
      pp(1:3)   = x_out(ip, 1:3)
      f_out(ip) = qshep_interpolate(interp, rf_in, pp(1:3))
    end do
    call kill_qshep(interp)

  case(1)
#ifdef R_TCOMPLEX
    message(1) = 'Believe it or not: cannot do 1D complex interpolation, only 2D or 3D.'
    call messages_fatal(1)
#else
    call spline_init(interp1d)
    call spline_fit(npoints_in, R_REAL(rx_in(:, 1)), rf_in, interp1d)
    do ip = 1, npoints_out
      f_out(ip) = spline_eval(interp1d, x_out(ip, 1))
    end do
    call spline_end(interp1d)
#endif
  end select

#ifdef SINGLE_PRECISION
  SAFE_DEALLOCATE_P(rf_in)
  SAFE_DEALLOCATE_P(rx_in)
#endif

  POP_SUB(X(mf_interpolate_points))
end subroutine X(mf_interpolate_points)


! ---------------------------------------------------------
!> Given a function ff defined on mesh, and a plane, it gives 
!! back the values of ff on the plane, by doing the appropriate
!! interpolation.
subroutine X(mf_interpolate_on_plane)(mesh, plane, ff, f_in_plane)
  type(mesh_t),       intent(in)  :: mesh
  type(mesh_plane_t), intent(in)  :: plane
  R_TYPE,             intent(in)  :: ff(:)
  R_TYPE,             intent(out) :: f_in_plane(plane%nu:plane%mu, plane%nv:plane%mv)

  integer :: iu, iv, ip
  R_DOUBLE, allocatable :: f_global(:)
  real(8) :: pp(3)
  type(qshep_t) :: interp
  real(8), allocatable :: xglobal(:, :)

  PUSH_SUB(X(mf_interpolate_on_plane))

  SAFE_ALLOCATE(xglobal(1:mesh%np_part_global, 1:MAX_DIM))
  do ip = 1, mesh%np_part_global
    xglobal(ip, 1:) = mesh_x_global(mesh, ip)
  end do

  SAFE_ALLOCATE(f_global(1:mesh%np_global))
#if defined HAVE_MPI
  call X(vec_gather)(mesh%vp, mesh%vp%root, f_global, ff)
#else
  f_global(1:mesh%np_global) = ff(1:mesh%np_global)
#endif

  call init_qshep(interp, mesh%np_global, f_global, xglobal(:, 1), xglobal(:, 2), xglobal(:, 3) )

  do iu = plane%nu, plane%mu
    do iv = plane%nv, plane%mv
      pp(1) = plane%origin(1) + iu*plane%spacing * plane%u(1) + iv * plane%spacing * plane%v(1)
      pp(2) = plane%origin(2) + iu*plane%spacing * plane%u(2) + iv * plane%spacing * plane%v(2)
      pp(3) = plane%origin(3) + iu*plane%spacing * plane%u(3) + iv * plane%spacing * plane%v(3)
      f_in_plane(iu, iv) = qshep_interpolate(interp, f_global, pp(1:3))
    end do
  end do

  call kill_qshep(interp)

  SAFE_DEALLOCATE_A(xglobal)
  SAFE_DEALLOCATE_A(f_global)
  POP_SUB(X(mf_interpolate_on_plane))
end subroutine X(mf_interpolate_on_plane)


! ---------------------------------------------------------
!> Given a function ff defined on mesh, and a line, it gives 
!! back the values of ff on the line, by doing the appropriate
!! interpolation.
subroutine X(mf_interpolate_on_line)(mesh, line, ff, f_in_line)
  type(mesh_t),       intent(in)  :: mesh
  type(mesh_line_t),  intent(in)  :: line
  R_TYPE,             intent(in)  :: ff(:)
  R_TYPE,             intent(out) :: f_in_line(line%nu:line%mu)

  integer :: iu, ip
  R_DOUBLE, allocatable :: f_global(:)
  real(8) :: pp(2)
  type(qshep_t) :: interp
  real(8), allocatable :: xglobal(:, :)

  PUSH_SUB(X(mf_interpolate_on_line))

  SAFE_ALLOCATE(xglobal(1:mesh%np_part_global, 1:MAX_DIM))
  do ip = 1, mesh%np_part_global
    xglobal(ip, 1:MAX_DIM) = mesh_x_global(mesh, ip)
  end do
  
  SAFE_ALLOCATE(f_global(1:mesh%np_global))
#if defined HAVE_MPI
  call X(vec_gather)(mesh%vp, mesh%vp%root, f_global, ff)
#else
  f_global(1:mesh%np_global) = ff(1:mesh%np_global)
#endif

  call init_qshep(interp, mesh%np_global, f_global, xglobal(:, 1), xglobal(:, 2))
  do iu = line%nu, line%mu
    pp(1) = line%origin(1) + iu * line%spacing * line%u(1)
    pp(2) = line%origin(2) + iu * line%spacing * line%u(2)
    f_in_line(iu) = qshep_interpolate(interp, f_global, pp(1:2))
  end do
  call kill_qshep(interp)

  SAFE_DEALLOCATE_A(f_global)
  SAFE_DEALLOCATE_A(xglobal)

  POP_SUB(X(mf_interpolate_on_line))
end subroutine X(mf_interpolate_on_line)


! ---------------------------------------------------------
!> This subroutine calculates the surface integral of a scalar
!! function on a given plane.
R_TYPE function X(mf_surface_integral_scalar) (mesh, ff, plane) result(dd)
  type(mesh_t),       intent(in) :: mesh
  R_TYPE,             intent(in) :: ff(:)  !< (mesh%np)
  type(mesh_plane_t), intent(in) :: plane

  R_TYPE, allocatable :: f_in_plane(:, :)

  PUSH_SUB(X(mf_surface_integral_scalar))

  if(mesh%sb%dim /= 3) then
    message(1) = 'INTERNAL ERROR at Xmf_surface_integral: wrong dimensionality.'
    call messages_fatal(1)
  end if

  SAFE_ALLOCATE(f_in_plane(plane%nu:plane%mu, plane%nv:plane%mv))

  call X(mf_interpolate_on_plane)(mesh, plane, ff, f_in_plane)

  dd = sum(f_in_plane(:, :) * plane%spacing**2)

  SAFE_DEALLOCATE_A(f_in_plane)
  POP_SUB(X(mf_surface_integral_scalar))
end function X(mf_surface_integral_scalar)


! ---------------------------------------------------------
!> This subroutine calculates the surface integral of a vector
!! function on a given plane.
R_TYPE function X(mf_surface_integral_vector) (mesh, ff, plane) result(dd)
  type(mesh_t), intent(in)       :: mesh
  R_TYPE,       intent(in)       :: ff(:, :)  !< (mesh%np, MAX_DIM)
  type(mesh_plane_t), intent(in) :: plane

  R_TYPE, allocatable :: fn(:)
  integer :: ip

  PUSH_SUB(X(mf_surface_integral_vector))

  SAFE_ALLOCATE(fn(1:mesh%np))
  do ip = 1, mesh%np
    fn(ip) = sum(ff(ip, :) * plane%n(:))
  end do

  dd =  X(mf_surface_integral_scalar)(mesh, fn, plane)

  SAFE_DEALLOCATE_A(fn)
  POP_SUB(X(mf_surface_integral_vector))
end function X(mf_surface_integral_vector)


! ---------------------------------------------------------
!> This subroutine calculates the line integral of a scalar
!! function on a given line.
R_TYPE function X(mf_line_integral_scalar) (mesh, ff, line) result(dd)
  type(mesh_t),      intent(in) :: mesh
  R_TYPE,            intent(in) :: ff(:)  !< (mesh%np)
  type(mesh_line_t), intent(in) :: line

  R_TYPE, allocatable :: f_in_line(:)

  PUSH_SUB(X(mf_line_integral_scalar))

  if(mesh%sb%dim /= 2) then
    message(1) = 'INTERNAL ERROR at Xmf_surface_integral: wrong dimensionality.'
    call messages_fatal(1)
  end if

  SAFE_ALLOCATE(f_in_line(line%nu:line%mu))

  call X(mf_interpolate_on_line)(mesh, line, ff, f_in_line)

  dd = sum(f_in_line(:) * line%spacing)

  SAFE_DEALLOCATE_A(f_in_line)
  POP_SUB(X(mf_line_integral_scalar))
end function X(mf_line_integral_scalar)


! ---------------------------------------------------------
!> This subroutine calculates the line integral of a vector
!! function on a given line.
R_TYPE function X(mf_line_integral_vector) (mesh, ff, line) result(dd)
  type(mesh_t),      intent(in) :: mesh
  R_TYPE,            intent(in) :: ff(:, :)  !< (mesh%np, MAX_DIM)
  type(mesh_line_t), intent(in) :: line

  R_TYPE, allocatable :: fn(:)
  integer :: ip

  PUSH_SUB(X(mf_line_integral_vector))

  SAFE_ALLOCATE(fn(1:mesh%np))
  do ip = 1, mesh%np
    fn(ip) = sum(ff(ip, 1:mesh%sb%dim) * line%n(1:mesh%sb%dim))
  end do

  dd = X(mf_line_integral_scalar)(mesh, fn, line)

  SAFE_DEALLOCATE_A(fn)
  POP_SUB(X(mf_line_integral_vector))
end function X(mf_line_integral_vector)


! ---------------------------------------------------------
!> Converts a spline that represents a radial function into a mesh.
subroutine X(mf_put_radial_spline)(mesh, spl, center, ff, add)
  type(mesh_t),        intent(in)    :: mesh
  type(spline_t),      intent(in)    :: spl
  FLOAT,               intent(in)    :: center(:)
  R_TYPE,              intent(inout) :: ff(:)
  logical, optional,   intent(in)    :: add

  integer :: ip
  FLOAT :: rr
  logical :: add_

  PUSH_SUB(X(mf_put_radial_spline))

  add_ = .false.

  if(present(add)) then
    add_ = add
  endif

  if( add_ ) then 

    do ip = 1, mesh%np
      rr = sqrt(sum((mesh%x(ip, 1:mesh%sb%dim) - center(1:mesh%sb%dim))**2))
      ff(ip) = ff(ip) + spline_eval(spl, rr)
    end do

  else

    do ip = 1, mesh%np
      rr = sqrt(sum((mesh%x(ip, 1:mesh%sb%dim) - center(1:mesh%sb%dim))**2))
      ff(ip) = spline_eval(spl, rr)
    end do

  end if

  POP_SUB(X(mf_put_radial_spline))
end subroutine X(mf_put_radial_spline)


! -----------------------------------------------------------------------------
!> This routine calculates the multipoles of a function ff,
!! defined in the following way:
!! multipole(1) is the trace of ff (defined to be positive; integral
!!   of ff).
!! multipole(2:4) contains the dipole: integral of ff times x, y or z.
!! multipole(5:9, is) contains the quadrupole, defined in the usual way using
!!   the spherical harmonics: multipole(5) = Integral [ ff * Y_{2,-2} ],
!!   multipole(6, is) = Integral [ f * Y_{2, -1} ].
!! And so on.
!! -----------------------------------------------------------------------------
subroutine X(mf_multipoles) (mesh, ff, lmax, multipole, cmplxscl_th)
  type(mesh_t),      intent(in)  :: mesh
  R_TYPE,            intent(in)  :: ff(:)
  integer,           intent(in)  :: lmax
  R_TYPE,            intent(out) :: multipole(:) !< ((lmax + 1)**2)
  FLOAT, optional,   intent(in)  :: cmplxscl_th !< the space complex scaling angle cmplxscl%theta 

  integer :: idim, ip, ll, lm, add_lm
  FLOAT   :: xx(MAX_DIM), rr, ylm
  R_TYPE, allocatable :: ff2(:)
  logical :: cmplxscl

  PUSH_SUB(X(mf_multipoles))
  cmplxscl = .false.
  if(present(cmplxscl_th)) cmplxscl = .true.


  ASSERT(ubound(ff, dim = 1) == mesh%np .or. ubound(ff, dim = 1) == mesh%np_part)

  SAFE_ALLOCATE(ff2(1:mesh%np))

  ff2(1:mesh%np) = ff(1:mesh%np)
  multipole(1) = X(mf_integrate)(mesh, ff2)
  
  if(.not. cmplxscl) then
    
    if(lmax > 0) then
      do idim = 1, 3
        ff2(1:mesh%np) = ff(1:mesh%np) * mesh%x(1:mesh%np, idim)
        multipole(idim+1) = X(mf_integrate)(mesh, ff2)
      end do
    end if

    if(lmax>1) then
      add_lm = 5
      do ll = 2, lmax
        do lm = -ll, ll
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, rr, coords=xx)
            call loct_ylm(1, xx(1), xx(2), xx(3), ll, lm, ylm)
            ff2(ip) = ff(ip) * ylm * rr**ll
          end do
          multipole(add_lm) = X(mf_integrate)(mesh, ff2)
          add_lm = add_lm + 1
        end do
      end do
    end if

  else

    if(lmax > 0) then
      do idim = 1, 3
        ff2(1:mesh%np) = ff(1:mesh%np) * mesh%x(1:mesh%np, idim)
        ff2(1:mesh%np) = ff2(1:mesh%np)* exp(M_zI * cmplxscl_th)
        multipole(idim+1) = X(mf_integrate)(mesh, ff2)
      end do
    end if

    if(lmax>1) then
      add_lm = 5
      do ll = 2, lmax
        do lm = -ll, ll
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, rr, coords=xx)
            call loct_ylm(1, xx(1), xx(2), xx(3), ll, lm, ylm)
            ff2(ip) = ff(ip) * ylm * rr**ll * exp(M_zI * cmplxscl_th * ll) 
          end do
          multipole(add_lm) = X(mf_integrate)(mesh, ff2)
          add_lm = add_lm + 1
        end do
      end do
    end if

  end if

  SAFE_DEALLOCATE_A(ff2)
  POP_SUB(X(mf_multipoles))
end subroutine X(mf_multipoles)


!--------------------------------------------------------------
!> add some part of the function in mesh1 to mesh2 (index based, so same spacing is assumed)
subroutine X(mf_add)(mesh1, start1, end1, func1, mesh2, start2, end2, func2, per_dim, include_bounds)
  type(mesh_t),      intent(in)    :: mesh1
  integer,           intent(in)    :: start1(1:3)    !< starting point in mesh1
  integer,           intent(in)    :: end1(1:3)      !< end point in mesh1
  R_TYPE,            intent(in)    :: func1(:)       !< copy from this
  type(mesh_t),      intent(in)    :: mesh2
  integer,           intent(in)    :: start2(1:3)    !< starting point in mesh2
  integer,           intent(in)    :: end2(1:3)      !< end point in mesh2
  R_TYPE,            intent(inout) :: func2(:)       !< add to this
  integer,           intent(in)    :: per_dim        !< the periodic dimension
  logical, optional, intent(in)    :: include_bounds !< also use the points outside the box?

  integer :: np1, np2, ip1, ip2, ix2, iy2, iz2
  integer :: ix1(1:3)

  PUSH_SUB(X(mf_add))

  if (optional_default(include_bounds, .false.)) then
    ASSERT(per_dim == 0)
    np1 = mesh1%np_part
    np2 = mesh2%np_part
  else
    np1 = mesh1%np
    np2 = mesh2%np
  end if

  ix1(1) = start1(1)
  do ix2 = start2(1), end2(1)

    ix1(2) = start1(2)
    do iy2 = start2(2), end2(2)

      ix1(3) = start1(3)
      do iz2 = start2(3), end2(3)

        ip1 = mesh1%idx%lxyz_inv(ix1(1), ix1(2), ix1(3))
        ip2 = mesh2%idx%lxyz_inv(ix2, iy2, iz2)

        if(ip2 > 0 .and. ip2 <= np2 .and. ip1 > 0 .and. ip1 <= np1) then
          func2(ip2) = func2(ip2) + func1(ip1)
        end if

        INCR(ix1(3), 1)
        if(per_dim == 3 .and. ix1(3) > end1(3)) ix1(3) = start1(3)
      end do

      INCR(ix1(2), 1)
      if(per_dim >= 2 .and. ix1(2) > end1(2)) ix1(2) = start1(2)
    end do

    INCR(ix1(1), 1)
    if(per_dim >= 1 .and. ix1(1) > end1(1)) ix1(1) = start1(1)
  end do

  POP_SUB(X(mf_add))
end subroutine X(mf_add)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
