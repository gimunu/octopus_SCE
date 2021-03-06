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
!! $Id: restart_inc.F90 10978 2013-07-11 15:28:46Z micael $


! ---------------------------------------------------------
subroutine X(restart_write_function)(dir, filename, mesh, ff, ierr)
  character(len=*), intent(in)  :: dir
  character(len=*), intent(in)  :: filename
  type(mesh_t),     intent(in)  :: mesh
  R_TYPE,           intent(in)  :: ff(:)
  integer,          intent(out) :: ierr

  PUSH_SUB(X(restart_write_function))

  call X(io_function_output)(restart_format, trim(dir), trim(filename), mesh, ff(:), unit_one, ierr, is_tmp=.true.)
  ! all restart files are in atomic units

  POP_SUB(X(restart_write_function))
end subroutine X(restart_write_function)


! ---------------------------------------------------------
!> In domain parallel case each process reads a part of the file.
!! At the end all the processes have the corresponding mesh part
subroutine X(restart_read_function)(dir, filename, mesh, ff, ierr, map)
  character(len=*), intent(in)    :: dir
  character(len=*), intent(in)    :: filename
  type(mesh_t),     intent(in)    :: mesh
  R_TYPE, target,   intent(inout) :: ff(:)
  integer,          intent(out)   :: ierr
  integer, optional, intent(in)   :: map(:)

  integer :: ip, np, offset
  R_TYPE, pointer :: read_ff(:)
  type(profile_t), save :: prof_io
  type(batch_t) :: ffb
  type(profile_t), save :: prof_comm

  PUSH_SUB(X(restart_read_function))
  
  nullify(read_ff)

  if(present(map) .and. mesh%parallel_in_domains) then 
    ! for the moment we do not do this directly
    call X(io_function_input) (trim(dir)//'/'//trim(filename)//'.obf', mesh, ff(1:mesh%np), ierr, is_tmp=.true., map = map)

    POP_SUB(X(restart_read_function))
    return
  end if

  if(present(map)) then
    call io_binary_get_info(trim(dir)//'/'//trim(filename)//'.obf', np, ierr)

    if(ierr /= 0) then
      POP_SUB(X(restart_read_function))
      return
    end if

    ASSERT(np > 0)
    SAFE_ALLOCATE(read_ff(1:np))
  else
    np = mesh%np
    read_ff => ff
  end if

  offset = 0
  !in the parallel case, each node reads a part of the file
  if(mesh%parallel_in_domains) then
    offset = mesh%vp%xlocal - 1
  end if

  ASSERT(associated(read_ff))

  call profiling_in(prof_io, "RESTART_READ_IO")

  call io_binary_read(trim(dir)//'/'//trim(filename)//'.obf', np, read_ff, ierr, offset = offset)
  call profiling_count_transfers(np, read_ff(1))
  call profiling_out(prof_io)

  if(mesh%parallel_in_domains) then
    call profiling_in(prof_comm, "RESTART_READ_COMM")
    ! this is the global index of the points we read

    ff(1:mesh%np) = read_ff(1:mesh%np)

    call batch_init(ffb, 1)
    call batch_add_state(ffb, ff)
    call X(mesh_batch_exchange_points)(mesh, ffb, backward_map = .true.)
    call batch_end(ffb)
    
    call profiling_out(prof_comm)
  end if

  if(present(map)) then
    ff(1:mesh%np_global) = M_ZERO
    do ip = 1, min(np, ubound(map, dim = 1))
      if(map(ip) > 0) ff(map(ip)) = read_ff(ip)
    end do
    
    SAFE_DEALLOCATE_P(read_ff)
  end if

  POP_SUB(X(restart_read_function))
end subroutine X(restart_read_function)


! ---------------------------------------------------------
subroutine X(restart_write_lr_rho)(lr, gr, nspin, restart_dir, rho_tag)
  type(lr_t),        intent(in)    :: lr
  type(grid_t),      intent(in)    :: gr
  integer,           intent(in)    :: nspin
  character(len=*),  intent(in)    :: restart_dir
  character(len=*),  intent(in)    :: rho_tag

  character(len=100) :: fname
  integer :: is, ierr

  PUSH_SUB(X(restart_write_lr_rho))

  call block_signals()
  do is = 1, nspin
    write(fname, '(a,i1,a)') trim(rho_tag)//'_', is
    call X(restart_write_function)(trim(tmpdir)//trim(RESTART_DIR), fname, gr%mesh, lr%X(dl_rho)(:, is), ierr)
  end do
  call unblock_signals()

  POP_SUB(X(restart_write_lr_rho))
end subroutine X(restart_write_lr_rho)


! ---------------------------------------------------------
subroutine X(restart_read_lr_rho)(dl_rho, gr, nspin, restart_subdir, rho_tag, ierr)
  R_TYPE,            intent(inout) :: dl_rho(:,:) !< (gr%mesh%np, nspin)
  type(grid_t),      intent(in)    :: gr
  integer,           intent(in)    :: nspin
  character(len=*),  intent(in)    :: restart_subdir
  character(len=*),  intent(in)    :: rho_tag
  integer,           intent(out)   :: ierr

  character(len=80) :: fname
  integer :: is, s_ierr

  PUSH_SUB(X(restart_read_lr_rho))

  ierr = 0
  do is = 1, nspin
    write(fname, '(a, i1,a)') trim(rho_tag)//'_', is
    call X(restart_read_function)(trim(restart_dir)//trim(restart_subdir), fname, gr%mesh,&
      dl_rho(:, is), s_ierr)
    if( s_ierr /=0 ) ierr = s_ierr
  end do


  if( ierr == 0 ) then 
    write(message(1),'(a)') 'Loaded restart density '//trim(rho_tag)
    call messages_info(1)

  else

    write(message(1),'(a)') 'Could not load restart '//trim(rho_tag)
    call messages_info(1)

  end if

  POP_SUB(X(restart_read_lr_rho))
end subroutine X(restart_read_lr_rho)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
