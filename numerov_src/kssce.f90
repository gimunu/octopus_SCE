program skssce
  !*****************************************************************
  !
  ! this program solves self-consistently the Kohn-Sham equations
  !  with the non-interacting potential given by the SCE or LDA
  !
  !*****************************************************************
  use mglobvars
  use mio
  use malloc
  use mgrid
  use mdensity
  use mpot
  use mnumerov
!  use mplot
  use menrg

  implicit none

  integer :: istate, i, j, k
  logical :: errflag

  ksetol = .1d-5
  ksdtol = .1d-1
  rmserr = 0.0d0
  dflag = .false.
  outplt = 10

  write(*,*)
  write(*,*) '========================================'
  write(*,*)
  write(*,*) '=====   welcome to Kohn-Sham SCE   ====='
  write(*,*)
  write(*,*) '========================================'
  write(*,*)

  write(*,'(a30$)') 'reading input'
  call input_read
  write(*,*) '.....done'

  write(*,'(a30)') 'allocating memory'
  call alloc_vars
  write(*,*) '.....done'

  write(*,'(a30)') 'creating grid'
  call calc_grid
  write(*,*) '.....done'

  write(*,'(a30$)') 'initializing denisty'
  call init_dens
  write(*,*) '.....done'
  call calc_dnorm
  write(*,*) 'Initial density normalized to : ', dens_norm


  write(*,'(a35$)') 'calculating external potential'
  call calc_vext
  write(*,*) '.....done'

!!$  write(*,'(a35$)') 'initializing plot window'
!!$  call plotini
!!$!  call plot
!!$  write(*,*) '.....done'

  write(*,'(a40$)') 'starting self-consistent porcedure ...'
  estates = 0.0d0
  rmserr = 0.0d0
  ccycl = 0
  nerrflag = .true.
  lconv = .false.

  forall ( i = 1:maxcycl)
     cyclarr(i) = float(i)
  end forall

  outl: do

!     call plot

     do k = 1, outplt

        ccycl = ccycl + 1
        write(*,*)
        write(*,*) '========================================'
        write(*,*)
        write(*,*) 'KS cycle: ', ccycl
        write(*,*) 'dflag   : ', dflag
        write(*,'(a10,f10.4)') 'rmserr = ', rmserr(ccycl-1)
        write(*,*)

        if ( .not. readpot) then

           write(*,'(a30)') 'calculating potential'
           call calc_pot
           write(*,*) '.....done'

        end if

        emin(:) = minval( pottab(1:npnts), mask = pottab .ne. 0.0d0)
        emax(:) = maxval( pottab(1:npnts), mask = pottab .ne. 0.0d0)

        do istate = 0, nstates - 1

!!$           write(*,'(a30$)') 'adding centrifugal potential'
!!$           if ( cpot) then
!!$
!!$              do i = 0, npnts
!!$                 pottab(i) = pottab(i) + (nl_tab(istate)%l + 0.5d0) * (nl_tab(istate)%l + 0.5d0) / &
!!$                      & ( xtab(i) * xtab(i))
!!$              end do
!!$              
!!$           end if
!!$           write(*,*) '.....done'

           write(*,'(a30,i5$)') 'performing numerov for state', istate
           write(*,*) '.....done'
           nnodes = istate
           call snumerov( dflag, errflag)
           write(*,*) '.....done'
           write(*,'(a30,f20.8)') 'energy of state is:', nenergy
           wfn(:,istate) = nwfn(:)
           energy(istate) = nenergy
           estates(ccycl) = estates(ccycl) + energy(istate)

           !adjust bounds for other states
           do j = istate, nstates - 1
              if ( nenergy .gt. emin(j)) emin(j) = nenergy
           end do

           if (.not. errflag) nerrflag = .false.

        end do

        write(*,'(a30,f20.8)') 'sum state energies is:', estates(ccycl)
        write(*,'(a30$)') 'processing density'
        call calc_tdens
        density = ( 1.0d0 - updens) * density + updens * tdensity
        call calc_densdiff( rmserr(ccycl))
        write(*,*) '.....done'

        if ( .not. nerrflag) then
           write(*,*)
           write(*,*) '=======    something went wrong    ======='
           write(*,*) '===========   during numerov   ==========='
           write(*,*)
           write(*,*) '=======   terminating  Kohn-Sahm   ======='
           write(*,*)
           exit outl
        end if

        if ( abs(estates(ccycl-1) - estates(ccycl)) .lt. ksetol .and. &
             & abs(estates(ccycl-2) - estates(ccycl)) .lt. ksetol) then
           if ( rmserr(ccycl) .lt. ksdtol) then
              if ( dflag) then
                 write(*,*)
                 write(*,*) '======   self-consistency reached   ======'
                 write(*,*)
                 write(*,*) '=======   terminating  Kohn-Sahm   ======='
                 write(*,*)
                 lconv = .true.
                 exit outl
              else
                 dflag = .true.
              end if
           else
              dflag = .true.
           end if

        end if

        if ( ccycl .ge. maxcycl) then
           write(*,*)
           write(*,*) '=========    maximum number of   ========='
           write(*,*) '=========   iterations reached   ========='
           write(*,*)
           write(*,*) '=======   terminating  Kohn-Sahm   ======='
           write(*,*)
           exit outl
        end if

     end do

  end do outl

!  call plot

  !calcuate LUMO and ks_gap for lksgap = .true.
  if ( lksgap) then 

     write(*,'(a30)') 'calculating KS gap'
     nnodes = nstates
     emin(nnodes) = energy(nnodes-1)
     emax(nnodes) = pottab(npnts)
     write(*,'(a30,i5$)') 'performing numerov for LUMO'
     call snumerov( dflag, errflag)
     write(*,*) '.....done'
     write(*,'(a30,f20.8)') 'energy of LUMO is:', nenergy
     ks_gap = nenergy - energy(nnodes-1)
     
  endif

  write(*,'(a30$)') 'writing restart input'
  call restart_wrt
  write(*,*) '.....done'

  if ( lcvxc) then
     write(*,'(a30$)') 'calculating v_Hartree'
     call calc_vhart
     write(*,*) '.....done'
  end if

  write(*,'(a30)') 'calculating energies'
  call calc_enrg
  write(*,*) '.....done'

  write(*,'(a35$)') 'calculating density normalization'
  call calc_dnorm
  write(*,*) '.....done'

  !result output on screen
  call result_wrt( 6)

  write(*,'(a30$)') 'writing output files'
  call out_wrt
  write(*,*) '.....done'

  write(*,'(a30$)') 'deallocating memory'
  call dealloc_vars
  write(*,*) '.....done'

  write(*,*)
  write(*,*) 
  write(*,*) '=======    press "RETURN" to     ======='
  write(*,*) '=======    terminate  program    ======='
  write(*,*)
  read(*,*)

!!$  write(*,'(a30$)') 'terminating dislin'
!!$  call plotfin
!!$  write(*,*) '.....done'

  stop

end program skssce
