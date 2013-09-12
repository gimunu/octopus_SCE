module mio
  !I/O related subrouties
  use mglobvars

  implicit none

  integer :: nargs, inarg, dummyits, errvar

  namelist /kssce/ hpot, cpot, omg, L, zeta, readpot, fname_pot, npart,&
       & lrencoul, wireb, maxcycl, updens, ksetol, ksdtol, &
       & readidens, fname_dens, fullout, dflag, ufdens, lcvxc, outplt, &
       & llda, lsce, lksgap, apot, d, n_at, lsoftcoul, asoftc, hooke, lcoul, &
       & lscelda, lhf, lsceldasic
  namelist /numerov/ etol, dtol, maxits, outwrt
  namelist /grid/ xmax, dx, loggrid, xmin

contains

  subroutine input_read
    !reading command line arguments and input file
    logical :: linp, irest
    character(len=32) :: charg, inpname

    linp = .false.
    lrest = .false.
    irest = .false.

    hpot = .false. 
    hooke = .false. 
    omg = 1.0d0
    L = 2.0d0
    cpot = .false.
    zeta = 1.0d0
    lingrid = .true.
    loggrid = .false.
    xmin = -8.0d0
    apot = .false.
    n_at = 1
    d = 1.0d0
    lsce = .false.
    lscelda = .false.
    lsceldasic = .false.
    lcvxc = .false.
    lhf = .false.
    llda = .false.
    fullout = .false.
    ufdens = 0.0d0

    !check for input.inp
    inpname = 'input.inp'
    inquire( file=trim(inpname), exist=linp)

    !process command line arguments
    inarg = 1
    maxits = 0

    do

       call get_command_argument(inarg, charg)

       !reset run
       if ( charg(1:2) == '-r') then

          lrest = .true.
          !check for restart.inp
          inquire( file='restart.inp', exist=irest)

          if ( .not.irest) then
             write(*,*)
             write(*,*) '======= restart.inp does not exist ======='
             write(*,*)
             write(*,*) '==========      terminating      ========='
             write(*,*)
             stop
          end if
          
          !maxits
       elseif ( charg(1:7) == '-maxits') then

          inarg = inarg + 1
          call get_command_argument(inarg, charg)
          read(charg,'(i32)') dummyits

          !input file from command line
       elseif ((charg(1:2) == '-i') .or. (charg(1:2) == '-I')) then

          inarg = inarg + 1
          call get_command_argument(inarg, charg)
          inpname = charg
          linp = .true.

       end if

       inarg = inarg + 1
       if ( inarg  .gt. 10) exit

    end do

    !open and read input file
    if ( (.not.linp) .and. (.not.lrest) ) then
       write(*,*)
       write(*,*) '========= input file does not exist ========='
       write(*,*) 'input file name:', inpname
       write(*,*)
       write(*,*) '===========      terminating      ==========='
       write(*,*)
       stop
    end if

    if ( lrest) then
       open(7, file='restart.inp')
    else
       open(7,file=trim(adjustl(inpname)))
    end if

    rewind(7)
    read(7,nml=kssce,iostat=errvar)
    rewind(7)
    read(7,nml=numerov,iostat=errvar)     
    rewind(7)
    read(7,nml=grid,iostat=errvar)     
    close(7)

    if ( dummyits .gt. 0) maxits = dummyits

    if ( L .ne. 2.0d0) then
       omg = 4 / L**2
    elseif ( omg .ne. 1.0d0) then
       L = 2.0d0 / sqrt( omg)
    end if

    if ( cpot) loggrid = .true.
    if ( loggrid) lingrid = .false.

    !input checks
!    if ( maxits == 0) then
!       maxits = 10
!       write(*,*) 'no maximum number of iterations given'
!       write(*,*) 'assuming maxits = ', maxits
!    end if

    if ( ( hooke .or. cpot) .and. npart .gt. 2.0d0) then

       write(*,*)
       write(*,*) '=========================================='
       write(*,*) '=====                                ====='
       write(*,*) '=====      particle number > 2.0     ====='
       write(*,*) '=====  this is not supported for 3D  ====='
       write(*,*) '=====                                ====='
       write(*,*) '=========================================='
       write(*,*)

       stop

    endif

    return

  end subroutine  input_read

  subroutine out_wrt

    integer :: istate, i
    character(len=32) :: command
    character(len=32) :: name
    character(len=2)  :: cstate

    !result write
    open( 13, file='results.dat')
    call result_wrt( 13)
    close( 13)

    !density write
    name = 'dens_it.dat'
    if ( hpot .or. apot) call write_tab( name, xtab, density)
    if ( hooke) call write_tab( name, xtab, density / (xtab * xtab))
    if ( cpot) call write_tab( name, xtab, density * xtab)

    name = 'dens.dat'
    if ( hpot .or. apot) call write_tab( name, xtab, tdensity)
    if ( hooke) call write_tab( name, xtab, tdensity / (xtab * xtab))
    if ( cpot) call write_tab( name, xtab, tdensity * xtab)

    !wave function write
    if ( fullout) then

       do istate = 0, nstates - 1

          write(cstate,'(i2)') istate
          name = 'wfn_'//trim(adjustl(cstate))//'.dat'

          if ( hooke .or. cpot) then
             !transormation from auxilliary WF to true radial WF taken into account
             call write_tab( name, xtab, wfn(:,istate) / xtab)
          else
             call write_tab( name, xtab, wfn(:,istate))
          endif
          
       end do

    end if
    
    !v_ks write
    name = 'vks.dat'
    call write_tab( name, xtab, pottab)

    !vext write
    name = 'vext.dat'
    call write_tab( name, xtab, vexttab)

    !v_sce write
    if ( lsce .or. lscelda .or. lsceldasic) then
       name = 'vsce.dat'
       call write_tab( name, xtab, vscetab)
!    write( command, * ) 'gnuplot plot_wfn.gn' 
!    call system( trim(command))
!    call sleep(1)
    end if
    
    !v_xc write
    if ( lcvxc) then
       name = 'vxc.dat'
       if ( lsce) call write_tab( name, xtab, vscetab - vhtab)
       if ( lscelda) call write_tab( name, xtab, vscetab + vxctab + vslctab - vhtab)
       if ( lsceldasic) call write_tab( name, xtab, vscetab + sicp * ( vxctab + vslctab) - vhtab)
    end if

    if ( llda) then
       name = 'vxc.dat'
       call write_tab( name, xtab, vxctab)
    end if

    if ( lhf) then
       name = 'vxc.dat'
       call write_tab( name, xtab, vhfxtab)
    end if

    !v_hartree write
    if ( lcvxc .or. llda .or. lhf) then
       name = 'vhartree.dat'
       call write_tab( name, xtab, vhtab)
    end if
    
    !v_corr_scelda write
    if ( lscelda) then
       name = 'vcorr.dat'
       call write_tab( name, xtab, vxctab + vslctab)
    end if

    !v_corr_sceldasic write
    if ( lsceldasic) then
       name = 'vcorr.dat'
       call write_tab( name, xtab, sicp * ( vxctab + vslctab))
    end if

    !comotion function write
    if ( fullout .and. ( lsce .or. lscelda)) then
       do istate = 1, ncomf

          write(cstate,'(i2)') istate
          name = 'comf_'//trim(adjustl(cstate))//'.dat'
          call write_tab( name, xtab, comf(:,istate))

       end do
    end if

    ! write density on full grid
    if ( hpot .or. apot) then

       open( 14, file='dens_final.dat')  
       
       do i = -npnts, -1
          write(14,'(f20.8,f20.8)') -xtab(-i), density(-i) * 0.5d0
       enddo
       
       do i = 0, npnts
          write(14,'(f20.8,f20.8)') xtab(i), density(i) * 0.5d0
       enddo

       close( 14)

    end if
    
    return

  end subroutine out_wrt

  subroutine write_tab( fname, tab1, tab2)

    character(len=32),intent(in) :: fname
    double precision,intent(in) :: tab1(0:npnts), tab2(0:npnts)

    integer :: ipt
    double precision :: tmp1, tmp2

    open( 9, file=trim(fname))
    
    if ( loggrid) then

       do ipt = 0,npnts
          write(9,'(2(es20.12))') tab1(ipt), tab2(ipt)
       end do

    else
       
       do ipt = 0,npnts

          tmp1= tab1(ipt)          
          tmp2= tab2(ipt)          
          !NaN detection
          if ( tmp2 .ne. tmp2) tmp2 = 1.0d-20
          write(9,'(2(es20.12))') tmp1, tmp2

       end do

    end if
    
    close(9)

    return

  end subroutine write_tab
  
  subroutine result_wrt( iunit)
    
    integer, intent(in):: iunit

    integer :: istate

    write(iunit,*)
    write(iunit,*) 'numerov converged    :', nerrflag
    write(iunit,*) 'KSSCE converged      :', lconv
    write(iunit,'(a25,f16.5)') 'density difference   :', rmserr(ccycl)
    write(iunit,'(a25,f16.8)') 'energy difference    :', abs(estates(ccycl-1) - estates(ccycl))
    write(iunit,*)
    write(iunit,*) '=========================================='
    write(iunit,*)
    write(iunit,'(a30,f16.8)') 'Number of electrons :', npart
    write(iunit,'(a30,f16.8)') 'density normalization :', dens_norm

    if ( cpot ) then

       write(iunit,'(a30,f8.4)') 'nuclear charge :', zeta

    elseif ( apot ) then

       write(iunit,'(a30,i3)') 'Number of atomic sites :', int(n_at)
       write(iunit,'(a30,f6.2)') 'atomic-site separation :', d
       write(iunit,'(a30,f6.2)') 'soft-Coulomb parameter :', asoftc

    elseif ( hpot ) then

       write(iunit,'(a30,f10.2)') 'wire effective length :', L
       write(iunit,'(a30,f10.2)') 'wire thickness (b)    :', wireb

    endif

    write(iunit,*)
    write(iunit,*) '=========================================='
    write(iunit,*)

    write(iunit,'(a14)') 'eigenvalues :'
    write(iunit,*)

    do istate = 0, nstates - 1
       write(iunit,'(a18,i5,f20.8)') 'state, energy :', istate, energy(istate)
    enddo

    if( lksgap) then
       write(iunit,*)
       write(iunit,'(a18,f20.8)') 'KS gap :', ks_gap
    endif

    write(iunit,*)
    write(iunit,*) '=========================================='
    write(iunit,*)
    write(iunit,'(a21,f20.12)') 'kinetic energy     :', ekin
    write(iunit,'(a21,f20.12)') 'energy of states   :', ewev

    if ( lsce) then
       write(iunit,'(a21,f20.12)') 'interaction energy :', esceint
       write(iunit,'(a21,f20.12)') 'v_SCE energy       :', escev
    end if

    if ( llda) then
       if (llda) write(iunit,'(a21,f20.12)') 'v_Hartree energy   :', eHartv
       write(iunit,'(a21,f20.12)') 'XC energy          :', eLDAint
       write(iunit,'(a21,f20.12)') 'v_LDA energy       :', eLDAv
    end if

    if ( lscelda) then
       write(iunit,'(a21,f20.12)') 'SCE+LDA energy     :', esceint + eLDAint + eSCELDAcint
       write(iunit,'(a21,f20.12)') 'v_SCE+LDA energy   :', escev + eLDAv + eSCELDAcv
    end if

    if ( lsceldasic) then
       write(iunit,'(a21,f20.12)') 'SCE+LDA+SIC energy :', esceint + eSCELDASICint
       write(iunit,'(a21,f20.12)') 'v_SCE+LDA+SIC energ:', escev + eSCELDASICv
    end if

    write(iunit,'(a21,f20.12)') 'external energy    :', eext
    write(iunit,*) '------------------------------------------'

    if ( lhf) then
       write(iunit,'(a21,f20.12)') 'v_Hartree energy   :', eHartv
       write(iunit,'(a21,f20.12)') 'X energy           :', -0.5d0 * eHartv
    end if

    if ( lsce) then

       write(iunit,'(a21,f20.12)') 'total energy v1    :', ekin + esceint + eext
       write(iunit,*) '(e_kin + e_ext + e_int)'
       write(iunit,'(a21,f20.12)') 'total energy v2    :', - escev + esceint + ewev
       write(iunit,*) '(- e_vsce + e_int + e_states)'
       write(iunit,*) '------------------------------------------'
       write(iunit,'(a21,f20.12)') 'diff e_tot         :', ekin + esceint + eext &
            & - ( - escev + esceint + ewev)

       if(apot) then
          write(iunit,'(a21,f20.8)') 'E_{ee + nuclei} :', ekin + esceint + eext + enuclei
       endif

    elseif ( llda) then

       write(iunit,'(a21,f20.12)') 'total energy v1    :', ekin + eLDAint + eHartv + eext
       write(iunit,*) '(e_kin + e_ext + e_vHart + e_XC)'
       write(iunit,'(a21,f20.12)') 'total energy v2    :', - eHartv + eLDAint - eLDAv + ewev
       write(iunit,*) '(- e_vXC - e_vHart + e_XC + e_states)'
       write(iunit,*) '------------------------------------------'
       write(iunit,'(a21,f20.12)') 'diff e_tot         :', ekin + eLDAint + eHartv + eext &
            & - ( - eHartv + eLDAint - eLDAv + ewev)

       if(apot) then
          write(13,'(a21,f20.8)') 'E_{ee + nuclei} :', ekin + eLDAint + eHartv + eext + enuclei
       endif

    elseif ( lscelda) then

       write(iunit,'(a21,f20.12)') 'total energy v1    :', ekin + esceint + eLDAint + eSCELDAcint + eext
       write(iunit,*) '(e_kin + e_ext + e_int_sce + e_int_XC - e_int_corr)'
       write(iunit,'(a21,f20.12)') 'total energy v2    :', - escev - eLDAv - eSCELDAcv + esceint + eLDAint &
            & + eSCELDAcint+ ewev
       write(iunit,*) '(- e_vsce - e_vXC + e_vcorr + e_int_sce'
       write(iunit,*) '     + e_int_XC - e_int_corr + e_states)'
       write(iunit,*) '------------------------------------------'
       write(iunit,'(a21,f20.12)') 'diff e_tot         :', ekin + esceint + eLDAint + eSCELDAcint + eext &
            & - ( - escev - eLDAv - eSCELDAcv + esceint + eLDAint &
            & + eSCELDAcint + ewev)

       if(apot) then
          write(iunit,'(a21,f20.8)') 'E_{ee + nuclei} :', ekin + esceint + eext + enuclei
       endif

    elseif ( lsceldasic) then

       write(iunit,'(a21,f20.12)') 'total energy v1    :', ekin + esceint + &
            & eSCELDASICint + eext
       write(iunit,*) '(e_kin + e_ext + e_int_sce + e_int_XC - e_int_corr)'
       write(iunit,'(a21,f20.12)') 'total energy v2    :', - escev - eSCELDASICv &
            & + esceint + eSCELDASICint + ewev
       write(iunit,*) '(- e_vsce - e_vXC + e_vcorr + e_int_sce'
       write(iunit,*) '     + e_int_XC - e_int_corr + e_states)'
       write(iunit,*) '------------------------------------------'
       write(iunit,'(a21,f20.12)') 'diff e_tot         :', ekin + esceint &
            & + eSCELDASICint + eext - ( - escev - eSCELDASICv &
            & + esceint + eSCELDASICint + ewev)

       if(apot) then
          write(iunit,'(a21,f20.8)') 'E_{ee + nuclei} :', ekin + esceint + eext + enuclei
       endif

    elseif ( lhf) then

       write(iunit,'(a21,f20.12)') 'total energy v1    :', ekin + 0.5 * eHartv + eext
       write(iunit,*) '(e_kin + e_ext + e_vHart + e_HFx)'
       write(iunit,'(a21,f20.12)') 'total energy v2    :', ewev - 0.5 * eHartv
       write(iunit,*) '( e_HFx + e_states)'
       write(iunit,*) '------------------------------------------'
       write(iunit,'(a21,f20.12)') 'diff e_tot         :', ekin + 0.5d0 * eHartv + eext &
            & - ( ewev - 0.5d0 * eHartv)

    else

       write(iunit,'(a21,f20.12)') 'total energy v1    :', ekin + eext
       write(iunit,*) '(e_kin + e_ext)'
       write(iunit,'(a21,f20.12)') 'total energy v2    :', ewev
       write(iunit,*) '(e_states)'
       write(iunit,*) '------------------------------------------'
       write(iunit,'(a21,f20.12)') 'diff e_tot         :', ekin + eext  - ewev

    end if

    write(iunit,*)
    write(iunit,*) '=========================================='

    return

  end subroutine result_wrt

  subroutine restart_wrt

    !write restart input
    integer :: i

    open( 8, file='restart.inp')
    write(8,nml=KSSCE)
    write(8,nml=numerov)
    write(8,nml=grid)
    close(8)

    return

  end subroutine restart_wrt
  
end module mio
