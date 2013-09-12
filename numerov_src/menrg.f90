module menrg
  !
  !module realted to enrgy evaluation
  !
  ! - kinetic energy T_S
  ! - sum of eigenvalues
  ! - SCE interaction energy
  ! - SCE potential * density integral
  ! - LDA energy by energy density
  ! - LDA potential * density integral
  ! - Hartree energy
  ! - Hartree exchange energy for npart <= 2.0
  ! - external energy
  ! - nuclear repulsion energy
  !
  use mglobvars
  use malloc

  implicit none

  integer, private :: ipt, istate

contains

  subroutine calc_enrg

    write(*,'(a30$)') 'kinetic energy'    
    call calc_ekin
    write(*,*) '.....done'
    write(*,'(a30$)') 'state energy'
    call calc_ewev
    write(*,*) '.....done'

    if ( lsce .or. lscelda .or. lsceldasic) then

       write(*,'(a30$)') 'SCE interaction energy'
       call calc_esceint
       write(*,*) '.....done'
       write(*,'(a30$)') 'v_SCE energy'
       call calc_intvdens( vscetab, escev)
       write(*,*) '.....done'

    end if

    if ( llda .or. lscelda) then

       write(*,'(a30$)') 'LDA interaction energy'
       call calc_eLDAint
       write(*,*) '.....done'
       write(*,'(a30$)') 'v_LDA energy'
       call calc_intvdens( vxctab, eLDAv)
       write(*,*) '.....done'

       if ( llda) then
          write(*,'(a30$)') 'v_Hartree energy'
          call calc_intvdens( vhtab, eHartv)
          eHartv = 0.5d0 * eHartv
          write(*,*) '.....done'
       end if
       
    end if

    if ( lscelda) then

       write(*,'(a30$)') 'SCE+LDA interaction energy'
       call calc_intvdens( eslctab, eSCELDAcint)
       write(*,*) '.....done'
       write(*,'(a30$)') 'v_SCE+LDA energy'
       call calc_intvdens( vslctab, eSCELDAcv)
       write(*,*) '.....done'

    end if

    if ( lsceldasic) then

       write(*,'(a30$)') 'SCE+LDA+SIC interaction energy'
       call calc_intvdens( eslsitab, eSCELDASICint)
       write(*,*) '.....done'
       write(*,'(a30$)') 'v_SCE+LDA+SIC energy'
       call calc_intvdens( vslsitab, eSCELDASICv)
       write(*,*) '.....done'

    end if

    if ( lhf) then

       write(*,'(a30$)') 'Hartree+Fock energy'
       call calc_intvdens( vhtab, eHartv)
       eHartv = ( 1.0d0 - 1.0d0 / npart) * eHartv
       write(*,*) '.....done'

    end if
    
    write(*,'(a30$)') 'external energy'
    call calc_intvdens(vexttab,eext)
    write(*,*) '.....done'

    if(apot) call calc_enuclei

    return

  end subroutine calc_enrg

  subroutine calc_ekin

    integer :: lbnd
    double precision :: wfnder2, stekin, yder2

    ekin = 0.0d0
    lbnd = 0
    if ( hooke .or. cpot) lbnd = 1

    !calculate orbital second derivative and ekin
    do istate = 0, nstates - 1

       stekin = 0.0d0

       if ( loggrid .and. cpot) then

          do ipt = lbnd, npnts

             !second derivative auxilliary function y
             if( xtab(ipt) .le. 0.5d0 * xmax) then
                yder2 = ( wfn(ipt+2,istate) &
                     & - 2.0d0 * wfn(ipt+1,istate) + wfn(ipt,istate)) &
                     & / ( dx * dx)
             else
                yder2 = ( wfn(ipt-2,istate) &
                     & - 2.0d0 * wfn(ipt-1,istate) + wfn(ipt,istate)) &
                     & / ( dx * dx)
             endif

             !second derivative auxilliary function chi
             wfnder2 = ( - 0.25d0 * wfn(ipt,istate) + yder2) &
                  & * xtab(ipt)**(-1.5d0)

             !product of orbital with second derivative chi and summing up
             stekin = stekin - sqrtxtab(ipt) * wfn(ipt,istate) * wfnder2 * &
                  & ( xtab(ipt+1) - xtab(ipt))

          end do

       else

          do ipt = lbnd, npnts

             !second derivative orbital
             if(float(ipt) .le. 0.5d0 * npnts) then
                wfnder2 = ( wfn(ipt+2,istate) &
                     & - 2.0d0 * wfn(ipt+1,istate) + wfn(ipt,istate)) &
                     & / ( dx * dx)
             else
                wfnder2 = ( wfn(ipt-2,istate) &
                     & - 2.0d0 * wfn(ipt-1,istate) + wfn(ipt,istate)) &
                     & / ( dx * dx)
             endif

             !product of orbital with second derivative orbital and summing up
             stekin = stekin - wfn(ipt,istate) * wfnder2 * dx

          end do

       end if

       !sum over states
       ekin = ekin + occns(istate) * stekin * 0.5d0

    end do

    return

  end subroutine calc_ekin

  subroutine calc_ewev
    !sum over state eigenvalues

    ewev = 0.0d0

    do istate = 0, nstates - 1
       ewev = ewev + occns(istate) * energy(istate)
    end do

    return

  end subroutine calc_ewev

  subroutine calc_esceint
    !calculation of SCE energy
    !Integral of density * interaction type

    integer :: ist
    double precision :: howb, vb_int, tbound, rcfac, delt, t, oodelt, &
         & vsoftc_int, vc_int

    esceint = 0.0d0

    if ( lrencoul) then

       tbound = 15.0d0
       howb = 0.5d0 / wireb
       rcfac = sqrt( pi) * howb

       do ipt = 0, npnts

          vb_int = 0.0d0

          do ist = 1, ncomf

             delt = abs( xtab(ipt) - comf(ipt,ist))
             oodelt = 1.0d0 / delt
             t = delt * howb

             if (t .le. tbound) then
                vb_int =  vb_int + rcfac * exp( t * t) * erfc( t)
             else
                vb_int =  vb_int + oodelt -2.0d0 * wireb * wireb * &
                     & oodelt**3 + 12.0d0 * wireb**4 * oodelt**5 &
                     & - 120.0d0 * wireb**6 * oodelt**7
             end if

          end do

          esceint = esceint + density(ipt) * vb_int * dx

       end do

    elseif ( lsoftcoul) then

       do ipt = 0, npnts

          vsoftc_int = 0.0d0

          do ist = 1, ncomf

             t = abs( xtab(ipt) - comf(ipt,ist) )   
             vsoftc_int =  vsoftc_int + (t*t + asoftc*asoftc)**(-0.5d0)

          end do

          esceint = esceint + density(ipt) * vsoftc_int * dx

       end do

    elseif ( lcoul) then

       if ( loggrid) then

          do ipt = 0, npnts

             vc_int = 0.0d0

             do ist = 1, ncomf

                t = abs( xtab(ipt) - comf(ipt,ist) )   
                vc_int =  vc_int + 1.0d0 / t

             end do

             esceint = esceint + density(ipt) * vc_int * ( xtab(ipt+1) - xtab(ipt)) &
                  & * xtab(ipt)

          end do

       else

          do ipt = 0, npnts

             vc_int = 0.0d0

             do ist = 1, ncomf

                t = abs( xtab(ipt) - comf(ipt,ist) )   
                vc_int =  vc_int + 1.0d0 / t

             end do

             esceint = esceint + density(ipt) * vc_int * dx

          end do

       end if

    else

       write(*,*)
       write(*,*) '********************************************'
       write(*,*) '***                                      ***'
       write(*,*) '***  no SCE energy calculation for this  ***'
       write(*,*) '***   type of interaction implemented    ***'
       write(*,*) '***                                      ***'
       write(*,*) '***     -----   terminating   -----      ***'
       write(*,*) '***                                      ***'
       write(*,*) '********************************************'
       write(*,*)
       call dealloc_vars
       stop

    end if
    
    esceint = esceint * 0.5d0

    return

  end subroutine calc_esceint

  subroutine calc_eLDAint

    double precision:: edens(0:npnts)
    !LDA energy computation from enrgy density
    eLDAint = 0.0d0
    edens = extab + ectab

    if ( hpot .or. apot) then

       do ipt = 0, npnts
          eLDAint = eLDAint + density(ipt) * edens(ipt)
       end do

       eLDAint = eLDAint * dx
    
    elseif ( hooke) then

       do ipt = 1, npnts
          eLDAint = eLDAint + density(ipt) * edens(ipt)
       end do

       eLDAint = eLDAint * dx

    elseif ( loggrid) then

       do ipt = 1, npnts
          eLDAint = eLDAint + density(ipt) * edens(ipt) * x2dxtab(ipt)
       end do

    end if

    return

  end subroutine calc_eLDAint

  subroutine calc_intvdens( vtab, enrg)
    !integral of density * function
    !function is usualy a potential

    double precision, intent(in) :: vtab(0:npnts)
    double precision, intent(out) :: enrg

    enrg = 0.0d0

    if ( loggrid) then

       do ipt = 0, npnts
          enrg = enrg + density(ipt) * vtab(ipt) * x2dxtab(ipt)
       end do

    else

       do ipt = 0, npnts
          enrg = enrg + density(ipt) * vtab(ipt) * dx
       end do

    end if

    return

  end subroutine calc_intvdens
  
  subroutine calc_enuclei
  !Energy contribution from the nuclei

    integer :: i, j
   
    enuclei = 0.0d0  

    do i = 1, n_at-1
       do j = 1, n_at-i
          enuclei = enuclei + 1.0d0 / ( SQRT( FLOAT(j) * d * FLOAT(j) * d &
               & + asoftc * asoftc))
       enddo
    enddo

    return
    
  end subroutine calc_enuclei

end module menrg
