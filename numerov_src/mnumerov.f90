module mnumerov

  use mglobvars
  use malloc

  implicit none

  integer :: ipt, nanodes, rnodes, ql

contains

  subroutine snumerov( dflag, errflag)

    logical, intent(in) :: dflag
    logical, intent(out) :: errflag

    integer :: icl, incr, cits, iits, j
    double precision :: sqrten, msqrten, monop

    errflag = .false.

    nenergy = ( emax(nnodes) + emin(nnodes)) * 0.5d0
    if ( nenergy > pottab(npnts-1)) nenergy = pottab(npnts-1)
    denergy = 1.0d0
    nwfn = 0.0d0
    if ( cpot .or. hooke) ql = nl_tab(nnodes)%l 

    call calc_multftab
    call calc_ftab

    cits = 0
    incr = outwrt

    !outer loop to write output after N_outwrt iterations
    outwritel: do

       !increment has to be decreased at the end of iterations
       if ((cits + outwrt) > maxits) incr = mod( maxits, outwrt)

       do iits = 1, incr

          !determine matching point
          icl = -1

          lipt: do ipt = npnts - 1, 1, -1

             if ( pottab(ipt) > nenergy) then

                do j = ipt, 1 , -1

                   if ( pottab(j) < nenergy) then
                      icl = j
                      exit lipt
                   endif

                end do

             end if

          end do lipt
          
          if ( icl < 0) then
             write(*,*) 
             write(*,*) '==== no matching point found ===='
             write(*,*)
             write(*,*) '========== terminating =========='
             write(*,*)
             return
          endif

!          write(*,*) 'icl', icl

!!$          !another way of determining the matching point depending on ftab
!!$        do
!!$
!!$           icl = -1
!!$
!!$           do ipt = npnts, 1, -1
!!$              !
!!$              ! beware: if f(i) is exactly zero the change of sign is not observed
!!$              ! the following line is a trick to prevent missing a change of sign 
!!$              ! in this unlikely but not impossible case:
!!$              if( ftab(ipt) == 0.0d0) ftab(ipt) = 1.d-20
!!$
!!$              if( ftab(ipt) /= sign( ftab(ipt), ftab(ipt-1))) then
!!$                 icl = ipt
!!$                 exit
!!$              end if
!!$
!!$           end do
!!$
!!$           write(*,*) 'icl', icl
!!$
!!$           if ( icl < 0 .or. icl >= npnts-2 ) then
!!$              !
!!$              ! classical turning point not found or too far away
!!$              ! no panic: it may follow from a bad choice of eup in
!!$              ! the first iterations. Update energy and emax and retry
!!$              emax(nnodes) = nenergy
!!$!              write(*,*) emax(nnodes) - emin(nnodes)
!!$              if (( emax(nnodes) - emin(nnodes)) .lt. etol) then
!!$                 write(*,*) 
!!$                 write(*,*) '==== no matching point found ===='
!!$                 write(*,*)
!!$                 write(*,*) '========== terminating =========='
!!$                 write(*,*)
!!$                 return
!!$              endif
!!$              nenergy = 0.5d0 * ( emax(nnodes) + emin(nnodes))
!!$              call calc_ftab
!!$              cycle
!!$
!!$           end if
!!$        end do

          !initialize wavefunction
          call init_wfn
          
          !outward integration
          cnodes = rnodes

          outwint: do ipt = 1, icl-1

             nwfn(ipt+1) = (( 12.0d0 - 10.0d0 * ftab(ipt)) * nwfn(ipt) &
                  & - ftab(ipt-1) * nwfn(ipt-1)) / ftab(ipt+1)

             if ( nwfn(ipt) /= sign(nwfn(ipt),nwfn(ipt+1))) then
                cnodes = cnodes + 1
             endif

          enddo outwint

          !adjust bounds of other states
          do j = rnodes, min(nstates-1,(2*cnodes-rnodes-1)), 2
             if ( nenergy .lt. emax(j)) emax(j) = nenergy
          end do

          do j = 2*cnodes-rnodes+1, nstates - 1, 2
             if ( nenergy .gt. emin(j)) emin(j) = nenergy
          end do          

          if ( cnodes == nanodes) then
             !inward integration
             !initial values for outer wfn
             if ( hpot .and. lingrid) then

                !by asymptotic expansion of wavefunction 1D in harmonic potential             
                nwfn(npnts) = xtab(npnts)**nnodes * exp( -omg * 0.5d0 * xtab(npnts) * xtab(npnts))
                nwfn(npnts-1) = xtab(npnts-1)**nnodes * exp( -omg * 0.5d0 * xtab(npnts-1) * xtab(npnts-1))

             elseif ( apot .and. lingrid) then

                !same for the 1D atomic potential
                !only kinetic operator considered
!                nwfn(npnts)   = exp( -sqrt( 2.0d0 * abs( nenergy)) * xtab(npnts))
!                nwfn(npnts-1) = exp( -sqrt( 2.0d0 * abs( nenergy)) * xtab(npnts-1))

                !with coulomb interaction contribtions
                monop = zeta - npart + sigma
                if ( llda) monop = monop - sigma
                sqrten = sqrt( 2.0d0 * abs(nenergy))
                nwfn(npnts)   = xtab(npnts)**( monop / sqrten) * exp( - sqrten * xtab(npnts))
                nwfn(npnts-1) = xtab(npnts-1)**( monop / sqrten) * exp( - sqrten * xtab(npnts-1))

             elseif ( cpot .and. loggrid) then

                !3D coulombic potential
                if( nenergy < 0.0d0) then

                   monop = zeta - npart + sigma
                   if ( llda) monop = monop - sigma
                   sqrten = sqrt( 2.0d0 * abs(nenergy))
                   nwfn(npnts) = xtab(npnts)**( monop / &
                        & sqrten - 0.5d0) * exp( - sqrten * xtab(npnts))
                   nwfn(npnts-1) = xtab(npnts-1)**( monop / &
                        & sqrten - 0.5d0) * exp( - sqrten * xtab(npnts-1))
                   
                else

                   nwfn(npnts) = 0.0d0
                   nwfn(npnts-1) = dx

                end if
                
             elseif ( hooke .and. lingrid) then

                !3D harmonic potential             
                msqrten = nenergy / omg - 1.5d0 + 1.0d0
                nwfn(npnts) = xtab(npnts)**msqrten * &
                     & exp( -omg * 0.5d0 * xtab(npnts) * xtab(npnts))
                nwfn(npnts-1) = xtab(npnts-1)**msqrten * &
                     & exp( -omg * 0.5d0 * xtab(npnts-1) * xtab(npnts-1))

             else

                write(*,*) '***********************************************'
                write(*,*) '*****                                     *****'
                write(*,*) '***** No wavefunction boundary conditions *****'
                write(*,*) '***** implemented for this combination of *****'
                write(*,*) '*****      external potetnial and grid    *****'
                write(*,*) '*****                                     *****'
                write(*,*) '*****             Terminating             *****'
                write(*,*) '*****                                     *****'
                write(*,*) '***********************************************'
                call dealloc_vars
                stop

             endif

             inwint: do ipt = npnts-1, icl+2, -1

                nwfn(ipt-1) = (( 12.0d0 - 10.0d0 * ftab(ipt)) * nwfn(ipt) &
                     & - ftab(ipt+1) * nwfn(ipt+1)) / ftab(ipt-1)

                if (nwfn(ipt-1) > 1.0d10) then
                   do j = npnts, ipt-1, -1
                      nwfn(j) = nwfn(j) / nwfn(ipt-1)
                   end do
                end if

             enddo inwint

             !last step inward integration
             wfninwcl = (( 12.0d0 - 10.0d0 * ftab(icl+1)) * nwfn(icl+1) &
                  & - ftab(icl+2) * nwfn(icl+2)) / ftab(icl)

             !scaling outer function to match
             scalf = nwfn(icl) / wfninwcl

             do ipt = icl + 1, npnts
                nwfn(ipt) = scalf * nwfn(ipt)
             enddo

             !normalize
             norm = 0.0d0

             if ( lingrid) then

                do j = 1, npnts 
                   norm = norm + ( nwfn(j-1) * nwfn(j-1) + nwfn(j) * nwfn(j)) * x2dxtab(j) * 0.5d0
                end do

             elseif ( loggrid) then
                
                do j = 1, npnts
                   norm = norm + nwfn(j)*nwfn(j) * x2dxtab(j)
                end do

             end if
             
             nwfn = nwfn / sqrt( norm)
             if ( hpot .or. apot) nwfn = nwfn * sqrt( 2.0d0)

             !calculate derivative jump
             djump = ( nwfn(icl+1) + nwfn(icl-1) - &
                  & ( 14.0d0 - 12.0d0 * ftab(icl)) * nwfn(icl)) / dx

             if ( denergy < etol) then

                if (dflag) then

                   if ( abs(djump) .le. dtol) then
                      exit outwritel
                   end if

                elseif ( abs(djump) .le. 0.1) then                

!                write(*,*)
!                write(*,*) '=======    convergency reached    ======='
!                write(*,*)
                   exit outwritel

                end if

             end if

             if (djump * nwfn(icl) > 0.0d0) then
                ! Energy is too high --> choose lower energy range
                emax(nnodes) = nenergy
             else
                ! Energy is too low  --> choose upper energy range
                emin(nnodes) = nenergy
             endif

          else
          !narrow the energy range for the other states by the
          !information obtained in this search
             
             if ( cnodes > nanodes) then
!!$             do j = 0, min(nstates-1,(2*cnodes+rnodes-1))
!!$                if ( nenergy .lt. emax(j)) emax(j) = nenergy

!!$             end do
                emax(nnodes) = nenergy
             elseif ( cnodes < nanodes) then
!!$             do j = 2*cnodes+rnodes+1, nstates - 1
!!$                if ( nenergy .gt. emin(j)) emin(j) = nenergy
!!$             end do
                emin(nnodes) = nenergy
             end if
          end if
        
          !
          ! eigenvalue update using perturbation theory
          !
          !           wfncusp = ( wfn(icl-1) * ftab(icl-1) + ftab(icl+1) * wfn(icl+1) & 
          !                & + 10.0d0 * ftab(icl) * wfn(icl)) / 12.0d0
          !           dfcusp = ftab(icl) * ( wfn(icl) / wfncusp - 1.0d0)
          !           de = dfcusp * wfncusp * wfncusp * 12.0d0 / dx 
          !           if (de > 0.0d0) emin = energy
          !           if (de < 0.0d0) emax = energy
          !
          ! prevent e to go out of bounds, i.e. e > eup or e < elw 
          ! (might happen far from convergence)
          !
          !           energy = max( min( energy + de, emax), emin)
          !           call cftab
          !           cycle

          tenergy = ( emax(nnodes) + emin(nnodes)) * .5d0
          if ( tenergy > pottab(npnts-1)) tenergy = pottab(npnts-1)
          denergy = abs( nenergy -  tenergy)
          nenergy = tenergy
!          write(*,*) nenergy, emin(nnodes), emax(nnodes), cnodes, djump

          call calc_ftab

       enddo

       cits = cits + incr

       if (cits .ge. maxits) then
          write(*,*)
          write(*,*) '==== max number iterations reached ===='
          write(*,*)
          write(*,*) '============= terminating ============='
          write(*,*)
          return
       end if

    end do  outwritel

    errflag = .true.

    return

  end subroutine snumerov

  subroutine init_wfn

    nwfn = 0.0d0

    if ( ( hpot .or. apot .or. hooke) .and. lingrid) then
       !initial values inner wavefunction
       if ( mod( nnodes, 2) == 1) then
          !nnodes odd
          nanodes = ( nnodes - 1) / 2 + 1
          rnodes = 1

          nwfn(0) = 0.0d0
          nwfn(1) = dx

       else
          !nnodes even
          nanodes = nnodes / 2
          rnodes = 0

          nwfn(0) = 1.0d0
          nwfn(1) = ( 12.0d0 - 10.0d0 * ftab(0)) / ( 2.0d0 * ftab(1))

       endif

       if (hooke) then

          nwfn(0) = 0.0d0
          nwfn(1) = nwfn(1) * xtab(1)

       end if
       
    elseif( cpot .and. loggrid) then
       !
       ! determination of the wave-function in the first two points 
       ! by using the asymptotic expansion 
       !
       nanodes = nl_tab(nnodes)%n - ql - 1
       rnodes = 0
!       write(*,*) 'ql', ql, nl_tab(nnodes)%n 
       nwfn(0) = xtab(0)**( ql + 0.5d0) * ( 1.0d0 - zeta * xtab(0) / ( ql + 1))
       nwfn(1) = xtab(1)**( ql + 0.5d0) * ( 1.0d0 - zeta * xtab(1) / ( ql + 1))

    else

       write(*,*) '***********************************************'
       write(*,*) '*****                                     *****'
       write(*,*) '***** No wavefunction boundary conditions *****'
       write(*,*) '***** implemented for this combination of *****'
       write(*,*) '*****      external potetnial and grid    *****'
       write(*,*) '*****                                     *****'
       write(*,*) '*****             Terminating             *****'
       write(*,*) '*****                                     *****'
       write(*,*) '***********************************************'
       call dealloc_vars
       stop

    end if

  end subroutine init_wfn

  subroutine calc_multftab

    double precision :: mult

    mult = dx / 6.0d0
    
    if ( loggrid) then

       do ipt = 0, npnts
          multftab(ipt) = x2dxtab(ipt) * mult
       end do

    end if
    
    multf = mult * dx

    return

  end subroutine calc_multftab

  subroutine calc_ftab

    double precision :: lphsq, cfpotf

    if ( cpot .or. hooke) then

       if ( loggrid) then

          lphsq = ( ql + 0.5d0)**2

          do ipt = 0, npnts
             ftab(ipt) = ( nenergy - pottab(ipt)) * multftab(ipt) &
                  & - lphsq * multf * 0.5d0
          enddo
          
       elseif( lingrid) then

          cfpotf = ql * (ql +1) * 0.5d0

          do ipt = 1, npnts
             ftab(ipt) = ( nenergy - pottab(ipt) - cfpotf / ( xtab(ipt) * xtab(ipt))) &
                  & * multf 
          enddo

          if ( hooke) ftab(0) = ftab(1)
          if ( hooke .and. ql == 0) ftab(0) =( nenergy - pottab(ipt)) * multf 

       end if
       
    else

       do ipt = 0, npnts
          ftab(ipt) = ( nenergy - pottab(ipt)) * multf
       enddo

    end if
    
    ftab = ftab + 1.0d0

    return

  end subroutine calc_ftab

end module mnumerov
