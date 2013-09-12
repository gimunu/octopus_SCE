module mdensity
  !module related to density
  !- creating first guess for density for KS SCE
  !- calculates density by aufbau principle from states obtained in numerov
  !- calculates rooted mean squared difference between two densities
  !- calculates norm of a density
  use mglobvars

  implicit none

  integer :: ipt

contains

  subroutine init_dens
    !calculates initial density(gaussian) 
    !for initial SCE potential construction

    integer :: iuf
    double precision :: dummy

    density = 0.0d0
    
    if (readidens) then
       
       write(*,*)
       write(*,*) '======   initial density read from file   ======'
       write(*,*)
       write(*,*) '======   make sure the grid is the same   ======'
       write(*,*)
       open( 10, file=fname_dens)

       do ipt = 0, npnts
          read(10,*) dummy, density(ipt)
       end do

       close(10)

       if ( hooke) then 

          density = density * xtab * xtab
          density(0) = 0.0d0
          
       end if       

       if ( cpot) then 

          density = density / xtab
          density(0) = 0.0d0
          
       end if       

    else
       
       if ( ufdens .gt. 0.0d0) then

          iuf = int( ufdens / dx)

!          do ipt = npnts - iuf, 0
!             density(ipt+iuf) = density(ipt)             
!          end do

          do ipt = 0, iuf - 1
             density(ipt) = 1.0d0             
          end do

          do ipt = iuf, npnts
             density(ipt) = exp( -omg * xtab(ipt-iuf) * xtab(ipt-iuf))
          enddo

          if ( hpot .or. apot) density = density / ( ufdens + 0.5d0 * sqrt( pi) / sqrt( omg)) * npart
          if ( hooke .or. cpot) then

             if ( lingrid) density = density * xtab * xtab! * 4.0d0 * pi
             if ( loggrid) density = density * xtab! * 4.0d0 * pi
             call calc_dnorm
             density = density / dens_norm * npart

          end if
          
       else

          do ipt = 0, npnts
             density(ipt) = exp( -omg * xtab(ipt) * xtab(ipt))
          enddo

          if ( hpot .or. apot) density = 2.0d0 * density * npart * sqrt( omg / pi)
          if ( hooke .or. cpot) then 

             density = density * npart / sqrt(pi) * omg**1.5d0 / 0.25d0
             if ( lingrid) density = density * xtab * xtab !* 4.0d0 * pi
             if ( loggrid) density = density * xtab !* 4.0d0 * pi
             
          end if

       end if
    
    end if

    return
    
  end subroutine init_dens

  subroutine calc_tdens

    integer :: istate

    tdensity = 0.0d0

    do istate = 0, nstates-1

       ttdensity = wfn(:,istate) * wfn(:,istate)
       tdensity = tdensity + occns(istate) * ttdensity

    end do
    
    return

  end subroutine calc_tdens

  subroutine calc_densdiff( rmserr)

    double precision,intent(out) :: rmserr

    rmserr = 0.0d0

    do ipt = 0, npnts
       rmserr = rmserr + ( density(ipt) - tdensity(ipt)) &
            & * ( density(ipt) - tdensity(ipt))
    end do

    rmserr = sqrt( rmserr)
    
    return

  end subroutine calc_densdiff

  subroutine calc_dnorm

    dens_norm = 0.0d0

    if ( lingrid) then

       do ipt = 1, npnts
          dens_norm = dens_norm + ( density(ipt-1) + density(ipt)) * x2dxtab(ipt)
       end do
    
       dens_norm = dens_norm * 0.5d0
       if ( hpot .or. apot) dens_norm = dens_norm * 0.5d0

    elseif ( loggrid) then

       do ipt = 1, npnts
          dens_norm = dens_norm + density(ipt) * x2dxtab(ipt)
       end do

    end if
    
    return

  end subroutine calc_dnorm
  
end module mdensity
