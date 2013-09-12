module mplot

  use mglobvars

  implicit none

  character(len=32), save :: dlleg1, dlleg2

contains
  
  subroutine plotini

    call page( 2100, 2970)
    call metafl ( 'xwin') 
    call disini
    call pagera
    call complx
    call winmod ( 'none')
    call titlin ( 'progress monitor', 1 )
    call legini ( dlleg1, 2, 16)
    call legini ( dlleg2, 2, 16)
    call legtit ( ' ')
    call axslen ( 1350, 950 )
    call leglin ( dlleg1, 'e_states', 1)
    call leglin ( dlleg1, 'rmserr', 2)
    call leglin ( dlleg2, 'dens.dat', 1)
    call leglin ( dlleg2, 'dens_it.dat', 2)

    return

  end subroutine plotini

  subroutine plot
    
    double precision :: mines, maxes, smines, smaxes, &
         & maxre, smaxre, smaxdens, smaxtdens, smaxd

    write(*,*)
    write(*,*) '========================================='
    write(*,*) '==============             =============='
    write(*,*) '==============   ploting   =============='
    write(*,*) '==============             =============='
    write(*,*) '========================================='
    write(*,*)

    call erase

    !first plot, the energy and density difference 'rmserr'
    call name ( 'step', 'X' )
    call name ( 'energy', 'Y' )
    call labdig ( -1, 'x' )
    call labdig ( 6, 'y' )
    call ticks ( 5, 'y' )
    call axspos ( 400, 1500 )
    
    mines = minval( estates, mask = estates .ne. 0.0d0) 
    maxes = maxval( estates, mask = estates .ne. 0.0d0)  

    if ( hpot .or. hooke) then
       smaxes = maxes * 1.00001d0
       smines = mines * 0.99999d0
    elseif ( apot .or. cpot) then
       smaxes = maxes * 0.99999d0
       smines = mines * 1.00001d0
    endif

    maxre = maxval( abs(rmserr))
    smaxre = maxre * 1.01d0

    call graf ( 1.0d0, dfloat( maxcycl), 1.0d0, &
         & dfloat(maxcycl) / 5.0d0, &
         & smines, smaxes, smines, &
         & ( smaxes - smines) / 5.0d0 )
    call labdig ( 4, 'y' )
    call yaxis( 0.0d0, smaxre, 0.0d0, smaxre/5.0d0, 950, 'rmserr', 0, 1750, 1500)
    call title

    if (  dflag) then
       call color ( 'red')
    else
       call color ( 'magenta')
    end if

    call curve ( cyclarr(1:ccycl), estates(1:ccycl), ccycl)
     
    if (  dflag) then
       call color ( 'orange')
    else
       call color ( 'yellow')
    end if

    call curve ( cyclarr(1:ccycl), rmserr(1:ccycl)*(smaxes-smines)/smaxre+smines, ccycl)  
    call color ( 'white')
    call legend ( dlleg1, 7)
    call endgrf

    !second plot, the density
    call name ( 'r', 'X' )
    call name ( 'density', 'Y' )
    call labdig ( 1, 'x' )
    call labdig ( 3, 'y' )
    call axspos ( 400, 2700 )
    !     mines = minval( density, mask = estates .gt. 0.0d0) * 0.99999d0
    smaxdens = maxval( density) * 1.1d0
    smaxtdens = maxval( tdensity) * 1.1d0
    smaxd = max( smaxdens, smaxtdens)

    call graf ( 0.0d0, xmax, 0.0d0, &
         & xmax / 5.0d0, &
         & 0.0d0, smaxd, 0.0d0, &
         & smaxd / 5.0d0 )
    call color ( 'red')
    call curve ( xtab, tdensity, npnts)  
    call color ( 'orange')
    call curve ( xtab, density, npnts)  
    call color ( 'white')
    call legend ( dlleg2, 7)
    call endgrf
!    call sendbf

    return

  end subroutine plot

  subroutine plotfin

      call disfin

      return

  end subroutine plotfin
  
end module mplot

