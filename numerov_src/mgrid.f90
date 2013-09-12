module mgrid

    use mglobvars

    implicit none

    integer, private :: ipt

contains

  subroutine calc_grid

    double precision :: x

    xtab = 0.0d0
    x2dxtab = 0.0d0
    
    if ( loggrid) then

       xtab(0) = exp(  xmin) / zeta
       sqrtxtab(0)= sqrt( xtab(0))
       x2dxtab(0) = xtab(0)**2.0d0

       do ipt = 1, npnts

          x = xmin + dx * ipt
          xtab(ipt) = exp(  x) / zeta
          sqrtxtab(ipt)= sqrt( xtab(ipt))
          x2dxtab(ipt) = ( xtab(ipt) - xtab(ipt-1)) * xtab(ipt)

       end do

    else

       do ipt=1,npnts
          xtab(ipt) = xtab(ipt-1) + dx
       end do

       sqrtxtab = 1.0d0
       if ( hpot .or. apot) x2dxtab = 2.0d0 * dx
       if ( hooke) x2dxtab = dx

    end if

    write(*,*) 'number grid points:', npnts

    return

  end subroutine calc_grid

end module mgrid
