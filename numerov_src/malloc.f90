module malloc

  use mglobvars

  implicit none

contains

  subroutine alloc_vars

    integer :: istate, in, il, errm

    frac_tol = 0.000000000001d0

    !determine number of comotion functions
    if ( (npart - int(npart)) .le. frac_tol .or. &
         & int(npart + 1.0d0) - npart .le. frac_tol) then

       write(*,*) 'integer particle number detected'
       npart = float(nint(npart))
       sigma = 1.0d0 

    else

       write(*,*) 'fractional praticle number detected'
       !set sigma
       sigma = npart - float( int( npart - frac_tol))
       write(*,*) 'sigma is', sigma

    end if   

    ncomf = int( npart - frac_tol)

    !deterimne number of states to be searched for
    nstates = int( ( npart + 1.0d0 - frac_tol) * 0.5d0 + 0.5d0)

    !determine number of grid points
    if ( loggrid) then
       npnts = ( log( zeta * xmax) - xmin) / dx
    else
       npnts = xmax / dx
    endif

    allocate( xtab(0:npnts), ftab(0:npnts), nwfn(0:npnts),&
         & energy(0:(nstates-1)), density(0:npnts), &
         & tdensity(0:npnts), ttdensity(0:npnts), &
         & wfn(0:npnts,0:(nstates-1)), occns(0:(nstates-1)), &
         & pottab(0:npnts), vexttab(0:npnts), &
         & sqrtxtab(0:npnts), x2dxtab(0:npnts), &
         & emin(0:nstates-1), emax(0:nstates-1), &
         & estates(maxcycl), rmserr(maxcycl), cyclarr(maxcycl))

    !initialize occupation
    occns = 2.0d0
    occns(nstates-1) = npart - 2.0d0 * float( nstates - 1)

    if ( cpot .or. hooke) then

       allocate( nl_tab(0:nstates-1))
       nl_tab(:)%n = 0
       nl_tab(:)%l = 0
       in = 0
       istate = 0

       out: do
          
          in = in + 1

          do il = 0, in - 1

             nl_tab(istate)%n = in
             nl_tab(istate)%l = il
             istate = istate + 1
             if ( istate .gt. nstates) exit out
             
          end do

       end do out
       
    end if

    if ( loggrid) allocate ( multftab(0:npnts), stat = errm)

    !vHartree
    if ( lcvxc .or. llda .or. lhf) then
       allocate( vhtab(0:npnts))
       if ( hooke .or. cpot) allocate( vh2c(0:npnts))
    end if

    !SCE
    if ( lsce .or. lscelda .or. lsceldasic) then
       allocate( vscetab(0:npnts), ne(0:npnts), sply2tab(1:npnts+1), u(1:npnts), &
            & comf(0:npnts,ncomf), shas(ncomf), shai(ncomf))
    end if
    
    !SCE+LDA
    if ( lscelda) allocate( vslctab(0:npnts), eslctab(0:npnts))

    !SIC+LDA+SIC
    if ( lsceldasic) allocate( vslsitab(0:npnts), eslsitab(0:npnts))
    
    !Hartee-Fock
    if ( lhf) allocate( vhfxtab(0:npnts), ne(0:npnts))

    !LDA
    if ( llda .or. lscelda .or. lsceldasic) then

       allocate( rstab(0:npnts), extab(0:npnts), vxtab(0:npnts), ectab(0:npnts), &
            & vctab(0:npnts), vxctab(0:npnts), densorsq(0:npnts), &
            & partial_extab(0:npnts), partial_ectab(0:npnts))

       if ( lrencoul) then

          !Parameters of the correlation energy functional 
          !(Casula et al., PRB 74, 245427 (2006))
          if (wireb .eq. 0.1d0) then
             xA=4.66d0
             xB=2.092d0
             xC=3.735d0
             xg=1.379d0
             xgp=1.837d0
             xD=23.63d0
             xE=109.9d0
          else
             write(*,*)
             write(*,*) '======================================================='
             write(*,*) '====  renormalized coulomb ee interaction selected ===='
             write(*,*) '====      but no LDA paramters set for this b      ===='
             write(*,*) '====               (check malloc.f90)              ===='
             write(*,*) '====                  terminating                  ===='
             write(*,*) '======================================================='
             call dealloc_vars
             stop
          endif

       elseif ( lsoftcoul) then

          if (asoftc .eq. 1.0d0) then
             xA = 18.40d0
             xB = 0.0d0
             xC = 7.501d0
             xD = 0.10185d0
             xE = 0.012827d0
             xalpha = 1.511d0
             xbeta = 0.258d0
             xm = 4.424
          else
             write(*,*)
             write(*,*) '================================================'
             write(*,*) '====  soft coulomb ee interaction selected  ===='
             write(*,*) '====   but no LDA paramters set for this b  ===='
             write(*,*) '====               (check malloc.f90)       ===='
             write(*,*) '====                  terminating           ===='
             write(*,*) '================================================'
             call dealloc_vars
             stop
          endif

       elseif ( lcoul .and. llda) then

          allocate( ne(0:npnts))

       end if

    end if

    return

  end subroutine alloc_vars
  
  subroutine dealloc_vars

    deallocate(xtab,ftab,nwfn,energy,density,tdensity,ttdensity,&
         & wfn,occns,pottab,sqrtxtab,x2dxtab,emin,emax,&
         & vexttab,estates,rmserr,cyclarr)
    
    if ( cpot .or. hooke) deallocate( nl_tab)
    if ( loggrid) deallocate( multftab)

    !vHartree
    if ( lcvxc .or. llda .or. lhf) then
       deallocate(vhtab)
       if ( hooke .or. cpot) deallocate( vh2c)
    end if

    !SCE
    if ( lsce .or. lscelda .or. lsceldasic) then
       deallocate(vscetab,ne,sply2tab,u,comf,shas,shai)
    end if

    !SCE+LDA
    if ( lscelda) deallocate( vslctab, eslctab)

    !SCE+LDA+SIC
    if ( lsceldasic) deallocate( vslsitab, eslsitab)

    !Hartree-Fock
    if ( lhf) deallocate( vhfxtab, ne)

    !LDA
    if ( llda .or. lscelda .or. lsceldasic) then
       deallocate( rstab, extab, vxtab, ectab, vctab, vxctab, densorsq, partial_extab, partial_ectab)
       if ( lcoul .and. llda) deallocate( ne)
    end if

    return

  end subroutine dealloc_vars

end module malloc


