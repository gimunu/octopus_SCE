module mglobvars

  implicit none
  
  !energies
  double precision :: ekin                  !kinetic energy
  double precision :: eSCEint               !SCE energy by interaction
  double precision :: eSCEv                 !SCE energy by potential
  double precision :: eLDAint               !LDA energy from energy density
  double precision :: eLDAv                 !LDA energy by potential
  double precision :: eHartv                !Hartree energy by potential
  double precision :: eHF                   !Hartree-Fock energy
  double precision :: eSCELDAcv             !SCE+LDA correction energy by potential
  double precision :: eSCELDAcint           !SCE+LDA correction energy by energy dens
  double precision :: eSCELDASICv           !SCE+LDA+SIC energy by potential
  double precision :: eSCELDASICint         !SCE+LDA+SIC energy by energy dens
  double precision :: eext                  !external pot energy
  double precision :: enuclei               !energy contribution from nuclei
  double precision :: ewev                  !sum of wighted eigen values
  double precision :: ks_gap                !HOMO-LUMO gap 
  double precision,allocatable :: cyclarr(:)!number of cycle array
  double precision,allocatable :: estates(:)!sum of state energies
  double precision,allocatable :: rmserr(:) !mean squared error dens

  !electron control
  double precision :: npart                 !particle number
  double precision :: sigma                 !1 - fractional electron number
  double precision :: frac_tol              !smallest sigma allowed

  !Kohn-Sham control
  integer :: outplt                         !intermediate output plot
  double precision :: ufdens                !initial density uniform in core
  integer :: nstates                        !number occupied states
  double precision, allocatable :: occns(:) !occupancy per state

  integer :: maxcycl                        !maximum number KS cycles
  integer :: ccycl                          !counter cycles
  double precision :: updens                !update parameter density
  double precision :: dens_norm             !denisy norm
  double precision :: ksetol                !dflag swith criteria etot
  double precision :: ksdtol                !convergency criteria dens
  logical :: dflag                          !minimizing derivative jump in numerov
  logical :: lconv                          !convergence flag

  double precision, allocatable :: energy(:)         !energy of states
  double precision, allocatable :: emin(:), emax(:)  !energy boundaries
  double precision, allocatable :: wfn(:,:)          !wavefunctions
  double precision, allocatable :: density(:)        !density, in 3D * r^2
  double precision, allocatable :: tdensity(:), ttdensity(:)
  double precision, allocatable :: densorsq(:)       !3D density not multiplied by r^2
  double precision, allocatable :: vexttab(:)        !external potential
  double precision, allocatable :: pottab(:)         !total potential

  !!harmonic external potential
  double precision :: omg                   !omega for harm osc and init dens
  double precision :: L                     !L conversion of omega

  !!1D harmonic
  logical :: hpot                           

  !!3D harmonic
  logical :: hooke                          

  !!3D coulomb external potential
  logical :: cpot                           !coulomb ext potential
  double precision :: zeta                  !nuclear charge

  type tnl_tab
     integer :: n, l
  end type tnl_tab
  
  type (tnl_tab), allocatable :: nl_tab(:)  !table of qunatum numbers

  !!1D soft coulomb external potential
  logical :: apot                           !atomic potential
  double precision :: d                     !separation between atomic sites
  double precision :: n_at                  !number of atoms
  double precision :: z1, z2                !atomic charges (for n_at=2)

  !!3D coulomb ee interaction
  logical :: lcoul                          !flag

  !!1D renormalized coulomb ee interaction
  logical :: lrencoul                       !renormalized coulomb int
  double precision :: wireb                 !wire thicknes

  !!1D soft coulomb ee interaction
  logical :: lsoftcoul                      !soft Coulomb int
  double precision :: asoftc                !coulomb-softening parameter

  logical :: lcvxc                          !v_xc will be calculated
  logical :: lksgap                         !KS gap will be calculated
  logical :: fullout                        !complete output flag

  !!SCE
  logical :: lsce                              !SCE potential flag
  integer :: ncomf                             !number comotion functions
  double precision, allocatable :: vscetab(:)  !SCE potential
  double precision, allocatable :: ne(:)       !the cumulant
  double precision, allocatable :: sply2tab(:) !second derivative spline
  double precision, allocatable :: u(:)        !spline u
  double precision, allocatable :: comf(:,:)   !comotion functions
  double precision, allocatable :: shas(:)     !shell boundary wrt sigma
  double precision, allocatable :: shai(:)     !shell boundary wrt integer

  !!SCE+LDA
  logical :: lscelda                           !SCE potential + LDA correction
  double precision, allocatable :: eslctab(:)  !correction energy density
  double precision, allocatable :: vslctab(:)  !correction potential

  !SCE+LDA+SIC
  logical :: lsceldasic                        !SCE potential + LDA correction
  double precision, allocatable :: eslsitab(:)  !correction energy density
  double precision, allocatable :: vslsitab(:)  !correction potential
                                               ! + self-interaction correction
  double precision :: sicp                     !SI parameter

  !!Hartree
  double precision, allocatable :: vhtab(:)    !Hartree potential
  double precision, allocatable :: vh2c(:)     !Hartree potential 2nd contribution

  !!Hartree-Fock
  logical :: lhf
  double precision, allocatable :: vhfxtab(:)  !HF exchange potential

  !!LDA
  logical :: llda                              !LDA exchange-correlation flag
  double precision :: xA,xB,xC,xg,xgp,xD,xE, & !parameters functional
       & xalpha,xbeta,xm                       !parameters soft coulomb
  double precision, allocatable :: rstab(:)    !wigner seiz radius
  double precision, allocatable :: extab(:)    !exchange energy density
  double precision, allocatable :: vxtab(:)    !exchange potential
  double precision, allocatable :: ectab(:)    !correlation energy density
  double precision, allocatable :: vctab(:)    !correlation potential
  double precision, allocatable :: partial_extab(:)!partial derivative of exchange energy density w.r.t. rho times rho
  double precision, allocatable :: partial_ectab(:)!partial derivative of correlation energy density w.r.t. rho times rho
  double precision, allocatable :: vxctab(:)   !exchange-correlation potential

  !!numerov control
  logical :: nerrflag                      !numerov convergence flag

  integer :: nnodes                        !number of nodes in wfn
  integer :: cnodes                        !node count
  integer :: maxits                        !maximum number iterations
  integer :: outwrt                        !intermediate output write

  double precision :: etol                 !energy tolerance
  double precision :: dtol                 !derivative jump tolerance
  double precision :: nenergy              !energy of state
  double precision :: denergy              !difference energy
  double precision :: tenergy              !temporal variable
  double precision :: wfninwcl             !wfn value of outer wfn 
                                           ! at matching point
  double precision :: scalf                !scaling factor outer wfn
  double precision :: norm                 !norm wfn
  double precision :: djump                !derivative discontinuity
  double precision :: multf                !multiplicator ftab

  double precision, allocatable :: nwfn(:)      !wavefunction
  double precision, allocatable :: ftab(:)      !f-value table
  double precision, allocatable :: multftab(:)  !multiplicator ftab
  double precision, allocatable :: xtab(:)      !grid array
  double precision, allocatable :: sqrtxtab(:)  !square root grid point
  double precision, allocatable :: x2dxtab(:)   !volume element

  !!grid control
  integer :: npnts           !number of grid points
  double precision :: xmax   !outer boundary
  double precision :: dx     !spacing grid points

  logical :: lingrid         !linear grid flag
  logical :: loggrid         !logarithmic grid flag
  double precision :: xmin   !loggrid min val

  !!read potential from external file
  logical :: readpot
  character(len=32) :: fname_pot

  !!read initial density from file
  logical :: readidens       !initial density flag
  character(len=32) :: fname_dens 

  !!restart control
  logical :: lrest           !restart flag

  double precision, parameter :: pi = acos( -1.0d0)

end module

