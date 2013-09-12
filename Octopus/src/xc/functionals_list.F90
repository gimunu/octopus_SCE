! Note: this file is generated automatically by build/mk_functionals_list.pl
!
!%Variable XCFunctional
!%Type integer
!%Section Hamiltonian::XC
!%Description
!% Defines the exchange and correlation functional to be used;
!% they should be specified as a sum of a correlation term and an
!% exchange term. Defaults:
!% <br> 1D: lda_x_1d + lda_c_1d_csc
!% <br> 2D: lda_x_2d + lda_c_2d_amgb
!% <br> 3D: lda_x + lda_c_pz_mod
!%Option oep_x                    901
!% OEP: Exact exchange
!%Option ks_inversion             801 
!% Inversion of KS potential
!%Option lda_xc_cmplx             701
!% LDA complex scaled exchange-correlation.
!%Option xc_half_hartree          917
!% Half-Hartree exchange for two electrons (supports complex scaling)
!%Option rdmft_xc_m               601
!% RDMFT Mueller functional
!%Option hgga_hxc_sce_1d           666
!% One dimensional strictly correlated electrons Hartree + XC-functional  
!%Option none                       0
!% Exchange and correlation set to zero.
!%End
