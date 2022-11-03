module physical_constants
  ! Default
  implicit none
  public
  save


  ! Mathematical constants
  real(8), parameter :: pi = acos(-1D0)


  ! Conversions
  real(8), parameter :: yr2sec  = 365*86400
  real(8), parameter :: sec2yr  = 1D0/yr2sec

  real(8), parameter :: Myr2sec = yr2sec*1D6
  real(8), parameter :: sec2Myr = sec2yr/1D6

  real(8), parameter :: hr2rad = pi/12
  real(8), parameter :: rad2hr = 1D0/hr2rad

  real(8), parameter :: deg2rad    = pi/180
  real(8), parameter :: rad2deg    = 1D0/deg2rad

  real(8), parameter :: arcmin2rad = deg2rad/60
  real(8), parameter :: rad2arcmin = 1D0/arcmin2rad


!------------------------------------------------------------------------------!


  ! Physical constants
  real(8), parameter :: c_cgs   = 2.99792D+10
  real(8), parameter :: c_mks   = 2.99792D+08
  real(8), parameter :: e_cgs   = 4.80321D-10
  real(8), parameter :: e_mks   = 1.60218D-19
  real(8), parameter :: eV_cgs  = 1.60218D-12
  real(8), parameter :: eV_mks  = 1.60218D-19
  real(8), parameter :: G_cgs   = 6.67259D-08
  real(8), parameter :: G_mks   = 6.67259D-11
  real(8), parameter :: h_cgs   = 6.62608D-27
  real(8), parameter :: h_mks   = 6.62608D-34
  real(8), parameter :: k_cgs   = 1.38066D-16
  real(8), parameter :: k_mks   = 1.38066D-23
  real(8), parameter :: mu_cgs  = 1.66054D-24
  real(8), parameter :: mu_mks  = 1.66054D-27
  real(8), parameter :: mH_cgs  = 1.67223D-24
  real(8), parameter :: mH_mks  = 1.67223D-27
  real(8), parameter :: me_cgs  = 9.10939D-28
  real(8), parameter :: me_mks  = 9.10939D-31
  real(8), parameter :: mHe_cgs = 4*mH_cgs
  real(8), parameter :: mHe_mks = 4*mH_mks
  real(8), parameter :: sT_cgs  = 6.65246D-25
  real(8), parameter :: sT_mks  = 6.65246D-29
  real(8), parameter :: sS_cgs  = 5.67051D-05
  real(8), parameter :: sS_mks  = 5.67051D-08


  ! Conversions
  real(8), parameter :: K2erg = k_cgs
  real(8), parameter :: erg2K = 1D0/K2erg

  real(8), parameter :: K2eV  = k_cgs/eV_cgs
  real(8), parameter :: eV2K  = 1D0/K2eV
  real(8), parameter :: K2keV = K2eV/1E3
  real(8), parameter :: keV2K = eV2K*1E3

  real(8), parameter :: eV2Hz = eV_cgs/h_cgs
  real(8), parameter :: Hz2eV = 1D0/eV2Hz

  real(8), parameter :: eV2erg = eV_cgs
  real(8), parameter :: erg2eV = 1D0/eV2erg


!------------------------------------------------------------------------------!


  ! Astrophysical constants
  real(8), parameter :: Msun_cgs = 1.98900D+33
  real(8), parameter :: Msun_mks = 1.98900D+30
  real(8), parameter :: Rsun_cgs = 6.96000D+10
  real(8), parameter :: Rsun_mks = 6.96000D+08
  real(8), parameter :: Lsun_cgs = 3.82600D+33
  real(8), parameter :: Lsun_mks = 3.82600D+26


  ! Cosmological constants
  real(8), parameter :: H0_cgs   = 3.24086D-18
  real(8), parameter :: H0_mks   = 3.24086D-18
  real(8), parameter :: rhoc_cgs = 1.87890D-29
  real(8), parameter :: rhoc_mks = 1.87890D-26
  real(8), parameter :: rhoc_ast = 2.77550D+11
  real(8), parameter :: Mpc2cm   = 3.08560D+24
  real(8), parameter :: Mpc2m    = 3.08560D+22
  real(8), parameter :: Mpc2km   = 3.08560D+19
  real(8), parameter :: cm2Mpc   = 1D0/Mpc2cm
  real(8), parameter :: m2Mpc    = 1D0/Mpc2m
  real(8), parameter :: km2Mpc   = 1D0/Mpc2km


!------------------------------------------------------------------------------!


  ! Physical constants in astrophysical units
  ! G_ast has units of Msun*Mpc^3/Gyr^2
  ! Use to get potential in (Mpc/Gyr)^2
  real(8), parameter :: G_ast = G_cgs &
                              * cm2Mpc**3*Msun_cgs*(1D9/sec2yr)**2


!------------------------------------------------------------------------------!


  ! Atomic wavelengths/frequencies
  real(8), parameter :: ang2cm    = 1D-08
  real(8), parameter :: ang2m     = 1D-10
  real(8), parameter :: lya_cgs   = 1.21567D-05
  real(8), parameter :: lya_mks   = 1.21567D-07
  real(8), parameter :: HI_eV     = 13.60
  real(8), parameter :: HeI_eV    = 24.59
  real(8), parameter :: HeII_eV   = 54.42
  real(8), parameter :: HI_K      = HI_eV  *eV2K
  real(8), parameter :: HeI_K     = HeI_eV *eV2K
  real(8), parameter :: HeII_K    = HeII_eV*eV2K
  real(8), parameter :: lam21_cgs = 21.106D0
  real(8), parameter :: f21_cgs   = 1420405751.768D0


end module physical_constants
