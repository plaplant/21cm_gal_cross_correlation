module cosmo_parameters
  ! Default
  implicit none


  ! Simulation parameters
  real(8), parameter :: box   = 2000  ! Mpc/h


  ! Cosmological parameters
  ! Planck 2018
  real(8), parameter :: omegam  = 0.315823D0
  real(8), parameter :: omegal  = 1-omegam
  real(8), parameter :: omegab  = 0.049387D0
  real(8), parameter :: hubble0 = 0.673212D0
  real(8), parameter :: sigma8  = 0.8120D0
  real(8), parameter :: nsinit  = 0.96605D0
  real(8), parameter :: wde     = -1.0D0


  ! Baryonic parameters
  real(8), parameter :: YHe         = 0.245401D0
  real(8), parameter :: XHy         = 1-YHe
  real(8), parameter :: Zsolar      = 0.02D0
  real(8), parameter :: gamma_ideal = 5D0/3


  ! CMB parameters
  real(8), parameter :: Tcmb0     = 2.7255D0
  real(8), parameter :: omegar    = 4.48D-7*(1+0.69D0)*Tcmb0**4/hubble0**2
  real(8), parameter :: aequality = omegar/omegam
  real(8), parameter :: zequality = 1D0/aequality - 1


  ! zreion parameters
  real(8), parameter :: b0_zre      = 1D0/1.686D0
  real(8), parameter :: zmean_zre   = 8D0
  real(8), parameter :: alpha_zre   = 0.564D0
  real(8), parameter :: kb_zre      = 0.185D0  ! h/Mpc
  real(8), parameter :: Rsmooth_zre = 1D0      ! Mpc/h


end module cosmo_parameters
