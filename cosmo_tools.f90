module cosmo_tools
  ! OpenMP
  use OMP_LIB

#ifdef HAVE_IFORT
  ! Intel
  use IFPORT
#endif


  ! Global
  use physical_constants
  use cosmo_parameters


  ! Default
  implicit none


  ! Local parameters
  integer(4), parameter :: N_scalefactor = 1000000


  ! Local arrays
  real(8), dimension(:,:), allocatable :: d_comoving,angular_distance


contains


  pure function E_of_a(a)
    ! Assumes we have a flat universe (OmegaK = 0)
    ! Default
    implicit none


    ! Function arguments
    real(8), intent(in) :: a
    real(8)             :: E_of_a


    ! Calculate Hubble parameter E(a)
    E_of_a = sqrt(OmegaM/a**3 + OmegaR/a**4 + OmegaL/a**(3*(1+wde)))
    return
  end function E_of_a


!------------------------------------------------------------------------------!


  subroutine calc_cosmo_distance
    ! Default
    implicit none


    ! Local variables
    integer(4) :: i
    real(8)    :: a,H,da


    ! Timing variables
#ifdef HAVE_IFORT
    character(8) :: ts1,ts2
    real(8)      :: tr1,tr2
    call time(ts1)
#else
    character(26) :: ts1,ts2
    integer(4)    :: ti1,ti2
    real(8)       :: tr1,tr2
    ti1 = time()
#endif
    tr1 = omp_get_wtime()


    ! Initialize
    allocate(d_comoving(2,N_scalefactor))
    d_comoving = 0


    ! Numerically integrate to calculate time
    ! chi[a] = Integrate[c/(ap)^2/H[ap], {ap, 1, a}]
    da = 1D0/N_scalefactor

    ! Comoving distance in Mpc/h
    d_comoving(1,N_scalefactor) = 1
    d_comoving(2,N_scalefactor) = 0

    do i=N_scalefactor-1,1,-1
       a = (i-0.5)*da
       H = 100*E_of_a(a)

       d_comoving(1,i) = a
       d_comoving(2,i) = d_comoving(2,i+1) + (c_cgs/1D5/H)*da/a**2
    enddo


    tr2 = omp_get_wtime()
#ifdef HAVE_IFORT
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called calc cosmo distance'
#else
    ti2 = time()
    call ctime(ti1,ts1)
    call ctime(ti2,ts2)
    write(*,'(f8.2,2a30,a)') tr2-tr1,ts1,ts2,'  Called calc cosmo distance'
#endif
    return
  end subroutine calc_cosmo_distance


!------------------------------------------------------------------------------!


  subroutine calc_angular_distance
    ! Assumes we have a flat universe (i.e., OmegaK = 0)
    ! Default
    implicit none


    ! Local variables
    integer(4) :: i
    real(8)    :: a,z


    ! Timing variables
#ifdef HAVE_IFORT
    character(8) :: ts1,ts2
    real(8)      :: tr1,tr2
    call time(ts1)
#else
    character(26) :: ts1,ts2
    integer(4)    :: ti1,ti2
    real(8)       :: tr1,tr2
    ti1 = time()
#endif
    tr1 = omp_get_wtime()


    ! Initialize
    allocate(angular_distance(2,N_scalefactor))
    angular_distance = 0


    ! Compute angular diameter distance
    ! For a universe with OmegaK = 0, D_A(z) = D_c(z)/(1+z)
    do i=1,N_scalefactor
       a = d_comoving(1,i)
       z = 1/a - 1

       angular_distance(1,i) = a
       angular_distance(2,i) = d_comoving(2,i)/(1 + z)
    enddo


    tr2 = omp_get_wtime()
#ifdef HAVE_IFORT
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called calc angular distance'
#else
    ti2 = time()
    call ctime(ti1,ts1)
    call ctime(ti2,ts2)
    write(*,'(f8.2,2a30,a)') tr2-tr1,ts1,ts2,'  Called calc angular distance'
#endif
    return
  end subroutine calc_angular_distance


end module cosmo_tools
