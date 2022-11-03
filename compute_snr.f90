program main
  ! OpenMP
  use OMP_LIB

#ifdef HAVE_IFORT
  ! Intel
  use IFPORT
#endif


  ! HDF5
  use ISO_C_BINDING
  use HDF5


  ! Global + tools
  use physical_constants, only: pi,c_cgs,Mpc2cm,rad2deg,f21_cgs
  use cosmo_parameters,   only: OmegaM,OmegaB,hubble0
  use general_tools
  use cosmo_tools


  ! Default
  implicit none


  !! Parameters
  ! Computation and constants
  integer(4), parameter :: N_cpu  = 8

  ! I/O
  character(*), parameter :: dir       = "data/" ! directory of fn_pk and fn_out
  character(*), parameter :: fn_pk     = "pk_arrays.hdf5"
  character(*), parameter :: fn_out    = "snr.hdf5"
  character(*), parameter :: fn_nbl    = "data/nbl_of_u.txt"
  character(*), parameter :: fn_btl    = "data/dndz_low.txt"
  character(*), parameter :: fn_bth    = "data/dndz_high.txt"
  character(*), parameter :: fn_ares   = "data/dndz_ares.txt"
  logical,      parameter :: use_bt    = .false. !true = BlueTides, false = ares
  logical,      parameter :: overwrite = .false.

  ! Survey parameters
  real(8), parameter :: area_dsq    = 500   ! Area of joint-observation
  real(8), parameter :: hls_dsq     = 2200  ! Area of Roman HLS survey
  real(8), parameter :: f_sky       = area_dsq/(4*pi*rad2deg**2)
  real(8), parameter :: hls_sky     = hls_dsq/(4*pi*rad2deg**2)
  real(8), parameter :: tobs_s      = 1D3*3600  ! Observing time in seconds
  real(8), parameter :: tsys_K      = 400       ! System temperature in Kelvin
  real(8), parameter :: bw_MHz      = 6         ! Bandwidth of observation

  ! Define parameter space to search of sigma_z and f_lae
  integer(4), parameter :: N_sigma = 4
  real(8),    parameter, dimension(N_sigma) :: sigma_z  = &
       (/ 0.001D0, 0.01D0, 0.1D0, 0.5D0 /)
  integer(4), parameter :: N_lae   = 3
  real(8),    parameter, dimension(N_lae)   :: lae_frac = &
       (/ 0.01D0, 0.1D0, 1D0 /)


  !! Variables
  integer(4) :: Nperp,Npara,Nredshift,Nkmax,Nprpmax,Npramax


  !! Arrays
  real(8), dimension(:),           allocatable :: zmid_vals
  real(8), dimension(:),           allocatable :: freq_array,omega_p,omega_pp
  real(8), dimension(:,:),         allocatable :: kperp_vals,kpara_vals,p21_var
  real(8), dimension(:,:),         allocatable :: nbl_of_u,bth,btl
  real(8), dimension(:,:),         allocatable :: dndz_a
  real(8), dimension(:,:,:,:),     allocatable :: pk_t21,pk_gal,pk_txg
  real(8), dimension(:,:,:,:,:,:), allocatable :: snr_array_kp


  ! Initialize
  call OMP_SET_NUM_THREADS(N_cpu)
  call calc_cosmo_distance
  call calc_angular_distance

  ! Read in baseline distribution
  call read_nbl_of_u

  ! Read in galaxy distribution
  if (use_bt) then
     call read_bt
  else
     call read_ares
  endif

  ! Read in beam properties
  call read_beam_area
  call calc_hera_variance

  ! Read in spectra
  call read_avg_pk

  ! Compute S/N
  call compute_snr

  ! Write out answer
  call write_snr


contains


  subroutine read_nbl_of_u
    ! Default
    implicit none


    ! Local variables
    integer(4)      :: un,i,nrows
    real(8)         :: u,nbl
    character(1000) :: tmp


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


    ! Scan file
    un = 11
    write(*,*) "Reading ",fn_nbl
    open(un,file=fn_nbl)
    nrows = 0
    read(un,*) tmp  ! throw away header
    do
       read(un,*,end=1) u,nbl
       nrows = nrows + 1
    enddo
1   continue

    ! allocate and initialize
    allocate(nbl_of_u(2,nrows))
    nbl_of_u = 0

    ! read in
    rewind(un)
    read(un,*) tmp
    do i=1,nrows
       read(un,*) u,nbl
       nbl_of_u(1,i) = u
       nbl_of_u(2,i) = nbl
    enddo
    close(un)


    tr2 = omp_get_wtime()
#ifdef HAVE_IFORT
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called read nbl of u'
#else
    ti2 = time()
    call ctime(ti1,ts1)
    call ctime(ti2,ts2)
    write(*,'(f8.2,2a30,a)') tr2-tr1,ts1,ts2,'  Called read nbl of u'
#endif
    return
  end subroutine read_nbl_of_u


!------------------------------------------------------------------------------!


  subroutine read_bt
    ! Default
    implicit none


    ! Local variables
    integer(4) :: un,i,nrows
    real(8)    :: z,dndz


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


    ! Read in data
    un = 12
    write(*,*) "Reading ",fn_btl
    open(un,file=fn_btl)
    nrows = 0
    do
       read(un,*,end=2) z,dndz
       nrows = nrows + 1
    enddo
2   continue
    allocate(btl(2,nrows))
    rewind(un)
    do i=1,nrows
       read(un,*) z,dndz
       btl(1,i) = z
       btl(2,i) = log(abs(dndz))  ! convert to positive dn/dz
    enddo
    close(un)

    un = 13
    write(*,*) "Reading ",fn_bth
    open(un,file=fn_bth)
    nrows = 0
    do
       read(un,*,end=3) z,dndz
       nrows = nrows + 1
    enddo
3   continue
    allocate(bth(2,nrows))
    rewind(un)
    do i=1,nrows
       read(un,*) z,dndz
       bth(1,i) = z
       bth(2,i) = log(abs(dndz))
    enddo
    close(un)


    tr2 = omp_get_wtime()
#ifdef HAVE_IFORT
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called read bt'
#else
    ti2 = time()
    call ctime(ti1,ts1)
    call ctime(ti2,ts2)
    write(*,'(f8.2,2a30,a)') tr2-tr1,ts1,ts2,'  Called read bt'
#endif
    return
  end subroutine read_bt


!------------------------------------------------------------------------------!


  subroutine read_ares
    ! Default
    implicit none


    ! Local variables
    integer(4) :: un,i,nrows
    real(8)    :: z,dndz


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


    ! Read in data
    un = 13
    write(*,*) "Reading ",fn_ares
    open(un,file=fn_ares)
    nrows = 0
    do
       read(un,*,end=4) z,dndz
       nrows = nrows + 1
    enddo
4   continue
    allocate(dndz_a(2,nrows))
    rewind(un)
    do i=1,nrows
       read(un,*) z,dndz
       dndz_a(1,i) = z
       dndz_a(2,i) = log(abs(dndz))
    enddo
    close(un)


    tr2 = omp_get_wtime()
#ifdef HAVE_IFORT
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called read ares'
#else
    ti2 = time()
    call ctime(ti1,ts1)
    call ctime(ti2,ts2)
    write(*,'(f8.2,2a30,a)') tr2-tr1,ts1,ts2,'  Called read ares'
#endif
    return
  end subroutine read_ares


!------------------------------------------------------------------------------!


  subroutine read_beam_area
    ! Default
    implicit none


    ! Local variables
    character(1000)  :: fn
    integer(4)       :: error,Nfreqs
    integer(HID_T)   :: dset_id,file_id,dspace_id
    integer(HSIZE_T) :: dims(1),maxdims(1)


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


    ! Prepare for reading
    fn = "data/beam_area.hdf5"
    write(*,*) "Reading ",trim(fn)
    call h5open_f(error)
    call h5fopen_f(fn, H5F_ACC_RDONLY_F, file_id, error)

    ! Get the size of the datasets
    call h5dopen_f(file_id, "/data/freqs", dset_id, error)
    call h5dget_space_f(dset_id, dspace_id, error)
    call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, error)
    call h5sclose_f(dspace_id, error)
    Nfreqs = dims(1)

    allocate(freq_array(Nfreqs))
    allocate(omega_p(   Nfreqs))
    allocate(omega_pp(  Nfreqs))
    freq_array = 0
    omega_p    = 0
    omega_pp   = 0

    ! Read in data
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, freq_array, dims, error)
    call h5dclose_f(dset_id, error)

    call h5dopen_f(file_id, "/data/omega_p", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, omega_p, dims, error)
    call h5dclose_f(dset_id, error)

    call h5dopen_f(file_id, "/data/omega_pp", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, omega_pp, dims, error)
    call h5dclose_f(dset_id, error)

    ! Close out
    call h5fclose_f(file_id, error)
    call h5close_f(error)


    tr2 = omp_get_wtime()
#ifdef HAVE_IFORT
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called read beam area'
#else
    ti2 = time()
    call ctime(ti1,ts1)
    call ctime(ti2,ts2)
    write(*,'(f8.2,2a30,a)') tr2-tr1,ts1,ts2,'  Called read beam area'
#endif
    return
  end subroutine read_beam_area


!------------------------------------------------------------------------------!


  subroutine read_avg_pk
    ! Default
    implicit none


    ! Local variables
    character(1000)  :: fn
    integer(4)       :: error
    integer(HID_T)   :: file_id,dset_id,dspace_id
    integer(HSIZE_T) :: dims1(1),maxdims1(1),dims2(2),maxdims2(2)
    integer(HSIZE_T) :: dims4(4),maxdims4(4)


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
    call h5open_f(error)

    fn = dir//fn_pk
    write(*,*) "Reading ",trim(fn)

    call h5fopen_f(fn, H5F_ACC_RDONLY_F, file_id, error)

    ! Read in redshift values
    call h5dopen_f(file_id, "/data/zmid_vals", dset_id, error)
    call h5dget_space_f(dset_id, dspace_id, error)
    call h5sget_simple_extent_dims_f(dspace_id, dims1, maxdims1, error)
    call h5sclose_f(dspace_id, error)
    Nredshift = dims1(1)
    allocate(zmid_vals(Nredshift))
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, zmid_vals, dims1, error)
    call h5dclose_f(dset_id, error)

    ! Read in kperp values
    call h5dopen_f(file_id, "/data/kperp_vals", dset_id, error)
    call h5dget_space_f(dset_id, dspace_id, error)
    call h5sget_simple_extent_dims_f(dspace_id, dims2, maxdims2, error)
    call h5sclose_f(dspace_id, error)
    Nperp = dims2(1)
    allocate(kperp_vals(Nperp, Nredshift))
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, kperp_vals, dims2, error)
    call h5dclose_f(dset_id, error)

    ! Read in kpara values
    call h5dopen_f(file_id, "/data/kpara_vals", dset_id, error)
    call h5dget_space_f(dset_id, dspace_id, error)
    call h5sget_simple_extent_dims_f(dspace_id, dims2, maxdims2, error)
    call h5sclose_f(dspace_id, error)
    Npara = dims2(1)
    allocate(kpara_vals(Npara, Nredshift))
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, kpara_vals, dims2, error)
    call h5dclose_f(dset_id, error)

    ! Read in T21
    call h5dopen_f(file_id, "/data/pk_t21", dset_id, error)
    call h5dget_space_f(dset_id, dspace_id, error)
    call h5sget_simple_extent_dims_f(dspace_id, dims4, maxdims4, error)
    call h5sclose_f(dspace_id, error)
    allocate(pk_t21(3,Npara,Nperp,Nredshift))
    allocate(pk_gal(3,Npara,Nperp,Nredshift))
    allocate(pk_txg(3,Npara,Nperp,Nredshift))

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, pk_t21, dims4, error)
    call h5dclose_f(dset_id, error)

    ! Read in galaxy
    call h5dopen_f(file_id, "/data/pk_gal", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, pk_gal, dims4, error)
    call h5dclose_f(dset_id, error)

    ! Read in T21 x galaxy
    call h5dopen_f(file_id, "/data/pk_txg", dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, pk_txg, dims4, error)
    call h5dclose_f(dset_id, error)

    ! Close out
    call h5fclose_f(file_id, error)
    call h5close_f(error)


    tr2 = omp_get_wtime()
#ifdef HAVE_IFORT
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called read avg pk'
#else
    ti2 = time()
    call ctime(ti1,ts1)
    call ctime(ti2,ts2)
    write(*,'(f8.2,2a30,a)') tr2-tr1,ts1,ts2,'  Called read avg pk'
#endif
    return
  end subroutine read_avg_pk


!------------------------------------------------------------------------------!


  pure function calc_ngal(z1,z2)
    ! Default
    implicit none


    ! Function parameters
    integer(4), parameter :: N_z = 10000


    ! Function arguments
    real(8), intent(in) :: z1,z2
    real(8)             :: calc_ngal


    ! Local variables
    integer(4) :: i
    real(8)    :: z,d1,d2,dm,ng,dz


    ! Integrate dN/dz between z1 and z2
    if (use_bt) then
       ! Use geometric mean of high and low bands
       ng = 0
       dz = abs(z2 - z1)/N_z
       do i=1,N_z
          z  = z1 + (i-0.5)*dz
          ! d1 and d2 have had log applied
          d1 = interpolate(z,btl,1,2,'lin')
          d2 = interpolate(z,bth,1,2,'lin')
          dm = (d1+d2)/2  ! taking geometric mean
          dm = exp(dm)    ! exponentiate back to number
          ng = ng + dm
       enddo
       calc_ngal = ng*dz
    else
       ! We have dn/dz directly; need to add h**3 factor
       ng = 0
       dz = abs(z2 - z1)/N_z
       do i=1,N_z
          z  = z1 + (i-0.5)*dz
          d1 = interpolate(z,dndz_a,1,2,'lin')
          ng = ng + exp(d1)
       enddo
       calc_ngal = ng*dz/hubble0**3
    endif

    return
  end function calc_ngal


!------------------------------------------------------------------------------!


  pure function t0_of_z(z)
    ! Default
    implicit none

    ! Function arguments
    real(8), intent(in) :: z
    real(8)             :: t0_of_z


    t0_of_z = 38.6D0*hubble0*(omegab/0.045D0) &
         *sqrt(0.27D0/omegam*(1+z)/10D0)
    return
  end function t0_of_z


!------------------------------------------------------------------------------!


  subroutine calc_hera_variance
    ! Default
    implicit none


    ! Local parameters
    logical, parameter :: write_pvar = .false.


    ! Local variables
    integer(4) :: Nfreqs,ifreq,Npol
    real(8)    :: freq,zval,aval,Hz
    real(8)    :: xval,yval,op,opp


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


    ! Compute the interferometer noise response for each frequency
    Nfreqs = size(freq_array)
    allocate(p21_var(2,Nfreqs))
    p21_var(1,:) = freq_array

    if (write_pvar) then
       open(12,file="p21_var.txt")
       write(12,"(a,a16,4a17,a25)") "#","freq [Hz]","z","a","X [(Mpc/h)/rad]",&
            "Y [(Mpc/h)/Hz]","P21 [mK**2 (Mpc/h)**3]"
    endif

    do ifreq=1,Nfreqs
       freq = freq_array(ifreq)                      ! Hz
       zval = f21_cgs/freq - 1
       aval = 1D0/(1+zval)
       Hz   = 100*E_of_a(aval)                       ! km/s/(Mpc/h)
       xval = interpolate(aval,d_comoving,1,2,'lin') ! (Mpc/h)/rad
       yval = (c_cgs/1D5)*(1+zval)/Hz/freq           ! (Mpc/h)/Hz
       op   = omega_p(ifreq)
       opp  = omega_pp(ifreq)
       Npol = 2  ! xx + yy

       ! Put it all together
       p21_var(2,ifreq) = op**2*xval**2*yval/(opp*Npol)  ! sr*s*(Mpc/h)**3

       if (write_pvar) then
          write(12,"(5es17.8,es25.8)") freq,zval,aval,xval,yval,&
               p21_var(2,ifreq)/tobs_s*tsys_K**2*1D6
       endif
    enddo

    if (write_pvar) then
       close(12)
    endif


    tr2 = omp_get_wtime()
#ifdef HAVE_IFORT
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called calc hera variance'
#else
    ti2 = time()
    call ctime(ti1,ts1)
    call ctime(ti2,ts2)
    write(*,'(f8.2,2a30,a)') tr2-tr1,ts1,ts2,'  Called calc hera variance'
#endif
    return
  end subroutine calc_hera_variance


!------------------------------------------------------------------------------!


  subroutine compute_snr
    ! Default
    implicit none


    ! Local parameters
    logical, parameter :: write_ngal = .false.


    ! Local variables
    integer(4) :: ik,imu,iperp,ipara,iz,Nprp,Npra,isig,ilae
    real(8)    :: kperp,kpara,k,mu,dk,uval,nbl
    real(8)    :: bnorm,zval,aval,t0,d0,f0,f1,f2
    real(8)    :: a1,a2,z1,z2,d1,d2,dd,dA,lam,Hz,w
    real(8)    :: sigma_chi,ngal,vol,sys_noise,sz,lf
    real(8)    :: sigma_a,sigma_b,sigma_c,pt,pg,px,pvar


    ! Local arrays
    real(8), dimension(:,:,:), allocatable :: snr_kp


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


    ! Decide size of output array
    Nkmax   = 0
    Nprpmax = 0
    Npramax = 0
    do iz=1,Nredshift
       Nprp = count(kperp_vals(:,iz) > 0)
       Npra = count(kpara_vals(:,iz) > 0)
       if (any(kperp_vals(:Nprp,iz) == 0)) then
          write(*,*) "invalid kperp ordering; stopping"
          stop
       endif
       if (any(kpara_vals(:Npra,iz) == 0)) then
          write(*,*) "invalid kpara ordering; stopping"
          stop
       endif
       if (any(kperp_vals(2:Nprp,iz) < kperp_vals(1,iz))) then
          write(*,*) "kperp smallest value not at 1st index; stopping"
          stop
       endif
       if (any(kpara_vals(2:Npra,iz) < kpara_vals(1,iz))) then
          write(*,*) "kpara smallest value not at 1st index; stopping"
          stop
       endif

       ! Okay, now we can actually figure out maximum binning
       dk    = min(kperp_vals(1,iz),kpara_vals(1,iz))
       kperp = maxval(kperp_vals(:Nprp,iz))
       kpara = maxval(kpara_vals(:Npra,iz))
       k     = sqrt(kperp**2 + kpara**2)
       Nkmax = max(Nkmax,nint(k/dk))

       Nprpmax = max(Nprpmax,Nprp)
       Npramax = max(Npramax,Npra)
    enddo
    write(*,*) "Nkmax:   ",Nkmax
    write(*,*) "Nprpmax: ",Nprpmax
    write(*,*) "Npramax: ",Npramax
    allocate(snr_array_kp(7,Npramax,Nprpmax,N_sigma,N_lae,Nredshift))
    snr_array_kp = 0

    ! Also allocate per-iteration arrays
    allocate(snr_kp(7,Npramax,Nprpmax))
    snr_kp = 0

    ! Do some up-front work
    bnorm = bw_MHz*1D6*tobs_s

    if (write_ngal) then
       ! Write out file
       if (use_bt) then
          open(12,file="ngal_bt.txt")
       else
          open(12,file="ngal_ares.txt")
       endif
    endif

    ! Loop over redshift
    do iz=1,Nredshift
       ! Get redshift
       zval = zmid_vals(iz)
       aval = 1D0/(1+zval)
       t0   = t0_of_z(zval)

       ! Get depth of observation
       d0  = interpolate(aval,d_comoving,1,2,'lin')
       f0  = f21_cgs*aval
       f1  = f0 - bw_MHz*1D6/2
       f2  = f0 + bw_MHz*1D6/2
       a1  = f1/f21_cgs
       a2  = f2/f21_cgs
       z1  = 1/a1 - 1
       z2  = 1/a2 - 1
       d1  = interpolate(a1,d_comoving,1,2,'lin')
       d2  = interpolate(a2,d_comoving,1,2,'lin')
       dd  = abs(d2-d1)
       dA  = interpolate(aval,angular_distance,1,2,'lin')
       lam = c_cgs/f0

       ! Get observation properties
       pvar = interpolate(f0,p21_var,1,2,'lin')

       ! Get the number of non-zero elements along each axis
       Nprp = count(kperp_vals(:,iz) > 0)
       Npra = count(kpara_vals(:,iz) > 0)
       dk   = min(kperp_vals(1,iz),kpara_vals(1,iz))

       ! Loop over sigma_z and lae_frac parameters
       do isig=1,N_sigma
          sz = sigma_z(isig)
          do ilae=1,N_lae
             lf = lae_frac(ilae)

             ! Get uncertainties for this combo of parameters
             Hz        = 100*E_of_a(aval)   ! h km/s/Mpc
             sigma_chi = (c_cgs/1D5)*sz/Hz  ! Mpc/h
             if (use_bt) then
                ! calc_ngal gives number; convert to number density
                ngal = max(calc_ngal(z1,z2),1D0)*f_sky/hls_sky*lf  ! number
                vol  = f_sky*4*pi*(dA/aval)**2*dd   ! (cMpc/h)**3
                ngal = ngal/vol  ! number density
                if (write_ngal) then
                   if (isig == 1 .and. ilae == 3) then
                      write(12,"(2f7.2,es13.5)") z1,z2,ngal
                   endif
                endif
             else
                ! calc_ngal gives number density
                ngal = calc_ngal(z1,z2)*lf
                if (write_ngal) then
                   if (isig == 1 .and. ilae == 3) then
                      write(12,"(2f7.2,es13.5)") z1,z2,ngal
                   endif
                endif
             endif

             ! Initialize
             snr_kp = 0

             ! Compute the variance for each bin
             !$omp parallel do                    &
             !$omp default(shared)                &
             !$omp private(iperp,ipara,ik,imu)    &
             !$omp private(kperp,kpara,uval,nbl)  &
             !$omp private(k,mu,pt,pg,px,sigma_a) &
             !$omp private(sigma_b,sigma_c)       &
             !$omp private(sys_noise,w)           &
             !$omp reduction(+:snr_kp)
             do iperp=1,Nprp
                kperp = kperp_vals(iperp,iz)

                ! get number of baselines for this kperp
                uval = d0*kperp/(2*pi)
                nbl  = interpolate(uval,nbl_of_u,1,2,'lin')
                if (nbl < 1D0) then
                   ! This is invalid; make nbl very small so the noise term is
                   ! very large
                   nbl = 1d-6
                endif

                ! compute interferometer noise
                sys_noise = tsys_K**2/(t0/1D3)**2 &
                     *pvar/tobs_s/nbl

                do ipara=1,Npra
                   kpara = kpara_vals(ipara,iz)

                   ! Remove factors of T0 from 21cm spectra
                   w   = pk_t21(1,ipara,iperp,iz)
                   pt  = pk_t21(2,ipara,iperp,iz)/t0**2
                   pg  = pk_gal(2,ipara,iperp,iz)
                   px  = pk_txg(2,ipara,iperp,iz)/t0

                   ! Compute contribution to the given bin
                   ! Order of bins:
                   ! count  pA  sigma_A^2  pB  sigma_B^2  pC  sigma_C^2
                   ! where A = 21cm x gal, B = 21cm x 21cm, C = gal x gal
                   sigma_b = (pt + sys_noise)**2
                   sigma_c = (pg + exp((kpara*sigma_chi)**2)/ngal)**2
                   sigma_a = (px**2 + sqrt(sigma_b*sigma_c))/2
                   snr_kp(1,ipara,iperp) = snr_kp(1,ipara,iperp) + w
                   snr_kp(2,ipara,iperp) = snr_kp(2,ipara,iperp) + w*px
                   snr_kp(3,ipara,iperp) = snr_kp(3,ipara,iperp) + w*sigma_a
                   snr_kp(4,ipara,iperp) = snr_kp(4,ipara,iperp) + w*pt
                   snr_kp(5,ipara,iperp) = snr_kp(5,ipara,iperp) + w*sigma_b
                   snr_kp(6,ipara,iperp) = snr_kp(6,ipara,iperp) + w*pg
                   snr_kp(7,ipara,iperp) = snr_kp(7,ipara,iperp) + w*sigma_c
                enddo
             enddo
             !$omp end parallel do

             ! Compute the average for the different bins
             do iperp=1,Nprpmax
                do ipara=1,Npramax
                   if (snr_kp(1,ipara,iperp) > 0) then
                      snr_kp(2:7,ipara,iperp) = snr_kp(2:7,ipara,iperp)&
                           /snr_kp(1,ipara,iperp)
                   else
                      snr_kp(2:7,ipara,iperp) = 0
                   endif
                enddo
             enddo

             ! Save in global array
             snr_array_kp(:,:,:,isig,ilae,iz) = snr_kp
          enddo
       enddo
    enddo

    if (write_ngal) then
       close(12)
    endif


    tr2 = omp_get_wtime()
#ifdef HAVE_IFORT
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called compute snr'
#else
    ti2 = time()
    call ctime(ti1,ts1)
    call ctime(ti2,ts2)
    write(*,'(f8.2,2a30,a)') tr2-tr1,ts1,ts2,'  Called compute snr'
#endif
    return
  end subroutine compute_snr


!------------------------------------------------------------------------------!


  subroutine write_snr
    ! Default
    implicit none


    ! Local variables
    character(1000)  :: fn
    logical          :: avail
    integer(4)       :: error,filter_info
    integer(HID_T)   :: dset_id,file_id,dspace_id
    integer(HID_T)   :: data_id,header_id,dcpl
    integer(HSIZE_T) :: dims1(1),dims2(2),dims6(6),chunk(6)


    ! Local arrays
    real(8), dimension(:,:), allocatable :: kprp,kpra


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


    ! Write out file
    fn = dir//fn_out
    write(*,*) "Writing ",trim(fn)

    call h5open_f(error)
    if (overwrite) then
       call h5fcreate_f(trim(fn), H5F_ACC_TRUNC_F, file_id, error)
    else
       call h5fcreate_f(trim(fn), H5F_ACC_EXCL_F, file_id, error)
       if (error /= 0) then
          write(*,*) "File exists; skipping..."
          ! skip to the end of the routine
          goto 100
       endif
    endif

    ! Create header group
    call h5gcreate_f(file_id, "/header", header_id, error)

    ! Write observational properties
    dims1 = (/ N_sigma /)
    call h5screate_simple_f(1, dims1, dspace_id, error)
    call h5dcreate_f(header_id, "sigma_z", H5T_NATIVE_DOUBLE, dspace_id, &
         dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, sigma_z, dims1, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)

    dims1 = (/ N_lae /)
    call h5screate_simple_f(1, dims1, dspace_id, error)
    call h5dcreate_f(header_id, "lae_frac", H5T_NATIVE_DOUBLE, dspace_id, &
         dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, lae_frac, dims1, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)

    dims1 = (/ 1 /)
    call h5screate_simple_f(1, dims1, dspace_id, error)
    call h5dcreate_f(header_id, "area_dsq", H5T_NATIVE_DOUBLE, dspace_id, &
         dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, area_dsq, dims1, error)
    call h5dclose_f(dset_id, error)

    call h5dcreate_f(header_id, "tobs_s", H5T_NATIVE_DOUBLE, dspace_id, &
         dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tobs_s, dims1, error)
    call h5dclose_f(dset_id, error)

    call h5dcreate_f(header_id, "tsys_K", H5T_NATIVE_DOUBLE, dspace_id, &
         dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, tsys_K, dims1, error)
    call h5dclose_f(dset_id, error)

    call h5dcreate_f(header_id, "bw_MHz", H5T_NATIVE_DOUBLE, dspace_id, &
         dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, bw_MHz, dims1, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)

    call h5gclose_f(header_id, error)

    ! Create data group
    call h5gcreate_f(file_id, "/data", data_id, error)

    ! Write small arrays
    dims1 = (/ Nredshift /)
    call h5screate_simple_f(1, dims1, dspace_id, error)
    call h5dcreate_f(data_id, "zmid_vals", H5T_NATIVE_DOUBLE, dspace_id, &
         dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, zmid_vals, dims1, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)

    allocate(kpra(Npramax,Nredshift))
    allocate(kprp(Nprpmax,Nredshift))
    kpra = kpara_vals(:Npramax,:)
    kprp = kperp_vals(:Nprpmax,:)
    dims2 = (/ Npramax, Nredshift /)
    call h5screate_simple_f(2, dims2, dspace_id, error)
    call h5dcreate_f(data_id, "kpara_vals", H5T_NATIVE_DOUBLE, dspace_id, &
         dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, kpra, dims2, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)

    dims2 = (/ Nprpmax, Nredshift /)
    call h5screate_simple_f(2, dims2, dspace_id, error)
    call h5dcreate_f(data_id, "kperp_vals", H5T_NATIVE_DOUBLE, dspace_id, &
         dset_id, error)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, kprp, dims2, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)

    ! Write snr data
    dims6 = (/ 7, Npramax, Nprpmax, N_sigma, N_lae, Nredshift /)
    call h5screate_simple_f(6, dims6, dspace_id, error)

    ! Check to see if gzip is available
    call h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F, avail, error)
    if (avail) then
       ! Check to see if we can encode
       call h5zget_filter_info_f(H5Z_FILTER_DEFLATE_F, filter_info, error)
       if (iand(filter_info, H5Z_FILTER_ENCODE_ENABLED_F) <= 0) then
          avail = .false.
       endif
    endif

    ! Set dset properties based on filter availablility
    call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, error)
    if (avail) then
       call h5pset_deflate_f(dcpl, 9, error)
       chunk = (/ 7, Npramax, Nprpmax, 1, 1, 1 /)
       call h5pset_chunk_f(dcpl, 6, chunk, error)
    endif

    call h5dcreate_f(data_id, "snr", H5T_NATIVE_DOUBLE, dspace_id, dset_id, &
         error, dcpl)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, snr_array_kp, dims6, error)
    call h5pclose_f(dcpl, error)
    call h5dclose_f(dset_id, error)
    call h5sclose_f(dspace_id, error)


    ! Close out
    call h5gclose_f(data_id, error)
    call h5fclose_f(file_id, error)
    call h5close_f(error)


100 continue
    tr2 = omp_get_wtime()
#ifdef HAVE_IFORT
    call time(ts2)
    write(*,'(f8.2,2a10,a)') tr2-tr1,ts1,ts2,'  Called write snr'
#else
    ti2 = time()
    call ctime(ti1,ts1)
    call ctime(ti2,ts2)
    write(*,'(f8.2,2a30,a)') tr2-tr1,ts1,ts2,'  Called write snr'
#endif
    return
  end subroutine write_snr


end program main
