# Data files
This folder contains the data files used to support the calculations presented
in the code and the notebook. Please refer to the main README for more
information on the general methods of the code.

## Specific file information
Below we list some of the important metadata for the following data files. Note
that throughout, we use Fortran-style column-major array indexing notation.

### `beam_area.hdf5`
This is a file which contains the beam area information for the HERA beam, which
is necessary for the interferometer noise calculation. The main datasets are:
- `data/freqs`: an array of shape `(Nfreqs,)` that contains the frequencies at
  which the various area values are evaluated, in Hertz.
- `data/omega_p`: an array of shape `(Nfreqs,)` that contains the quantity
  $\Omega_p$ (i.e., the average beam area), in units of rad.
- `data/omega_pp`: an array of shape `(Nfreqs,)` that contains the quantity
  $\Omega_{pp}$ (i.e., the average beam-squared area), in units of rad$^2$.

### `dndz_ares.txt`
This is the average galaxy number density per unit redshift $dn/dz$ taken from
<span style="font-variant:small-caps;">ares</span>. The first column is
redshift, and the second column is the comoving number density, in Mpc$^{-3}$.
Note that all of the values of $dn/dz$ are negative, because the number of
galaxies decreases with increasing redshift (i.e., there are more galaxies at
lower redshift values). Note also that the calculations in `compute_snr.f90` use
little-$h$ units, so these values must be converted to $(\mathrm{Mpc}/h)^{-3}$
when read in.

### `dndz_high.txt`
This is the total number of galaxies in the BlueTides simulation for the assumed
survey parameters of the Roman HLS from [Waters et
al. 2016](https://ui.adsabs.harvard.edu/abs/2016MNRAS.463.3520W/abstract). Note
that despite the name of the file, this is actually $dN/dz$, i.e., the change in
the _total number_ of galaxies. This means that it must be converted to a number
density according to the projected survey volume, as is done in
`compute_snr.f90`.  The first column is redshift, and the second column is
galaxy number. These numbers correspond to the higher edge of the number of
galaxies observed at a 5$\sigma$ detection threshold (the blue band in Figure 3
of the above referenced paper). The calculation in `compute_snr.f90` takes the
geometric mean between this optimistic value and the "pessimistic" values of
`dndz_low.txt`.

### `dndz_low.txt`
The same as above for `dndz_high.txt`, but for the lower edge of the number of
galaxies observed at 5$\sigma$. The calculation in `compute_snr.f90` takes the
geometric mean between this pessimistic value and the "optimistic" values of
`dndz_high.txt`

### `nbl_of_u.txt`
This is the number of baselines in HERA that observe a particular $u$-mode as a
function of $u$. The first column is the value of $u$, and the second column is
the number of baselines $N_\mathrm{bl}(u)$. Note that we can relate $u$ and
$k_\perp$ as $u(\nu) = D_M(\nu) k_\perp / 2\pi$, where $D_M$ is the transverse
comoving distance. The relationship between frequency and redshift is given by
the redshifting of the 21cm line.

### `pk_arrays.hdf5`
This is a file which contains the raw auto- and cross-spectra measured from our
suite of semi-numeric reionization simulations. These values have been averaged
over 30 realizations to produce relatively smooth spectra, which have been
binned as a function of $k_\perp$ and $k_\parallel$. The main datasets are:
- `data/kpara_vals`: a two-dimensional array of size `(Nkpara, Nz)`, where `Nz`
  is the number of redshifts and `Nkpara` is the number of $k_\parallel$
  values. The values of $k_\parallel$ themselves are in units of $h$Mpc$^{-1}$.
- `data/kperp_vals`: a two-dimensional array of size `(Nkperp, Nz)`, where `Nz`
  is the number of redshifts and `Nkperp` is the number of $k_\perp$ values. The
  values of $k_\perp$ themselves are in units of $h$Mpc$^{-1}$.
- `data/pk_gal`: the galaxy auto-power spectrum. The array has size `(3, Nkpara,
  Nkperp, Nz)`, with the same meanings as above. The first column corresponds to
  the number count of how many modes from the raw simulation box went into the
  estimation of this particular mode. The second column is the power averaged
  over those modes. The third column is the empirical standard deviation of the
  power in those modes. Note that this is the power spectrum $P(k_\perp,
  k_\para)$, and as such has units of $(\mathrm{Mpc}/h)^{-3}$.
- `data/pk_t21`: the 21cm auto-power spectrum. The size of the array and
  convention are the same as for the galaxy auto-power spectrum. Note that the
  power spectrum also has "temperature units" associated with it, so that the
  units of this array are $\mathrm{mK}^2 (\mathrm{Mpc}/h)^{-3}$.
- `data/pk_txg`: the 21cm-galaxy cross-spectrum. The size of the array and
  convention are the same as above. The power spectrum has units of $\mathrm{mK}
  (\mathrm{Mpc}/h)^{-3}$.
- `data/zmid_vals`: an array of size `(Nz,)` which contains the central redshift
  value for each of the redshift windows recorded.

### `snr.hdf5`
This is a file which contains the output of the S/N calculation performed in
`compute_snr.f90`. We include a sample output file here for the user's
convenience, in case they do not have access to an HDF5-capable Fortran
compiler. The values in these files can be used to compute the overall S/N
values reported in the paper. The main datasets are:
- `header/area_dsq`: the joint sky coverage assumed for the cross-correlation
  calculation, in deg$^2$.
- `header/bw_MHz`: the bandwidth of 21cm observations assumed, in MHz.
- `header/lae_frac`: an array of size `(Nlae,)` which records the different
  Ly$\alpha$ emitter observation fraction $f_\mathrm{LAE}$ assumed for the S/N
  calculation.
- `header/sigma_z`: an array of size `(Nsigma,)` which records the different
  uncertainties of the observed galaxy redshift $\sigma_z$.
- `header/tobs_s`: the observational time of the 21cm observation, in seconds.
- `header/tsys_K`: the system temperature of the 21cm observation, in Kelvin.
- `data/kpara_vals`: an array of size `(Nkpara, Nz)` where `Nkpara` is the
  number of $k_\parallel$ values observed and `Nz` is the number of redshifts.
  The values of $k_\parallel$ are in $h$Mpc$^{-1}$.
- `data/kperp_vals`: an array of size `(Nkperp, Nz)` where `Nkperp` is the
  number of $k_\perp$ values observed and `Nz` is the number of redshifts.  The
  values of $k_\perp$ are in $h$Mpc$^{-1}$.
- `data/snr`: an array of size `(7, Nkpara, Nkperp, Nsigma, Nlae, Nz)`, where
  the quantities have their same meanings as above. The 7 columns have values
  as follows:
  - 1: the count of values that contributed to the particular bin
  - 2: the cross-power spectrum value $P_{21,\mathrm{gal}}$, units
    $(\mathrm{Mpc}/h)^{-3}$
  - 3: the variance of this quantity $\sigma_{21,\mathrm{gal}}^2$ units
    $(\mathrm{Mpc}/h)^{-6}$
  - 4: the auto-power spectrum value $P_{21}$, units $(\mathrm{Mpc}/h)^{-3}$
  - 5: the variance of this quantity $\sigma_{21}^2$, units
    $(\mathrm{Mpc}/h)^{-6}$
  - 6: the auto-power spectrum value $P_{gal}$, units $(\mathrm{Mpc}/h)^{-3}$
  - 7: the variance of this quantity $\sigma_{gal}^2$, units
    $(\mathrm{Mpc}/h)^{-6}$
- `data/zvals_mid`: an array of size `(Nz,)` which contains the central redshift
  value for each of the redshift windows recorded.
