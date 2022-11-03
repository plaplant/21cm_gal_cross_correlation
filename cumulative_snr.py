# -*- coding: utf-8 -*-

import os
import h5py
import numpy as np
import pandas as pd
from numba import njit
from astropy import units, constants
from astropy.cosmology import Planck18
from scipy.interpolate import RectBivariateSpline

spectrum_type = "cross"  # must be: "cross", "21cm", or "gal"
cut_redshift = True

fid_dir = "data/"
dir_list = [fid_dir]
dir_labels = ["Fiducial"]

warm_igm_fn = "snr.hdf5"
fn_list = [warm_igm_fn]
fn_labels = ["saturated"]

wedge_cases = ["horizon", 1.0, 0.5]
wedge_labels = ["horizon", "1.0", "0.5"]

# constants
f21_cgs = 1420405751.768
bw_MHz = 6
dlnkperp = 1e-1
dlnkpara = 1e-1
k_max = 1  # h/Mpc

# define wedge slope
def wedge_slope(z):
    dcomoving = Planck18.comoving_distance(z)
    return float(
        dcomoving * Planck18.H(z) / (constants.c * (1 + z))
    )

# define HERA FoV
hera_FoV = 10  # degrees
def lmax_of_z(z):
    """
    Compute the maximum comoving distance for transverse k-modes at redshift z.

    This assumes the field of view of HERA is already defined in degrees. The
    result is in comoving Mpc/h and assumes a Planck18 cosmology.

    Parameters
    ----------
    z : float
        Redshift of interest.

    Returns
    -------
    float:
        The transverse scale corresponding to HERA's FoV in comoving Mpc/h.
    """
    a = 1.0 / (1 + z)
    dA = Planck18.angular_diameter_distance(z).to("Mpc").value
    theta_max = hera_FoV * np.pi / 180.0
    return theta_max * dA / a * Planck18.h

@njit
def calc_nmodes(v_survey, k_perp, k_para, dlnkperp, dlnkpara):
    """
    Compute the number of modes measured in a given Fourier bin.

    Parameters
    ----------
    v_survey : float
        The volume of the survey. Units of comoving length cubed.
    k_perp : float
        The k_perpendicular coordinate. Units of inverse comoving length.
    k_para : float
        The k_parallel coordinate. Units of inverse comoving length.
    dlnkperp : float
        The logarithmic width of the bin in the k_perp direction. Dimensionless.
    dlnkpara : float
        The logarithmic width of the bin in the k_parallel direction.
        Dimensionless.

    Returns
    -------
    float
        The number of Fourier modes that sample this cylindrical bin given the
        proposed survey volume.
    """
    return v_survey / (2 * np.pi)**2 * k_perp**2 * k_para * dlnkperp * dlnkpara

@njit(boundscheck=True)
def nmodes_loop(kperp_vals, kpara_vals, v_survey, nmodes):
    for i in range(kperp_vals.shape[0]):
        kperp = kperp_vals[i]
        for j in range(kpara_vals.shape[0]):
            kpara = kpara_vals[j]

            nmodes[i, j] = calc_nmodes(v_survey, kperp, kpara, dlnkperp, dlnkpara)

    return

@njit(boundscheck=True)
def quick_loop(kperp_vals, kpara_vals, spn, spectrum_type, slp, kmin_hera):
    cs = 0.0
    nkperp = kperp_vals.shape[0]
    nkpara = kpara_vals.shape[0]
    for i in range(nkperp):
        kperp = kperp_vals[i]
        for j in range(nkpara):
            kpara = kpara_vals[j]

            if spectrum_type != "gal":
                # test for wedge cut
                if kpara < kperp * slp:
                    continue
                if kperp < kmin_hera:
                    continue
            cs += spn[i, j]**2

    return cs

# define survey parameters
hls_dsq = 2200
hera_dsq = 1000
cross_dsq = 500
if spectrum_type == "gal":
    area_dsq = hls_dsq
elif spectrum_type == "21cm":
    area_dsq = hera_dsq
elif spectrum_type == "cross":
    area_dsq = cross_dsq

# we assume a square observation for the purposes of mode-counting
if spectrum_type == "gal":
    theta = np.sqrt(area_dsq) * np.pi / 180.0
else:
    area = min(hera_FoV, np.sqrt(area_dsq))
    theta = area * np.pi / 180.0
n_patch_total = max(int(area_dsq / (hera_FoV**2)), 1)
print("theta: ", theta)
print("n_patch: ", n_patch_total)

# initialize
table_rows = []

# loop over histories and IGM properties
for direc, dir_label in zip(dir_list, dir_labels):
    for fn, fn_label in zip(fn_list, fn_labels):
        if dir_label not in ["Fiducial", "Short", "Late"] and fn_label == "cold":
            continue
        if "Deep" in dir_label:
            n_patch = 1.0
        else:
            n_patch = n_patch_total
        # read in data
        snr_fn = os.path.join(direc, fn)
        print(f"reading {snr_fn}...")
        with h5py.File(snr_fn, "r") as h5f:
            zmid_vals = h5f["data/zmid_vals"][()]
            kperp_vals = h5f["data/kperp_vals"][()]
            kpara_vals = h5f["data/kpara_vals"][()]
            snr_array = h5f["data/snr"][()]
            sigma_vals = h5f["header/sigma_z"][()]
            lae_frac = h5f["header/lae_frac"][()]

        # compute the cumulative S/N ratio
        nsigma = sigma_vals.size
        nlae = lae_frac.size

        for isig in range(nsigma):
            sig_val = sigma_vals[isig]
            for ilae in range(nlae):
                f_lae = lae_frac[ilae]
                for slope_val, wedge_label in zip(wedge_cases, wedge_labels):
                    # initialize S/N
                    cs = 0.0
                    for iz, zval in enumerate(zmid_vals):
                        # get observational factors
                        aval = 1.0 / (1 + zval)
                        Lperp_max = (
                            theta
                            * Planck18.angular_diameter_distance(zval).to("Mpc").value
                            / aval
                            * Planck18.h
                        )
                        kperp_min = 2 * np.pi / Lperp_max

                        f0 = f21_cgs * aval
                        f1 = f0 - bw_MHz * 1e6 / 2
                        f2 = f0 + bw_MHz * 1e6 / 2
                        a1 = f1 / f21_cgs
                        a2 = f2 / f21_cgs
                        z1 = 1 / a1 - 1
                        z2 = 1 / a2 - 1
                        d1 = Planck18.comoving_distance(z1).to("Mpc").value
                        d2 = Planck18.comoving_distance(z2).to("Mpc").value
                        Lpara_max = np.abs(d2 - d1) * Planck18.h
                        kpara_min = 2 * np.pi / Lpara_max

                        v_survey = Lperp_max**2 * Lpara_max

                        if cut_redshift:
                            if zval < 7.2:
                                continue
                        if slope_val == "horizon":
                            slp = wedge_slope(zval)
                        else:
                            slp = slope_val
                        hera_lmax = lmax_of_z(zval)
                        kmin_hera = 2 * np.pi / hera_lmax

                        # build splines
                        nxvals = np.count_nonzero(kperp_vals[iz, :])
                        nyvals = np.count_nonzero(kpara_vals[iz, :])
                        if spectrum_type == "cross":
                            signal_array = snr_array[iz, ilae, isig, :nxvals, :nyvals, 1]
                            var_array = snr_array[iz, ilae, isig, :nxvals, :nyvals, 2]
                        elif spectrum_type == "21cm":
                            signal_array = snr_array[iz, ilae, isig, :nxvals, :nyvals, 3]
                            var_array = snr_array[iz, ilae, isig, :nxvals, :nyvals, 4]
                        else:
                            signal_array = snr_array[iz, ilae, isig, :nxvals, :nyvals, 5]
                            var_array = snr_array[iz, ilae, isig, :nxvals, :nyvals, 6]

                        if np.any(np.isinf(var_array)):
                            # replace infinite values
                            varmax = np.amax(
                                var_array, initial=0.0, where=np.isfinite(var_array)
                            )
                            var_array = np.where(
                                np.isfinite(var_array), var_array, varmax
                            )
                        spn_spline = RectBivariateSpline(
                            kperp_vals[iz, :nxvals],
                            kpara_vals[iz, :nyvals],
                            np.log(np.abs(signal_array / np.sqrt(var_array))),
                        )

                        # define kbins given survey area
                        nkperp = int(np.log(k_max / kperp_min) / dlnkperp)
                        nkpara = int(np.log(k_max / kpara_min) / dlnkpara)

                        kperp = np.asarray([
                            np.exp((i + 0.5) * dlnkperp) * kperp_min
                            for i in range(nkperp)
                        ])
                        kpara = np.asarray([
                            np.exp((i + 0.5) * dlnkpara) * kpara_min
                            for i in range(nkpara)
                        ])

                        # do array broadcasting
                        spn = np.exp(spn_spline(kperp, kpara, grid=True))
                        nmodes = np.zeros(
                            (kperp.shape[0], kpara.shape[0]), dtype=np.float64
                        )
                        nmodes_loop(kperp, kpara, v_survey, nmodes)
                        spn *= np.sqrt(nmodes)
                        if spectrum_type != "gal":
                            spn *= np.sqrt(n_patch)

                        loop_output = quick_loop(
                            kperp, kpara, spn, spectrum_type, slp, kmin_hera
                        )

                        cs += loop_output

                    # add to running list of table rows
                    snr = np.sqrt(cs)
                    data = [dir_label, fn_label, f_lae, sig_val, wedge_label, snr]
                    table_rows.append(data)

# make a data frame
columns = ["History", "IGM", "f_LAE", "sigma_z", "wedge", "S/N"]
df = pd.DataFrame(data=table_rows, columns=columns)

# save it out
output_file = f"data/cum_snr_nmodes_{spectrum_type}.hdf5"
print(f"saving {output_file}...")
df.to_hdf(output_file, key="snr", mode="w")

# print to screen
rows = [
    ["Fiducial", "saturated", 0.1, 0.01, "horizon"],
    ["Fiducial", "saturated", 0.1, 0.01, "1.0"],
    ["Fiducial", "saturated", 0.1, 0.01, "0.5"],
    ["Fiducial", "saturated", 0.1, 0.001, "horizon"],
    ["Fiducial", "saturated", 0.1, 0.1, "horizon"],
    ["Fiducial", "saturated", 1, 0.01, "horizon"],
    ["Fiducial", "saturated", 0.01, 0.01, "horizon"],
]
indices = []
for row in rows:
    selection = df[
        (df["History"] == row[0])
        & (df["IGM"] == row[1])
        & (df["f_LAE"] == row[2])
        & (df["sigma_z"] == row[3])
        & (df["wedge"] == row[4])
    ]
    indices.append(selection.index[0])
# select all indices
print(df.iloc[indices])
