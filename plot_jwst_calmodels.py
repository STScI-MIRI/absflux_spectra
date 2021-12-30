from __future__ import absolute_import, print_function, division

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import ScalarFormatter

from astropy.table import Table
import astropy.units as u

from rebin_ndarray import bin_ndarray


FNU = u.erg / (u.cm ** 2 * u.s * u.Hz)
FLAM = u.erg / (u.cm ** 2 * u.s * u.AA)


def set_params(lw=1.5, universal_color="#262626", fontsize=16):
    """Configure some matplotlib rcParams.
    Parameters
    ----------
    lw : scalar
        Linewidth of plot and axis lines. Default is 1.5.
    universal_color : str, matplotlib named color, rgb tuple
        Color of text and axis spines. Default is #262626, off-black
    fontsize : scalar
        Font size in points. Default is 12
    """
    rc("font", size=fontsize)
    rc("lines", linewidth=lw)
    rc("patch", linewidth=lw, edgecolor="#FAFAFA")
    rc(
        "axes",
        linewidth=lw,
        edgecolor=universal_color,
        labelcolor=universal_color,
        axisbelow=True,
    )
    rc("image", origin="lower")
    rc("xtick.major", width=lw)
    rc("xtick.minor", width=lw)
    rc("xtick", color=universal_color)
    rc("ytick.major", width=lw)
    rc("ytick.minor", width=lw)
    rc("ytick", color=universal_color)
    rc("grid", linewidth=lw)
    rc(
        "legend",
        loc="best",
        numpoints=1,
        scatterpoints=1,
        handlelength=1.5,
        fontsize=fontsize,
        columnspacing=1,
        handletextpad=0.75,
    )


def initialize_parser():
    """For running from command line, initialize argparse with common args"""
    ftypes = [
        "png",
        "jpg",
        "jpeg",
        "pdf",
        "ps",
        "eps",
        "rgba",
        "svg",
        "tiff",
        "tif",
        "pgf",
        "svgz",
        "raw",
    ]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s",
        "--savefig",
        action="store",
        default=False,
        choices=ftypes,
        help="Save figure to a file",
    )
    parser.add_argument(
        "-w",
        "--waverange",
        choices=["all", "nir", "mir"],
        default="all",
        help="Wavelength range to display",
    )
    return parser


def trunc_rebin(x, rfac):

    npts = len(x)
    # truncate the the array so it is an integer multiple of rebinned size
    x = x[0 : int(npts / rfac) * rfac]
    return bin_ndarray(x, (int(len(x) / rfac),), operation="mean")


if __name__ == "__main__":

    parser = initialize_parser()
    args = parser.parse_args()

    astar = Table.read("1743045_mod_003.fits")
    gstar = Table.read("p330e_mod_004.fits")
    wdstar = Table.read("gd71_mod_011.fits")
    # wdstar = Table.read('10lac_mod_002.fits')

    #    x = astar['WAVELENGTH']*1e-4
    #    indxs, = np.where((x > 0.6) & (x < 29.))
    #    x = x[indxs]
    #    delta = x[1:] - x[0:-1]
    #    print(x)
    #    print(x[1:]/delta)
    #    print(np.mean(x[1:]/delta))
    #    exit()

    awave = astar["WAVELENGTH"].quantity.value * 1e-4 * u.micron
    aflux = astar["FLUX"].quantity.value * FLAM
    aflux_mJy = aflux.to(u.mJy, u.spectral_density(awave))

    gwave = gstar["WAVELENGTH"].quantity.value * 1e-4 * u.micron
    gflux = gstar["FLUX"].quantity.value * FLAM
    gflux_mJy = gflux.to(u.mJy, u.spectral_density(gwave))

    wdwave = wdstar["WAVELENGTH"].quantity.value * 1e-4 * u.micron
    wdflux = wdstar["FLUX"].quantity.value * FLAM
    wdflux_mJy = wdflux.to(u.mJy, u.spectral_density(wdwave))

    # rebin to a resolution of 3000 (assume input spectrum is R=300,000)
    rfac = 50
    r100fac = 1000
    astar_wave = trunc_rebin(astar["WAVELENGTH"] * 1e-4, rfac)
    astar_flux = trunc_rebin(aflux_mJy, rfac)
    astarr100_wave = trunc_rebin(astar["WAVELENGTH"] * 1e-4, r100fac)
    astarr100_flux = trunc_rebin(aflux_mJy, r100fac)
    gstar_wave = trunc_rebin(gstar["WAVELENGTH"] * 1e-4, rfac)
    gstar_flux = trunc_rebin(gflux_mJy, rfac)
    gstarr100_wave = trunc_rebin(gstar["WAVELENGTH"] * 1e-4, r100fac)
    gstarr100_flux = trunc_rebin(gflux_mJy, r100fac)
    rfac = 10
    r100fac = 200
    wdstar_wave = trunc_rebin(wdstar["WAVELENGTH"] * 1e-4, rfac)
    wdstar_flux = trunc_rebin(wdflux_mJy, rfac)
    wdstarr100_wave = trunc_rebin(wdstar["WAVELENGTH"] * 1e-4, r100fac)
    wdstarr100_flux = trunc_rebin(wdflux_mJy, r100fac)

    fontsize = 18

    set_params(lw=2.0, fontsize=fontsize)

    xsize = 15.0
    ysize = 9.0
    fig, ax = plt.subplots(nrows=3, figsize=(xsize, ysize), sharex=True)

    if args.waverange == "all":
        ptype = "linear"
        kxrange = [0.6, 29.0]
        a_yrange = np.array([0.0, 30.])
        g_yrange = np.array([0.0, 120.])
        wd_yrange = np.array([4., 10.])
    elif args.waverange == "nir":
        ptype = "linear"
        kxrange = [0.6, 5.1]
        a_yrange = np.array([0, 8.0]) * 1e-15
        g_yrange = np.array([0.0, 3.25]) * 1e-14
        wd_yrange = np.array([1.4, 2.75]) * 1e-15
    elif args.waverange == "mir":
        ptype = "linear"
        kxrange = [4.9, 29.0]
        a_yrange = np.array([0.0, 30.])
        g_yrange = np.array([0.0, 120.])
        wd_yrange = np.array([5., 10.])

    astar["WAVELENGTH"] *= 1e-4
    cax = ax[1]
    cax.plot(
        astar["WAVELENGTH"],
        (awave ** 2) * aflux_mJy,
        "k-",
        label="R = 300,000",
        alpha=0.25,
    )
    cax.plot(astar_wave, (astar_wave ** 2) * astar_flux, "b-", label="R ~ 3,000")
    cax.plot(
        astarr100_wave, (astarr100_wave ** 2) * astarr100_flux, "m-", label="R ~ 150"
    )
    cax.set_xscale("log")
    cax.set_xlim(kxrange)
    cax.set_yscale(ptype)
    cax.set_ylim(a_yrange)
    cax.text(0.7, 25., "A dwarf (J1743045)")
    cax.set_ylabel(r"$\lambda^2 F(\nu)$")
    cax.tick_params("both", length=10, width=2, which="major")
    cax.tick_params("both", length=5, width=1, which="minor")
    cax.legend(loc="lower right")

    gstar["WAVELENGTH"] *= 1e-4
    cax = ax[2]
    cax.plot(
        gstar["WAVELENGTH"],
        (gwave ** 2) * gflux_mJy,
        "k-",
        label="R = 300,000",
        alpha=0.25,
    )
    cax.plot(gstar_wave, (gstar_wave ** 2) * gstar_flux, "b-", label="R ~ 3,000")
    cax.set_yscale(ptype)
    cax.set_ylim(g_yrange)
    cax.text(0.65, 85.0, "solar analog (GSPC P330-E)")
    cax.plot(
        gstarr100_wave, (gstarr100_wave ** 2) * gstarr100_flux, "m-", label="R ~ 150"
    )
    cax.set_xlabel(r"wavelength [$\mu m$]")
    cax.set_ylabel(r"$\lambda^2 F(\nu)$")
    cax.tick_params("both", length=10, width=2, which="major")
    cax.tick_params("both", length=5, width=1, which="minor")
    cax.legend(loc="lower right")
    cax.xaxis.set_major_formatter(ScalarFormatter())
    cax.xaxis.set_minor_formatter(ScalarFormatter())
    cax.set_xticks(
        [0.7, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10, 15, 20, 30],
        minor=True,
    )

    wdstar["WAVELENGTH"] *= 1e-4
    cax = ax[0]
    cax.plot(
        wdstar["WAVELENGTH"],
        (wdwave ** 2) * wdflux_mJy,
        "k-",
        label="R = 30,000",
        alpha=0.25,
    )
    cax.plot(wdstar_wave, (wdstar_wave ** 2) * wdstar_flux, "b-", label="R ~ 3,000")
    cax.plot(
        wdstarr100_wave, (wdstarr100_wave ** 2) * wdstarr100_flux, "m-", label="R ~ 150"
    )
    cax.set_yscale(ptype)
    cax.set_ylim(wd_yrange)
    cax.text(0.7, 9.0, "hot star (GD 71)")
    cax.set_ylabel(r"$\lambda^2 F(\nu)$")
    cax.tick_params("both", length=10, width=2, which="major")
    cax.tick_params("both", length=5, width=1, which="minor")
    cax.legend(loc="lower right")

    fig.tight_layout(h_pad=0.15)

    # save the plot
    basename = "jwst_abscal_exspec_" + args.waverange
    if args.savefig:
        fig.savefig("{}.{}".format(basename, args.savefig))
    else:
        plt.show()
