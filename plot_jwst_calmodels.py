from __future__ import absolute_import, print_function, division

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import ScalarFormatter

from astropy.table import Table

from rebin_ndarray import bin_ndarray


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

    astar = Table.read("1743045_mod_002.fits")
    gstar = Table.read("p330e_mod_002.fits")
    wdstar = Table.read("gd71_mod_010.fits")
    # wdstar = Table.read('10lac_mod_002.fits')

    #    x = astar['WAVELENGTH']*1e-4
    #    indxs, = np.where((x > 0.6) & (x < 29.))
    #    x = x[indxs]
    #    delta = x[1:] - x[0:-1]
    #    print(x)
    #    print(x[1:]/delta)
    #    print(np.mean(x[1:]/delta))
    #    exit()

    # rebin to a resolution of 3000 (assume input spectrum is R=300,000)
    rfac = 50
    astar_wave = trunc_rebin(astar["WAVELENGTH"] * 1e-4, rfac)
    astar_flux = trunc_rebin(astar["FLUX"], rfac)
    gstar_wave = trunc_rebin(gstar["WAVELENGTH"] * 1e-4, rfac)
    gstar_flux = trunc_rebin(gstar["FLUX"], rfac)
    rfac = 10
    wdstar_wave = trunc_rebin(wdstar["WAVELENGTH"] * 1e-4, rfac)
    wdstar_flux = trunc_rebin(wdstar["FLUX"], rfac)

    fontsize = 18

    set_params(lw=2.0, fontsize=fontsize)

    xsize = 15.0
    ysize = 9.0
    fig, ax = plt.subplots(nrows=3, figsize=(xsize, ysize), sharex=True)

    if args.waverange == "all":
        ptype = "linear"
        kxrange = [0.6, 29.0]
        a_yrange = np.array([0, 9.0]) * 1e-15
        g_yrange = np.array([0.0, 3.5]) * 1e-14
        wd_yrange = np.array([1.4, 2.75]) * 1e-15
    elif args.waverange == "nir":
        ptype = "linear"
        kxrange = [0.6, 5.1]
        a_yrange = np.array([0, 8.0]) * 1e-15
        g_yrange = np.array([0.0, 3.25]) * 1e-14
        wd_yrange = np.array([1.4, 2.75]) * 1e-15
    elif args.waverange == "mir":
        ptype = "linear"
        kxrange = [4.9, 29.0]
        a_yrange = np.array([6.5, 8.75]) * 1e-15
        g_yrange = np.array([2.25, 3.4]) * 1e-14
        wd_yrange = np.array([2.5, 3.0]) * 1e-15

    astar["WAVELENGTH"] *= 1e-4
    cax = ax[1]
    cax.plot(
        astar["WAVELENGTH"],
        (astar["WAVELENGTH"] ** 4) * astar["FLUX"],
        "k-",
        label="R = 300,000",
        alpha=0.25,
    )
    cax.plot(astar_wave, (astar_wave ** 4) * astar_flux, "b-", label="R ~ 3,000")
    cax.set_xscale("log")
    cax.set_xlim(kxrange)
    cax.set_yscale(ptype)
    cax.set_ylim(a_yrange)
    cax.text(0.7, 7e-15, "A star (1743045)")
    cax.set_ylabel(r"$\lambda^4 F(\lambda)$")
    cax.tick_params("both", length=10, width=2, which="major")
    cax.tick_params("both", length=5, width=1, which="minor")
    cax.legend(loc="lower right")

    gstar["WAVELENGTH"] *= 1e-4
    cax = ax[2]
    cax.plot(
        gstar["WAVELENGTH"],
        (gstar["WAVELENGTH"] ** 4) * gstar["FLUX"],
        "k-",
        label="R = 300,000",
        alpha=0.25,
    )
    cax.plot(gstar_wave, (gstar_wave ** 4) * gstar_flux, "b-", label="R ~ 3,000")
    cax.set_yscale(ptype)
    cax.set_ylim(g_yrange)
    cax.text(0.7, 2.5e-14, "G star (p330e)")
    cax.set_xlabel(r"wavelength [$\mu m$]")
    cax.set_ylabel(r"$\lambda^4 F(\lambda)$")
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
        (wdstar["WAVELENGTH"] ** 4) * wdstar["FLUX"],
        "k-",
        label="R = 30,000",
        alpha=0.25,
    )
    cax.plot(wdstar_wave, (wdstar_wave ** 4) * wdstar_flux, "b-", label="R ~ 3,000")
    cax.set_yscale(ptype)
    cax.set_ylim(wd_yrange)
    cax.text(0.7, 2.5e-15, "WD star (gd71)")
    cax.set_ylabel(r"$\lambda^4 F(\lambda)$")
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
