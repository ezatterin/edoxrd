"""
Set of functions to plot and compare CTRs and RSMs. Updated 7.11.
"""

import numpy as np
import matplotlib.pyplot as plt
from edoxrd.read import read_rsm_data

def plt_rsm(sub, film, normalisation='sub', scale='log', coordinates='hkl'):

    """
    Plot the RSM of a substrate and film around a reflection on the same
    figure.

    Parameters
    ----------
    sub: string
        The substrate scan, in .asc format.
    film: string
        The film scan, in .asc format.
    d: string
        Directory containing sub and film.
    normalisation: string
        Either 'film' or 'sub'. Normalise the plotted intensity with respect to
        the film or substrate maximum, respecively. Default: 'sub'.
    scale: string
        Either 'log' or 'lin'. Logarithmic or Linear intensity scale for
        the plot, respectively. Default: 'log'.
    coordinates: string
        Either 'hkl' or 'ttomega'. Reciprocal or Real space coordinate mapping
        for the measured intensity.

    Returns
    -------
    Plot. Can be overridden.

    Example
    -------
    >>> plt_rsm('e16014_01_-103_KTO_RSM_2-Theta.asc',
    'e16014_02_-103_PTO_RSM_2-Theta.asc',d)

    """

    # Read data
    h, l, i = read_rsm_data(sub, scale=scale, coordinates=coordinates)
    H, L, I = read_rsm_data(film, scale=scale, coordinates=coordinates)

    # plot
    if normalisation == 'sub':
        plt.pcolormesh(h,l,i,vmin=I.min(),vmax=I.max())
        plt.pcolormesh(H,L,I)
        plt.axis([H[0].min(),H[:,0].max(),l[0].min(),L[:,0].max()])
        plt.axes().set_aspect('equal')
        plt.colorbar()


    elif normalisation == 'film':
        plt.pcolormesh(h,l,i)
        plt.pcolormesh(H,L,I,vmin=i.min(),vmax=i.max())
        plt.axis([H[0].min(),H[:,0].max(),l[0].min(),L[:,0].max()])
        plt.axes().set_aspect('equal')
        plt.colorbar()


    if coordinates == 'ttomega':
        plt.xlabel(r'$/omega$'); plt.ylabel(r'$2/theta')
    elif coordinates == 'hkl':
        plt.xlabel('H'); plt.ylabel('L')

def plt_prof(sample, H, L, scale='log'):

    """
    Plots profile of an RSM plot. TODO!
    """

    h,l,I = read_rsm_data(sample, scale=scale)

    a = abs(l-L)
    b = abs(h-H)

    mask_a = ((a > 0) & (a < 1e-4))
    mask_b = ((b > 0) & (b < 1e-4))

    lst_a, lst_b = [], []

    for number in a[mask_a]:
        lst_a.append([(int(np.where(a==number)[0])),int(np.where(a==number)[1])])
    for number in b[mask_b]:
        lst_b.append([(int(np.where(b==number)[0])),int(np.where(b==number)[1])])

    idx_l = np.array(lst_a)
    idx_h = np.array(lst_b)

    I_h = I[idx_l[:,0], idx_l[:,1]][::-1]
    I_l = I[idx_h[:,0], idx_h[:,1]][::-1]

    h_h = h[idx_l[:,0], idx_l[:,1]][::-1]
    l_l = l[idx_h[:,0], idx_h[:,1]][::-1]

    fig, ax = plt.subplots(1,3, figsize=(13,5))
    ax[0].pcolormesh(h,l,I); ax[0].set_ylabel('L'); ax[0].set_xlabel('H')
    ax[0].plot([h.min(),h.max()],[L,L],'r-')
    ax[0].plot([H, H], [l.min(),l.max()],'r-')
    ax[1].plot(h_h, I_h); ax[1].set_ylabel(r'$I$'); ax[1].set_xlabel('H')
    ax[2].plot(l_l, I_l); ax[2].set_ylabel(r'$I$'); ax[2].set_xlabel('L')
    plt.tight_layout()
    plt.show()
