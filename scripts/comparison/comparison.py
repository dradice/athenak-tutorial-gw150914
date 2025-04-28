"""
Align and plot the SXS and ATK waveforms targeted to
GW150914.

If the SXS waveform is not already in the current directory, it will be downloaded
from the SXS catalog.

The ATK waveform is read from a local file.

run with:
python atksxs_305.py --t_window 600 1900 --path_atk ./ --path_sxs ./
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import h5py

try:
    from PyART.catalogs import sxs
    from PyART.waveform import Waveform
    from PyART.utils import utils, wf_utils
except ImportError:
    raise ImportError("PyART is not installed. Please install from \
                      git@github.com:RoxGamba/PyART.git \
                      and try again."
                      )
import argparse

matplotlib.rc('text', usetex=True)

def align_eobnr(wf1, wf2, time_window=[1000, 2500], ref=None):
    """
    Align two waveforms based on the (2,2) mode.
    """
    # align the (2,2) mode
    h22_1   = wf1.hlm[(2,2)]
    h22_2   = wf2.hlm[(2,2)]
    
    # extract waveform 1
    A_1, phi_1    = h22_1['A'], h22_1['p']
    u_1           = wf1.u
    imrg_1        = np.argmax(A_1)
    u_1_mrg       = u_1[imrg_1]
    Momg_1        = np.zeros_like(phi_1)
    Momg_1[1:]    = np.diff(phi_1)/np.diff(u_1)

    # extract waveform 2
    A_2, phi_2   = h22_2['A'], h22_2['p']
    u_2          = wf2.u
    imrg_2       = np.argmax(A_2)
    u_2_mrg      = u_2[imrg_2]
    Momg_2       = np.zeros_like(phi_2)
    Momg_2[1:]   = np.diff(phi_2)/np.diff(u_2)

    # shift mergers to same point
    u_1 = u_1 - u_1_mrg + u_2_mrg
    u_2 = u_2

    plt.plot(u_1, phi_1, label='A_1')
    plt.plot(u_2, phi_2, label='A_2')
    plt.show()

    if ref is not None:
        # align at reference time w.r.t merger
        t_mrg_1 = u_1[imrg_1]
        t_mrg_2 = u_2[imrg_2]
        i_ref_1 = utils.find_nearest(t_mrg_1+ref, u_1)
        i_ref_2 = utils.find_nearest(t_mrg_2+ref, u_2)
        dphi = phi_1[i_ref_1] - phi_2[i_ref_2]
        tau  = u_2[i_ref_2] - u_1[i_ref_1]
    else:
        # common time array
        u_new = np.linspace(max(u_1[0], u_2[0]), min(u_1[-1], u_2[-1]), 20000)

        win_start = u_2[utils.find_nearest(time_window[0], u_2)]
        win_end   = u_2[utils.find_nearest(time_window[1], u_2)]

        tau , dphi , _ = wf_utils.Align(u_new, win_end , win_end-win_start , u_1, phi_1, u_2, phi_2)
    
    return tau, dphi

if __name__ == '__main__':

    # parse command line arguments
    parser = argparse.ArgumentParser(description='Align and plot the SXS and ATK waveforms targeted to GW150914.')
    parser.add_argument('--t_window', type=float, nargs=2, default=[600, 1900], help='Time window for alignment')
    parser.add_argument('--path_atk', type=str, default='./', help='Path to the ATK waveform file')
    parser.add_argument('--path_sxs', type=str, default='./', help='Path to the SXS waveform file')
    args = parser.parse_args()

    time_window = args.t_window
    path_atk    = args.path_atk
    path_sxs    = args.path_sxs

    # SXS
    wf = sxs.Waveform_SXS(ID='0305', path=path_sxs, download=True, downloads=['hlm', 'metadata', 'horizons'], nu_rescale=False)

    # read ATK waveform
    with h5py.File(path_atk, 'r') as f:
        data = f['r0100']['h[2,2]'][()]
        t    = f['r0100']['t'][()]
    
    # create an empty waveform object
    wf_atk = Waveform()
    wf_atk._u = t
    wf_atk._hlm[(2,2)] = {'A': np.abs(data), 'p': -1*np.unwrap(np.angle(data))}

    colors     = ['#3A8DC5', '#C5723A']
    colors_dp  = ['#BAD827', '#27D8D8']

    t_mrg,_,_,_     = wf.find_max(kind='global')
    t_mrg_atk,_,_,_ = wf_atk.find_max(kind='global')

    # align the 22 mode
    tau, dphi = align_eobnr(wf, wf_atk, time_window=time_window)

    # shift NR
    wf._u     = wf.u - t_mrg + tau
    wf_atk._u = wf_atk.u - t_mrg_atk
    time_window -= t_mrg_atk
    mode = (2,2)

    fig, ((ax1, ax2), (pax1, pax2)) = plt.subplots(2, 2,figsize=(10,3.5), sharex='col', gridspec_kw={'width_ratios': [2, 1], 'height_ratios': [1.5, 1]})

    # plot the NR wfs
    for this_ax in [ax1, ax2]:
        this_ax.plot(wf.u, wf.hlm[mode]['z'],  color=colors[1], lw=1.5, label=r'$\tt{SXS:BBH:0305}$')

    this_dphi = dphi/2*mode[1]

    # remove phase diff
    common_time = np.linspace(max(wf.u[0], wf_atk.u[0]), min(wf.u[-1], wf_atk.u[-1]), 20000)
            
    # rel amp diff
    Anr  = np.interp(common_time, wf.u,   wf.hlm[mode]['A'])
    Aeob = np.interp(common_time, wf_atk.u, abs(wf_atk.hlm[mode]['A']))

    # phase
    wf_interp  = np.interp(common_time, wf.u,   wf.hlm[mode]['p'])
    eob_interp = np.interp(common_time, wf_atk.u, np.unwrap(wf_atk.hlm[mode]['p']))
    dphi_intp  = (wf_interp - eob_interp)
    dphi_intp -= this_dphi

    # plotting
    for this_ax in [ax1, ax2]:
        this_ax.plot(wf_atk.u, 
                        abs(wf_atk.hlm[mode]['A'])*np.cos(np.unwrap(wf_atk.hlm[mode]['p']) + this_dphi), 
                        color=colors[0], lw=1.5, ls='--', label='AthenaK'
                        )
        this_ax.axvline(0, color='k', lw=1.5, ls='--')
        this_ax.axvline(time_window[0], color='k', lw=1.5, ls=':')
        this_ax.axvline(time_window[1], color='k', lw=1.5, ls=':')
            
    for this_ax in [pax1, pax2]:
        this_ax.plot(common_time, dphi_intp,      color=colors_dp[0], linestyle='-', label=r'$\Delta\phi$')
        this_ax.plot(common_time, (Anr-Aeob)/Anr, color=colors_dp[1], linestyle='-', linewidth=1, label=r'$\Delta A/A$')
        this_ax.set_xlabel(r'$t-t_{\rm mrg} [M]$')
        this_ax.axvline(0, color='k', lw=1.5, ls='--')

    pax1.set_ylim(-0.5, 0.5)
    pax2.set_ylim(-3, 3)
    ax1.legend( frameon=False, ncol=2)
    pax1.legend( frameon=False, ncol=2)
    ax1.set_ylabel(rf'$\Re[h_{{{mode[0]}{mode[1]}}}/\nu$')
    ax1.set_xlim(-2000, 80)
    ax2.set_xlim(-200,  80)
    fig.tight_layout()
        
    plt.show()

