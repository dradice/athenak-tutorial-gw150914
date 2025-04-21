import argparse
from glob import glob
import h5py
import numpy as np
import os
import re

class Struct:
    pass

def bump_window(x, delta=0):
    """
    Double side, exponential windowing function on [0,1]
    """
    return np.piecewise(x,
            [(x == 0),
             (x > 0) & (x < delta),
             (x >= delta) & (x <= 1.0 - delta),
             (x > 1.0 - delta) & (x < 1),
             (x == 1)],
            [lambda x: 0.0,
             lambda x: np.exp(-1.0/(1.0 - ((x-delta)/delta)**2) + 1.0),
             lambda x: 1.0,
             lambda x: np.exp(-1.0/(1.0 - ((x-(1.0-delta))/delta)**2) + 1.0),
             lambda x: 0.0]
            )

def read_athenak_waveform(fname):
    """
    Reads an AthenaK waveform file and returns a tuple
        (t, psi4)
    where psi4 is a dictionary of l,m modes
    """
    with open(fname, "r") as f:
        psi4 = {}
        header = f.readline().split()
        for col in header[2:]:
            token = col.split(":")[1]
            l = int(token[0])
            m = int(token[1:])
            psi4[(l,m)] = None

    rawdata = np.loadtxt(fname, unpack=True, skiprows=1)
    # Resample data uniformly
    t = np.linspace(rawdata[0][0], rawdata[0][-1], rawdata[0].shape[0])
    for lm, d in zip(psi4.keys(), rawdata[1:]):
        psi4[lm] = np.interp(t, rawdata[0], d)

    return t, psi4

def fixed_freq_int_2(signal, fcut, dt=1):
    """
    Fixed frequency double time integration

    From Reisswig and Pollney, CQG 28 195015 (2011)

    * signal : a numpy array with the target signal
    * fcut   : the cutoff frequency
    * dt     : the sampling of the signal
    """
    from scipy.fftpack import fft, ifft, fftfreq
    from math import pi

    f = fftfreq(signal.shape[0], dt)

    idx_p = np.logical_and(f >= 0, f < fcut)
    idx_m = np.logical_and(f <  0, f > -fcut)

    f[idx_p] = fcut
    f[idx_m] = -fcut

    return ifft(-fft(signal)/(2*pi*f)**2)

class Waveforms:
    """
    AthenaK waveforms, all radii and l,m modes
    """
    def __init__(self, path_to_data, verbose=False):
        """
        Reads Psi4 from an AthenaK simulation

        * path_to_data : path to folder containing the data
        """
        self.path = path_to_data
        self.verbose = verbose

        # Check which radii are available
        radii = set([])
        for fnm in glob(os.path.join(self.path, r"rpsi4_real_*")):
            match = re.match(r"rpsi4_real_(\d+).txt", os.path.basename(fnm))
            if match is not None:
                rad = match.group(1)
                radii.add(rad)
        radii = sorted(list(radii), key=lambda x: float(x))

        # Read data for each radius
        self.rad = {}
        for r in radii:
            self.rad[r] = Struct()
            t, psi4 = self.__read_radius(r)
            self.rad[r].t = t
            self.rad[r].psi4 = psi4

    def compute_strain(self, fcut, twinfac):
        """
        Compute the strain for all l,m modes using the FFI method

        * fcut : cutoff frequency for the m=2 mode
                 (the frequencies for the other modes are obtained by rescaling)
        * win  : window size
        """
        if self.verbose:
            print("Computing strain... ", end="")
        for r in self.rad.keys():
            t = self.rad[r].t
            twin = twinfac * float(r)
            win = bump_window(t/t.max(), twin/t.max())
            self.rad[r].h = {}
            for lm in self.rad[r].psi4.keys():
                fc = 2*fcut/max(abs(lm[1]), 1)
                self.rad[r].h[lm] = fixed_freq_int_2(win*self.rad[r].psi4[lm], fc, dt=t[1]-t[0])
        if self.verbose:
            print("done!")

    def write_strain(self, fname):
        """
        Writes the strain data into an HDF5 file
        """
        if self.verbose:
            print(f"Writing {fname}... ", end="")
        with h5py.File(fname, "w") as f:
            for r in self.rad.keys():
                g = f.create_group(f"/r{r}")
                g.create_dataset("t", data=self.rad[r].t)
                for lm in self.rad[r].h.keys():
                    g.create_dataset(f"h[{lm[0]},{lm[1]}]", data=self.rad[r].h[lm])
        if self.verbose:
            print("done!")

    def __read_radius(self, rad):
        """
        Reads the data for a specific radius

        This function is not meant to be called by users

        * rad: radius to parse (string)

        Returns a dictionary psi4[(l,m)]
        """
        if self.verbose:
            print("Reading " + os.path.join(self.path, f"rpsi4_real_{rad}.txt") + "...", end='')
        t, psi4_re = read_athenak_waveform(os.path.join(self.path, f"rpsi4_real_{rad}.txt"))
        if self.verbose:
            print("done!")

        if self.verbose:
            print("Reading " + os.path.join(self.path, f"rpsi4_imag_{rad}.txt") + "...", end='')
        t, psi4_im = read_athenak_waveform(os.path.join(self.path, f"rpsi4_imag_{rad}.txt"))
        if self.verbose:
            print("done!")

        psi4 = {}
        for lm in psi4_re.keys():
            psi4[lm] = psi4_re[lm] + 1j*psi4_im[lm]
        return t, psi4

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("path", metavar="/path/to/waveforms", nargs=1,
                        help="Path to the waveform data to process")
    parser.add_argument("--fcut", default=0.0025, type=float,
                        help="Cutoff frequency for FFI")
    parser.add_argument("--plot", action=argparse.BooleanOptionalAction,
                        help="Plots a given component of the strain")
    parser.add_argument("-l", default=2, type=int, help="Which l-mode to plot")
    parser.add_argument("-m", default=2, type=int, help="Which m-mode to plot")
    parser.add_argument("--twin", default=0.5, type=float,
                        help="The window cutoff is twin * R_extract")
    parser.add_argument("--verbose", action=argparse.BooleanOptionalAction,
                        help="Print more output")
    args = parser.parse_args()
    dset = Waveforms(args.path[0], args.verbose)
    dset.compute_strain(args.fcut, args.twin)
    dset.write_strain(os.path.join(args.path[0], "strain.h5"))

    t = dset.rad["0200"].t
    h = dset.rad["0200"].h[args.l,args.m]
    A = np.abs(h)
    phi = np.unwrap(np.angle(h))
    omega = np.zeros_like(phi)
    omega[1:] = np.diff(phi)/np.diff(t)

    if args.plot:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(t, np.real(h), 'r', label=r"${\rm Re}\ h$")
        plt.plot(t, np.imag(h), 'b', label=r"${\rm Im}\ h$")
        plt.plot(t, A, 'k', label="$A$")
        plt.plot(t, -A, 'k')
        plt.plot(t, omega, color="g", label=r"$\omega$")
        plt.ylim(-0.5,0.5)
        plt.legend()
        plt.show()
