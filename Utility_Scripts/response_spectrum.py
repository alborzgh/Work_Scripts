import numpy as np
import matplotlib.pyplot as plt

def response_spectrum(a, dt, damp=0.05, periods=None, nperiod=100, kind='log'):
    """
    This function builds response spectra from acceleration time history. The 
    results are the pseudo-spectral acceleration and velocities.
    
    damp    = the damping ratio for which the response spectrum is calculated for
    nperiod = number of periods at which response spectrum is calculated
    kind    = linear/log/mix - defines the distribution of the points
    """

    if periods is None:
        # define range of considered periods by power of 10
        minpower = -3.0
        maxpower = 1.0

        # create vector of considered periods
        if kind == "lin" or kind == "linear":
            p = np.linspace(10**minpower, 10**maxpower, nperiod)
        elif kind == "log" or kind == "logarithmic":
            p = np.logspace(minpower, maxpower, nperiod)
        elif kind == "mix" or kind == "mixed":
            p = np.logspace(minpower, 0, nperiod)
            p = np.append(p, np.linspace(1, 10**maxpower, nperiod))
            nperiod *= 2
        else:
            raise ValueError('kind can only be lin, log or mix.')
    else:
        p = periods
        nperiod = len(periods)
        
    # vector of circular freq
    w = 2.0 * np.pi / p

    # fast fourier transform of acceleration
    afft = np.fft.fft(a, len(a))
    freq = np.fft.fftfreq(len(a), d = dt)
    wbar = 2.0 * np.pi * freq

    # calculate the transfer function
    temp1 = np.tile(w*w, [len(a), 1])
    temp2 = np.tile(wbar*wbar, [nperiod, 1]).transpose()
    temp3 = 2.0 * 1j * damp * np.outer(wbar, w)
    Hw = 1.0 / (temp2 - temp1 - temp3)

    # calculate the response of the system (displacement)
    u =  Hw * np.tile(afft, [nperiod, 1]).transpose()

    # transform the results to the time domain
    utime = np.fft.ifft(u, axis=0)

    # calculate other spectral parameters
    S_d = np.amax(np.abs(np.real(utime)), axis=0)
    S_v = w * S_d
    S_a = w * S_v

    return p, S_d, S_v, S_a

def main():
    time = np.linspace(0,10,10000)
    acc  = np.sin(2.0 * np.pi * 2 * time)

    p, dmax, vmax, amax = response_spectrum(acc, time[1] - time[0])

    plt.figure()
    plt.subplot(3,1,1)
    plt.plot(p, amax)
    plt.subplot(3,1,2)
    plt.plot(p, vmax)
    plt.subplot(3,1,3)
    plt.plot(p, dmax)
    plt.show()

if __name__ == "__main__":
    main()
