""" 80 columns
--------------------------------------------------------------------------------
"""
import scipy as sp
import scipy.signal as sps
import scipy.io
import matplotlib.pyplot as plt


gen_data_dir = "HWE1/"


def load_gen_data():
    var_names = ["fc", "sr", "sr", "t", "tr", "ut"]
    out = {}
    for name in var_names:
        out.update({name: sp.array(
            scipy.io.loadmat(gen_data_dir + name + ".mat")[name][0])})
    return out


def time_spacing(t):
    return sp.average(t[1:] - t[:-1])


def sampling_frequency(t):
    dt = time_spacing(t)
    return 1/dt if dt != 0 else sp.inf


def nbr_samples(signal):
    return len(signal)


def frequency_vector(n_samples, f_s):
    """
    Calculates the zero-centered frequency vector

    Parameters
    ----------
    n_samples : int
        The number of sample points.
    f_s : float
        The sampling frequency.
    """
    return sp.arange(-n_samples/2, n_samples/2-1)*f_s/n_samples


def interpolate_signal(signal, time, interp_factor):
    """
    interpolate_signal(signal, time, interp_factor)

    This function interpolates the signal signal(time) by the interpolation
    factor interp_factor, using zero-padding in the frequency domain. The
    time vector time must be evenly spaced.
    The output si is the interpolated signal at the imterpolation times ti.

    Parameters
    ----------
    signal : numpy.ndarray
        The signal that is to be modulated
    time : numpy.ndarray
        time vector
    interp_factor : numpy.ndarray
        Interpolation factor

    Returns
    -------
    tuple
        numpy.ndarray
            The interpolated signal.
        numpy.ndarray
            The time vector interpolated signal.
    """
    n_samples = nbr_samples(signal)  # number of signal samples
    dt = time_spacing(time)  # signal time spacing
    signal_len = n_samples*interp_factor   # interpolated signal length
    # time vector for interpolate signal
    ti = time[0] + sp.arange(n_samples)*dt/interp_factor
    si = sps.resample(signal, signal_len)  # interpolated signal
    return si, ti


def matched_filter(unfiltered_signal, signal_filter):
    """
    matched_filter(unfiltered_signal,signal_filter)

    Calculates the filtered signal filtered_signal. unfiltered_signal is the
    signal to be filtered, and signal_filter is the (matched) filter.
    Hint:
    - The time reverse xrev=x(-t) of a vector x(t) can be obtained by the
      command xrev = x(end:-1:1)

    Parameters
    ----------
    unfiltered_signal : numpy.ndarray
        The signal that is to be filtered.
    signal_filter : numpy.ndarray
        The (matched) filter

    Returns
    -------
    numpy.ndarray
        The filtered signal
    """
    filtered_signal = sps.convolve(unfiltered_signal, signal_filter)
    # Remove the edge-parts of the filtered signal
    begin = int(len(signal_filter)/2)
    end = begin + len(unfiltered_signal)
    filtered_signal = filtered_signal[begin:end]
    return filtered_signal


def quadrature_demodulation(signal, time):
    """
    quadrature_demodulation(signal, time)

    This function performs quadrature demodulation on the signal signal(time).
    The output p is the complex phasor, p=I-iQ, where I and Q are the in-phase
    and quadrature components of the signal.
    Note: It is assumed that the apperent carrier frequency is f_s/4, thus the
    program lacks generality.

    Parameters
    ----------
    signal : numpy.ndarray
        The signar that is to be demodulated
    time : numpy.ndarray
        time vector

    Returns
    -------
    numpy.ndarray
        The complex phasor, p=I-iQ.
    """
    f_s = sampling_frequency(time)
    # Assume that the apperent carrier frequency is f_s/4
    fc_app = f_s/4
    f = frequency_vector(nbr_samples(signal), f_s)

    iq_td = signal*sp.exp(1j*2*sp.pi*time*fc_app)
    # FFT of the I-channel, before LP-filtering
    i_fd = sp.fftshift(sp.fft(sp.real(iq_td)))
    i_fd *= abs(f) > fc_app  # LP-filtering for I
    # FFT of the Q-channel, before LP-filtering
    q_fd = sp.fftshift(sp.fft(sp.imag(iq_td)))
    q_fd *= abs(f) > fc_app  # LP-filtering for Q
    # Combine the I and Q channels in the frequency domain
    iq_fd = i_fd - 1j*q_fd
    p = 2*sp.ifft(sp.ifftshift(iq_fd))  # Complex phasor P
    return p


def signal_plot(signal, time):
    """
    [signal_fd,f]=signal_plot(signal,time)

    Plots the real part of the signal as a function of time, as well as the
    magnitude of the Fourier transform of signal (denoted signal_fd) as a
    function of frequency f. signal_fd and f are also given as output
    variables.
    The time vector (in seconds) is assumed to be evenly spaced.
    """
    f_s = sampling_frequency(time)
    n_samples = nbr_samples(signal)
    signal_fd = sp.fftshift(sp.fft(signal))/n_samples
    f = frequency_vector(n_samples, f_s)  # zero-centered

    fig = plt.figure()

    # Plot the real part of the signal in the time domain
    ax1 = fig.add_subplot([2, 1, 1])
    ax1.plot(time*1e6, sp.real(signal))
    ax1.set_xlabel('Time [\mus]')
    ax1.set_ylabel('Real part of the signal')

    # Plot the magnitude of the Fourier transform
    ax2 = fig.add_subplot([2, 1, 2])
    ax2.plot(f*1e-6, abs(signal_fd)/max(abs(signal_fd)))
    ax2.set_xlabel('Frequency [MHz]')
    ax2.set_ylabel('Magnitude of the Fourier transform.')
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    pass
