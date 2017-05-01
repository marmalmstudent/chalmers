""" 80 columns
--------------------------------------------------------------------------------
"""
import numpy as np
import numpy.fft as npfft
import scipy as sp
import scipy.optimize as spo
import scipy.signal as sps
import scipy.io
import matplotlib.pyplot as plt
import ha1utils


gen_data_dir = "HWE1/"


def load_gen_data():
    """
    Loads the generated data contained in the .mat files.

    Returns
    -------
    tuple
        float
            the carrier frequency.
        numpy.ndarray
            The received signal.
        numpy.ndarray
            The transmitted signal.
        numpy.ndarray
            The time vector specifying as which time points the transmitted
            signal is sampled at.
        numpy.ndarray
            The time vector specifying as which time points the received
            signal is sampled at.
        numpy.ndarray
            The complex envelope of the transmitted signal.
    """
    var_names = ["fc", "sr", "st", "t", "tr", "ut"]
    out = {}
    for name in var_names:
        out.update({name: np.array(
            scipy.io.loadmat(gen_data_dir + name + ".mat")[name][0])})
    return out


def fft(signal):
    """
    Short for npfft.fftshift(npfft.fft(signal))
    The transformed signal needs to be normalized.
    """
    return npfft.fftshift(npfft.fft(signal))


def time_spacing(t):
    """
    Calculates the average time spacing of the signal.

    Parameters
    ----------
    t : numpy.ndarray
        The sampling time vector.

    Returns
    -------
    float
        The average time_spacing of the signal.
    """
    return sp.average(t[1:] - t[:-1]) if isinstance(t, np.ndarray) else 0


def sampling_frequency(t):
    """
    Calculates the sampling frequency based on the sampling time.

    Parameters
    ----------
    t : numpy.ndarray
        The sampling time vector.

    Returns
    -------
    float
        The sampling frequency
    """
    dt = time_spacing(t)
    return 1/dt if dt != 0 else sp.inf


def nbr_samples(signal):
    """
    Calculates the number of sampling points the signal has.

    Parameters
    ----------
    signal : numpy.ndarray
        The sampled signal.

    Returns
    -------
    int
        The number of sampling points.
    """
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
    return sp.arange(-n_samples/2, n_samples/2)*f_s/n_samples


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
    i_fd = npfft.fftshift(npfft.fft(sp.real(iq_td)))
    i_fd *= abs(f) <= fc_app  # LP-filtering for I
    # FFT of the Q-channel, before LP-filtering
    q_fd = npfft.fftshift(npfft.fft(sp.imag(iq_td)))
    q_fd *= abs(f) <= fc_app  # LP-filtering for Q
    # Combine the I and Q channels in the frequency domain
    iq_fd = i_fd - 1j*q_fd
    p = 2*sp.ifft(np.fft.ifftshift(iq_fd))  # Complex phasor P
    return p


def signal_plot(signal, time):
    """
    signal_plot(signal,time)

    Plots the real part of the signal as a function of time, as well as the
    magnitude of the Fourier transform of signal (denoted signal_fd) as a
    function of frequency f. signal_fd and f are also given as output
    variables.
    The time vector (in seconds) is assumed to be evenly spaced.
    """
    f_s = sampling_frequency(time)
    n_samples = nbr_samples(signal)
    signal_fd = npfft.fftshift(npfft.fft(signal))/n_samples
    f = frequency_vector(n_samples, f_s)  # zero-centered

    fig = plt.figure()

    # Plot the real part of the signal in the time domain
    ax1 = fig.add_subplot(211)
    ax1.plot(time*1e6, sp.real(signal))
    ax1.set_xlabel('Time [$\mu$s]')
    ax1.set_ylabel('Real part of the signal')

    # Plot the magnitude of the Fourier transform
    ax2 = fig.add_subplot(212)
    ax2.plot(f*1e-6, abs(signal_fd)/max(abs(signal_fd)))
    ax2.set_xlabel('Frequency [MHz]')
    ax2.set_ylabel('Magnitude of the Fourier transform.')
    plt.grid(True)
    # plt.show()


def find_bandiwdth(signal, center_idx, threshold=1/np.sqrt(2)):
    """
    Calculates the bandwidth of the signal based on the given threshold.

    Parameters
    ----------
    signal : numpt.ndarray
        The signal.
    center_idx : int
        The index in signal array where the frequency is found.
    threshold float
        The fractional amplitude threshold which determines the bandwidth.
    """
    bw_cond = signal[center_idx]*threshold
    upper = lower = center_idx
    for i in range(center_idx, len(signal)):
        if signal[i] < bw_cond:
            upper = i
            break
    for i in range(center_idx, -1, -1):
        if signal[i] < bw_cond:
            lower = i
            break
    return lower, upper


def task1(tx_envelope, tx_sample_time, tx_signal, rx_sample_time, rx_signal):
    # plot transmitted envelope
    # signal_plot(tx_envelope, tx_sample_time)

    # find bandwidth to estimate amplitude of signal
    center_idx = ha1utils.find_ordered_zero(
        tx_sample_time, int(len(tx_sample_time)-1), int(0))
    # assum center_idx != -1

    """ pulse width
    u(t) = 1/sqrt(t_p)*rect(t/t_p)*exp(-i*pi*gma*t^2)
    choose t = 0 (-> rect(t/t_p) = 1) and solve |u(t)| = |1/sqrt(t_p)|
    """
    pulse_width = 1/np.abs(tx_envelope[center_idx])**2
    print("Pulse width: %.1f us" % (pulse_width*1e6))

    # bandwidth
    n_sampl = nbr_samples(tx_envelope)
    f_sampl = sampling_frequency(tx_sample_time)
    freq = frequency_vector(n_sampl, f_sampl)
    env_fft = abs(fft(tx_envelope))
    bandwidth = np.transpose(ha1utils.findbw(
        env_fft/max(env_fft), 1/np.sqrt(2), 5e-2))
    lower = bandwidth[0]
    upper = bandwidth[1]
    bw = (upper-lower) * f_sampl/n_sampl
    print("Bandwidth: %.1f MHz" % (bw[0]*1e-6))

    # bandwidth product
    print("Bandwidth product: %.1f" % (pulse_width*bw[0]))

    # nyquist sampling rate (should be 2*f_max)
    nyquist_r = 2*freq[upper[0]]
    f_s = sampling_frequency(tx_sample_time)
    if f_s < nyquist_r:
        sampling_text = "\033[91m>%s<\033[0m" % "undersampled"
    else:
        sampling_text = "\033[92m>%s<\033[0m" % "oversampled"
    print("Nyquist rate: %.1f MHz; Sampling rate: %.1f MHz."
          % (nyquist_r*1e-6, f_s*1e-6))
    print("\t-> The signal is %s" % sampling_text)

    sig_fft = abs(fft(tx_signal))
    bandwidth = np.transpose(ha1utils.findbw(
        sig_fft/max(sig_fft), 1/np.sqrt(2), 5e-2))
    lower = bandwidth[0]
    upper = bandwidth[1]
    center = np.array((upper+lower)/2, dtype=np.uint32)
    f_center = freq[center[1]+1]  # positive part, assume only 2 bands.
    print("Apparent carrier frequency: %.1f MHz" % (f_center*1e-6))
    print("Carrier frequency: %.1f MHz" % (carrier_freq*1e-6))
    # apparent carrier frequency and carrier frequency different because
    # f_ac = f_c - floor(f_c/f_s) which is a result from an fftshift

    # plot transmitted signal
    # signal_plot(tx_signal, tx_sample_time)

    # plot received signal
    # signal_plot(rx_signal, rx_sample_time)
    # can not detect any target because of te noise


def task2(rx_sample_time, rx_signal, tx_signal):
    rx_quad_demod = quadrature_demodulation(rx_signal, rx_sample_time)
    """
    signal_plot(tx_envelope, tx_sample_time)
    signal_plot(rx_quad_demod, rx_sample_time)
    plt.show()
    """
    print(len(tx_sample_time), len(rx_sample_time), len(rx_sample_time)/len(tx_sample_time))
    t_s = time_spacing(tx_sample_time)
    t = np.arange(tx_sample_time[0], rx_sample_time[-1], t_s)
    s_tx = np.zeros(len(t), dtype=np.complex128)
    s_rx = np.zeros(len(t), dtype=np.complex128)
    """
    s_tx = np.hstack((np.flipud(tx_signal.conjugate()), np.zeros(len(t)-len(tx_signal),
                                          dtype=np.complex128)))
    s_rx = np.hstack((rx_quad_demod, np.zeros(len(t)-len(rx_quad_demod),
                                              dtype=np.complex128)))
    """
    middle_1 = int((np.average(rx_sample_time)-tx_sample_time[0])/t_s)
    s_tx[middle_1-int(len(tx_sample_time)/2):middle_1+int(len(tx_sample_time)/2)+1] = np.flipud(tx_signal.conjugate())
    s_rx[middle_1-int(len(rx_sample_time)/2)-1:] = rx_quad_demod

    rx_filteded = matched_filter(s_rx, s_tx)
    # plt.plot(t, np.real(rx_filteded))
    # plt.show()


if __name__ == "__main__":
    gen_data = load_gen_data()
    carrier_freq = gen_data["fc"]
    rx_signal = gen_data["sr"]
    tx_signal = gen_data["st"]
    tx_sample_time = gen_data["t"]
    rx_sample_time = gen_data["tr"]
    tx_envelope = gen_data["ut"]
    print(len(rx_sample_time), rx_sample_time[-1])
    print(len(tx_sample_time), tx_sample_time[0])
    
    # task1(tx_envelope, tx_sample_time, tx_signal, rx_sample_time, rx_signal)
    task2(rx_sample_time, rx_signal, tx_signal)
