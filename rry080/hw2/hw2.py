import matplotlib.pyplot as plt
import numpy as np
import numpy.fft as npfft
import scipy as sp
import scipy.signal as sps
import scipy.io as spi


def createHfilter():
    R0_near = 10446  # Distance to near range [m]
    delta_R = 4  # Pixel size in range direction [m]
    nRows = 350  # Number of rows
    lbda = 0.0566  # Wavelength [m]
    n = np.arange(-700, 701)  # Update n after you determine the filter length!
    PRF_v = 2.32  # (Put your PRF_v here.)  # PRF divided by aircraft speed [m^-1]
    Hfilter = np.zeros((nRows, len(n)))
    for row in range(nRows):
        R = 0  # (R should be dependent on row number.)
        Hfilter[row] = 0  # (Put your filter expression here.)
    return Hfilter


def fftrows(signal):
    return npfft.fft(signal, 4096, 1)


def ifftrows(signal):
    return npfft.ifft(signal, 4096, 1)


def filterplot(H):
    H = H[1]
    t1 = np.shape(H)
    f = np.linspace(-.5, .5, t1[2])

    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax1.plot(f, 20*np.log10(abs(npfft.fftshift(npfft.fft(H)))))
    ax1.set_xlabel('Uniform frequnency')
    ax1.set_ylabel('Power (dB)')
    ax1.set_title('Filter properties')
    ax1.grid(True)

    ax2 = fig.add_subplot(312)
    a = np.unwrap(np.angle(npfft.fftshift(npfft.fft(H))))
    ax2.plot(f, a)
    ax2.set_xlabel('Uniform frequency')
    ax2.set_ylabel('Angle (rad)')
    ax2.axis([-.5, .5, min(a), 0])
    ax2.grid(True)

    ax3 = fig.add_subplot(313)
    a = sps.group_delay((H, 1))
    ax3.plot(np.linspace(-.5, .5, 512), npfft.fftshift(a) - t1(2)/2)
    ax3.set_xlabel('Uniform frequency')
    ax3.set_ylabel('Group delay')
    ax3.axis([-.5, .5, -500, 500])
    ax3.grid(True)


def imageplot(signal):
    plt.imshow(np.sqrt(np.sqrt(abs(signal))), cmap='gray')
    # brighten(.2)
    plt.show()


def sarfilter(H, signal):
    """
    image = sarfilter(H, signal);

    Apply the SAR filter H to the radar data signal
    using Fourier transforms.
    H should be a matched filter.

    Patrik Dammert 1998-02-23

    Smal modifications made to avoid alliasing effects.
    Björn Hallberg 2002-12-27
    get right linear convolution from ifft.
    Wiebke Aldenhoff 2016-05-24
    """
    rows, cols = np.shape(H)
    if (rows != 350) or (cols > 2000):
        print('Wrong size of filter input!')
        return

    image = np.zeros((350, 4096))

    # Choose between
    # 1 - Fourier transform of whole matrix (uses MUCH memory)
    # 2 - Fourier transform of row-by-row (uses little memory)
    choice = 2

    if choice == 1:
        # This implementation will be affected by alliasing
        image = ifftrows(fftrows(H)*fftrows(signal))
    elif choice == 2:
        for s in range(350):
            temp = npfft.ifft(npfft.fft(H[s], 8192)*npfft.fft(signal[s], 8192))
            # center portion of the convolution
            a = np.floor(np.size(H, axis=1)/2) + 1
            b = np.floor(sp.size(H, axis=2)/2) + 4096
            image[s] = temp[a:b]
    return image


def spectrumplot(signal, nlines):
    """
    spectrumplot(signal,nlines)

    M-file for easy access to estimates of
    a signal's spectrum.

    Signal is the raw data and nlines is the number
    of lines you want in the estimate of the
    spectrum.

    Version 1.0 Patrik Dammert 950523
    Version 1.1 Anders Berg    100827 clean-up
    """
    powof2 = 256
    size2 = int(np.size(signal[:nlines]) / powof2)
    if powof2*size2 != np.size(signal):
        print("resize of shape (%d, %d) not available for data of length %d"
              % (powof2, size2, np.size(signal)))
    signal_reshaped = np.transpose(np.reshape(
        np.transpose(signal[:nlines]),
        (powof2, size2)))
    print(np.shape(signal_reshaped), np.shape(signal))

    # figure;
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(np.linspace(-155, 155, powof2),
            np.mean(abs(npfft.fftshift(
                npfft.fft(signal_reshaped, n=size2, axis=1)))**2))
    ax.grid(True)
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Power')
    ax.seT_title('Power Spectral Density of azimuth spectrum')


def testfun(signal):
    return np.average(npfft.fft(signal, n=4096, axis=0), axis=0)


if __name__ == "__main__":
    data = spi.loadmat("sarlab.mat")
    signal = data["signal"]
    oneline = data["oneline"]
    spectrumplot(signal, len(signal))
