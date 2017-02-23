"""
This code is adapted for 80 columns
|------------------------------------------------------------------------------|
"""
import sys
try:
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    import ha5utils
except ImportError:
    sys.stderr.write("Could not execute program: " + str(sys.exc_info()[1])
                     + ".\nPlease install the module and try again")
    sys.stderr.flush()
    sys.exit()
pi = np.pi
c_0 = 2.99792458e8
mu_0 = pi*4e-7
eps_0 = 1/(mu_0*c_0**2)
eta_0 = 377


class Handin5(object):
    def __init__(self):
        """
        Class instantiation
        """
        self.txtTmpl = " "
        self.taskNbr = None
        self.figTxt = None  # the figure text object
        self.line = None  # plot
        self.line_2 = None  # plot 2
        self.lbda_0 = None
        self.lbdaZero = 550e-9  # center wavelength
        self.lbdaMax = 700e-9
        self.lbdaMin = 400e-9

    def simData(self):
        """
        This method updates the data that is to be plotted onto the figure.
        """

    def updatefig(self, simData):
        """
        Handles setting the new image/graph and info text for the simulation.
        Used by the animation function.

        Parameters
        ----------
        simData : tuple
        """
        if (self.taskNbr == 1):
            pass

    def initializePlot(self):
        """
        Initializes the 2D plot of the E-field amplitude vs propagated
        distance.
        """
        ymax = 1
        self.figTxt = ax1.text(self.lbdaMin*1e9, 1.05*ymax,
                               "", color="#000000")
        self.line, = ax1.plot([], [], linestyle="-", color="black")
        ax1.set_ylim(0, ymax)
        ax1.set_xlabel("$wavelength [nm]$", size=16)
        ax1.set_ylabel("$Reflectance$", size=16)
        ax1.set_xlim(self.lbdaMin*1e9, self.lbdaMax*1e9)

    def initializePlot_2(self):
        """
        Initializes the 2D plot of the E-field amplitude vs propagated
        distance.
        """
        ymax = 1
        self.line_2, = ax2.plot([], [], linestyle="-", color="black")
        ax2.set_ylim(0, ymax)
        ax2.set_xlabel("$wavelength [nm]$", size=16)
        ax2.set_ylabel("$Reflectance$", size=16)
        ax2.set_xlim(self.lbdaMin*1e9, self.lbdaMax*1e9)

    def PropTX(self, d, n, lambda_zero):
        """
        PropTX(self, d, n, lambda_zero)

        Transfer matrix for homogeneous material with thickness d
        and refractive index n, at (vacuum) wavelength lambda_zero.

        Parameters
        ----------
        d : float
            thickness of material
        n : numpy.ndarray
            refractive index of material (for each wavelength)
        lambda_zero : numpy.ndarray
            vacuum wavelength

        Returns
        -------
        numpy.ndarray
            The transfer matrix for the homogeneous media
        """
        k_zero = 2*pi/lambda_zero
        return np.array([[np.exp(1j*k_zero*n*d), np.zeros(len(lambda_zero))],
                         [np.zeros(len(lambda_zero)), np.exp(-1j*k_zero*n*d)]],
                        dtype=np.complex128)

    def BorderTX(self, nLeft, nRight):
        """
        BorderTX(self, nLeft, nRight)

        Transfer matrix for border between two media, with refractive index
        nLeft to the left and nRight to the right.

        Parameters
        ----------
        nLeft : numpy.ndarray
            refractive index of the left material.
        nRight : numpy.ndarray
            refractive index of the right material.

        Returns
        -------
        numpy.ndarray
            The transfer matrix for border between the two media.
        """
        return np.array([[nRight+nLeft, nRight-nLeft],
                         [nRight-nLeft, nRight+nLeft]],
                        dtype=np.complex128)/(2*nRight)

    def setupDBR(self, nStart, nHigh, nLow, nEnd, nbrPairs,
                 highestIdxFist=True, matrices=[]):
        """
        Sets up the matrices for propagating from the staring material to the
        end material through a dielectric mirror with given parameters.

        Parameters
        ----------
        nStart : float
            Refractive index of the material before the DBR.
        nHigh : float
            Refractive index of the material in the DBR with high refractive
            index
        nLow : float
            Refractive index of the material in the DBR with low refractive
            index
        nEnd : float
            Refractive index of the material after the DBR.
        nbrPairs : float
            Number of pairs of high and low refractive index in the DBR.
        highestIdxFist : bool
            True if the mirror starts with the high refractive index.
            False if the misrror starts with the low refractive index.
            Default value is True.
        matrices : list
            A list containing numpy.ndarray types that are the transfer
            matrices up to this point.
            Default value is [].

        Returns
        -------
        list
            A list containing numpy.ndarray types that are the transfer
            matrices.
        """
        nMirrFirst = nHigh
        nMirrLast = nLow
        if (not highestIdxFist):
            nMirrFirst = nLow
            nMirrLast = nHigh
        dMirrFirst = self.lbdaZero/(4*nMirrFirst)
        dMirrLast = self.lbdaZero/(4*nMirrLast)
        # propagate from border between start material and first layer to just
        # before border between first and second pair.
        matrices.extend([self.BorderTX(nStart, nMirrFirst),
                         self.PropTX(dMirrFirst, nMirrFirst, self.lbda_0),
                         self.BorderTX(nMirrFirst, nMirrLast),
                         self.PropTX(dMirrLast, nMirrLast, self.lbda_0)])
        # propagate from the border between first and second pair to border
        # between last pair and end material.
        for i in range(nbrPairs-1):
            matrices.extend([self.BorderTX(nMirrLast, nMirrFirst),
                             self.PropTX(dMirrFirst, nMirrFirst, self.lbda_0),
                             self.BorderTX(nMirrFirst, nMirrLast),
                             self.PropTX(dMirrLast, nMirrLast, self.lbda_0)])
        # propagate through the border between last pair and end material.
        matrices.extend([self.BorderTX(nMirrLast, nEnd)])
        return matrices

    def calcTX(self, matrices):
        """
        Calculates the reflected power (as a percentage of incident power)
        given the set of matrices

        Parameters
        ----------
        matrices : list
            A list of the transfer matrices starting with the final transfer
            matrix.

        Returns
        -------
        numpy.ndarray
            The reflected power as a percentage of them incident power.
        """
        nPoints = len(self.lbda_0)
        dims = np.array([2, 2, nPoints], dtype=np.int32)
        matTot = np.array([[np.ones(nPoints), np.zeros(nPoints)],
                           [np.zeros(nPoints), np.ones(nPoints)]],
                          dtype=np.complex128)
        for i in range(len(matrices), 0, -1):
            matTot = ha5utils.vecMatDot(dims, matrices[i-1], matTot)
        return np.abs(matTot[1, 0]/matTot[1, 1])**2

    def findMaxIdx(self, arr, startIdxHint=0):
        """
        TODO: Cythonize
        Finds the index in the given array where the maximum value of the array
        is found.

        Parameters
        ----------
        arr : array_like
            The array where the index of the maximum value is sought.
        startIdxHint : int
            A hint as to what index in them array the search should start from.
            Default value is 0.

        Returns
        -------
        list
            A list containing the indices where the maximum value of the array
            is found.
        """
        m = np.max(arr)
        for idx, val in enumerate(arr[startIdxHint:]):
            if (val == m):
                return startIdxHint + idx

    def findBandWidth(self, arr, peakValIdx, bwBounds):
        """
        TODO: Cythonize
        Attempts to find the bandwidth given the center point of the band and
        the lower bound of the band.

        Parameters
        ----------
        arr : array_like
            The array where the index of the maximum value is sought.
        peakValIdx : int
            The index peak value in the band (i.e. the center point of the
            band) in arr.
        bwBounds : float
            The lower bound of the band, i.e. any values in arr are considered
            to be outside the band.

        Returns
        -------
        """
        if (arr[peakValIdx] < bwBounds):
            # peak value is lower than them bounds
            return 0
        upperBound = peakValIdx
        lowerBound = peakValIdx
        for i in range(peakValIdx, len(arr)):
            if (arr[i] < bwBounds):
                upperBound = i
        for i in range(peakValIdx, 0, -1):
            if (arr[i] < bwBounds):
                lowerBound = i
        return upperBound - lowerBound

    def initTask1(self):
        """
        Initializes and executes task 1
        """
        self.taskNbr = 1
        self.lbdaZero = 550e-9  # center wavelength
        self.lbda_0 = np.arange(self.lbdaMin, self.lbdaMax, 1e-10,
                                dtype=np.float64)
        nStart = 1
        nSubst = 1.5  # substrate refr i
        nAR = np.sqrt(nSubst)  # anti-reflection coating refr i
        dAR = self.lbdaZero/(4*nAR)  # AR thickness
        # propagate from from start material through AR-coating to end material
        matrices = [self.BorderTX(nStart, nAR),
                    self.PropTX(dAR, nAR, self.lbda_0),
                    self.BorderTX(nAR, nSubst)]
        refl = self.calcTX(matrices)
        self.initializePlot()
        self.line.set_data(self.lbda_0*1e9, refl)
        maxidx = self.findMaxIdx(refl)[0]
        ax1.set_ylim(0, refl[maxidx])

    def initTask2(self):
        """
        Initializes and executes task 2
        """
        self.taskNbr = 2
        self.txtTmpl = "Maximum reflected power: %.2f percent at %.1f nm"
        self.lbdaZero = 633e-9
        self.lbda_0 = np.arange(self.lbdaMin, self.lbdaMax, 1e-10,
                                dtype=np.float64)
        nStart = 1
        nEnd = 1
        nHigh = 1.7
        nLow = 1.5
        matrices = self.setupDBR(nStart=nStart, nHigh=nHigh, nLow=nLow,
                                 nEnd=nEnd, nbrPairs=20)
        refl = self.calcTX(matrices)

        maxidx = self.findMaxIdx(refl)
        txts = self.txtTmpl % (refl[maxidx]*100, self.lbda_0[maxidx]*1e9)

        # plot data
        self.initializePlot()
        self.figTxt.set_text(txts)
        self.line.set_data(self.lbda_0*1e9, refl)

    def initTask3(self):
        """
        Initializes and executes task 3
        """
        self.taskNbr = 3
        self.txtTmpl = "Max refl pwr, w/o error: %.2f percent at %.1f nm; " +\
                       "Bandwidth @ %.1f nm: %.1f nm\n" +\
                       "Max refl pwr, w/ error: %.2f percent at %.1f nm; " +\
                       "Bandwidth @ %.1f nm: %.1f nm"
        self.lbdaZero = 980e-9
        self.lbdaMax = 1300e-9
        self.lbdaMin = 700e-9
        step = 1e-11
        self.lbda_0 = np.arange(self.lbdaMin, self.lbdaMax, step,
                                dtype=np.float64)
        nStart = 3.2
        nEnd = 1
        nHigh = 3.6
        nLow = 3.1
        
        matrices = self.setupDBR(nStart=nStart, nHigh=nHigh, nLow=nLow,
                                 nEnd=nEnd, nbrPairs=30)
        refl = self.calcTX(matrices)
        maxidx = self.findMaxIdx(refl, startIdxHint=5000)
        bw = self.findBandWidth(arr=refl, peakValIdx=maxidx,
                                bwBounds=0.99)*step

        matrices_2 = self.setupDBR(nStart=nStart, nHigh=nHigh, nLow=nLow,
                                   nEnd=nHigh, nbrPairs=15)
        # propagate from border between low index and high index with double
        # width to the border between low index with normal width and the next
        # pair.
        matrices.extend([self.PropTX(self.lbdaZero/(2*nHigh),
                                     nHigh, self.lbda_0),
                         self.BorderTX(nHigh, nLow),
                         self.PropTX(self.lbdaZero/(2*nLow),
                                     nLow, self.lbda_0)])
        matrices_2 = self.setupDBR(nStart=nLow, nHigh=nHigh, nLow=nLow,
                                   nEnd=nEnd, nbrPairs=14,
                                   matrices=matrices_2)
        refl_2 = self.calcTX(matrices_2)
        maxidx_2 = self.findMaxIdx(refl_2, startIdxHint=5000)
        bw_2 = self.findBandWidth(arr=refl_2, peakValIdx=maxidx_2,
                                  bwBounds=0.99)*step
        # plot data
        self.initializePlot()
        self.initializePlot_2()
        ax1.set_ylabel("$Reflectance$ $w/o$ $error$")
        ax2.set_ylabel("$Reflectance$ $w/$ $error$")
        txts = self.txtTmpl % (refl[maxidx]*100, self.lbda_0[maxidx]*1e9,
                               self.lbda_0[maxidx]*1e9, bw*1e9,
                               refl[maxidx_2]*100, self.lbda_0[maxidx_2]*1e9,
                               self.lbda_0[maxidx_2]*1e9, bw_2*1e9)
        self.figTxt.set_text(txts)
        self.line.set_data(self.lbda_0*1e9, refl)
        self.line_2.set_data(self.lbda_0*1e9, refl_2)


if __name__ == "__main__":
    """
    Called when the script runs.
    """
    fig = plt.figure()  # new figure
    hi5 = None
    if (len(sys.argv) > 1):
        args = sys.argv[1:]
        if (args[0] == "1"):
            ax1 = fig.add_subplot(111)
            hi5 = Handin5()  # the simulation class
            hi5.initTask1()  # method for task 1
        elif (args[0] == "2"):
            ax1 = fig.add_subplot(111)
            hi5 = Handin5()  # the simulation class
            hi5.initTask2()  # method for task 2
        elif (args[0] == "3"):
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)
            hi5 = Handin5()  # the simulation class
            hi5.initTask3()  # method for task 3
        else:
            sys.stdout.write("Usage: python <filename.py> <task_nbr>")
            sys.stdout.flush()
            sys.exit()
    else:
        sys.stdout.write("Usage: python <filename.py> <task_nbr>")
        sys.stdout.flush()
        sys.exit()
    """
    # the function handling the the animation itself
    ani = animation.FuncAnimation(fig, hi5.updatefig, hi5.simData,
                                  interval=20, blit=False, repeat=False)
    """
    plt.show()  # plot figure