""" calcs.py
This code is adapted for 80 columns
|------------------------------------------------------------------------------|
"""
import sys
try:
    import numpy as np
    import matplotlib.pyplot as plt
    import ha5utils as ha5u
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
        self.txtTmpl = " "
        self.figTxt = None  # the figure text object
        self.line_0 = None
        self.line_1 = None
        self.lbda_0 = None
        self.lbdaZero = 980e-9  # center wavelength
        self.lbdaMax = 1300e-9
        self.lbdaMin = 700e-9
        self.step = 1e-10  # lambda step separaton

    def initializePlot(self):
        """
        Initializes the 2D plot of the reflected power vs wavelength.
        """
        ymax = 1
        self.figTxt = ax0.text(self.lbdaMin*1e9, 1.05*ymax,
                               "", color="#000000")
        self.line_0, = ax0.plot([], [], linestyle="-", color="black")
        ax0.set_ylim(0, ymax)
        ax0.set_xlabel("$wavelength [nm]$", size=16)
        ax0.set_ylabel("$Reflectance$", size=16)
        ax0.set_xlim(self.lbdaMin*1e9, self.lbdaMax*1e9)

    def initializePlot_2(self):
        """
        Initializes the 2D plot of the E-field amplitude vs propagated
        distance.
        """
        ymax = 1
        self.line_1, = ax1.plot([], [], linestyle="-", color="black")
        ax1.set_ylim(0, ymax*1.01)
        ax1.set_xlabel("$wavelength [nm]$", size=16)
        ax1.set_ylabel("$Reflectance$", size=16)
        ax1.set_xlim(self.lbdaMin*1e9, self.lbdaMax*1e9)

    def PropTXM(self, d, n):
        """
        PropTXM(self, d, n)

        Transfer matrix for homogeneous material with thickness d
        and refractive index n.

        Parameters
        ----------
        d : float
            thickness of material
        n : numpy.ndarray
            refractive index of material (for each wavelength)

        Returns
        -------
        numpy.ndarray
            The transfer matrix for the homogeneous media
        """
        k_zero = 2*pi/self.lbda_0
        zrs = 0
        if (np.size(k_zero) > np.size(d)):
            zrs = np.zeros(np.size(k_zero))
        elif (np.size(k_zero) < np.size(d)):
            zrs = np.zeros(np.size(d))
        # if np.size(k_zero) == np.size(d) they are likely both equal to 1.
        return np.array([[np.exp(1j*k_zero*n*d), zrs],
                         [zrs, np.exp(-1j*k_zero*n*d)]],
                        dtype=np.complex128)

    def BorderTXM(self, nLeft, nRight):
        """
        BorderTXM(self, nLeft, nRight)

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
                 matrices, highestIdxFist=True):
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
        matrices : list
            A list containing numpy.ndarray types that are the transfer
            matrices up to this point.
        highestIdxFist : bool
            True if the mirror starts with the high refractive index.
            False if the misrror starts with the low refractive index.
            Default value is True.

        Returns
        -------
        list
            A list containing numpy.ndarray types that are the transfer
            matrices.
        """
        if (nbrPairs < 1):
            return matrices

        nFirstMirr = nHigh
        nLastMirr = nLow
        if (not highestIdxFist):
            nFirstMirr = nLow
            nLastMirr = nHigh

        dFirstMirr = self.lbdaZero/(4*nFirstMirr)
        dLastMirr = self.lbdaZero/(4*nLastMirr)

        # propagate from border between start material and first layer to just
        # before border between first and second pair.
        matrices.extend([self.BorderTXM(nStart, nFirstMirr),
                         self.PropTXM(dFirstMirr, nFirstMirr),
                         self.BorderTXM(nFirstMirr, nLastMirr),
                         self.PropTXM(dLastMirr, nLastMirr)])
        # propagate from the border between first and second pair to border
        # between last pair and end material.
        for i in range(1, nbrPairs):
            matrices.extend([self.BorderTXM(nLastMirr, nFirstMirr),
                             self.PropTXM(dFirstMirr, nFirstMirr),
                             self.BorderTXM(nFirstMirr, nLastMirr),
                             self.PropTXM(dLastMirr, nLastMirr)])
        # propagate through the border between last pair and end material.
        matrices.extend([self.BorderTXM(nLastMirr, nEnd)])
        return matrices

    def calcTXM(self, matrices, nPoints):
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
            The reflected power as a percentage of the incident power.
        """
        dims = np.array([2, 2, nPoints], dtype=np.int32)
        matTot = np.eye(2, dtype=np.complex128)
        for i in range(len(matrices)):
            matTot = ha5u.vecMatDot(dims, matrices[i], matTot)
        return np.abs(matTot[1, 0]/matTot[1, 1])**2

    def initTask1(self):
        """
        Displays then reflected power relative to then incident power for a
        wavelength spectrum. The setup uses an anti-reflective coating designed
        for 550e nm.
        """
        self.lbdaZero = 550e-9
        self.lbdaMax = 700e-9
        self.lbdaMin = 400e-9
        self.lbda_0 = np.arange(self.lbdaMin, self.lbdaMax, self.step,
                                dtype=np.float64)
        nStart = 1
        nSubst = 1.5  # substrate refr i
        nAntiRelf = np.sqrt(nSubst)  # anti-reflection coating refr i
        dAntiRelf = self.lbdaZero/(4*nAntiRelf)  # AR thickness

        # propagate from from start material through AR-coating to end material
        matrices_1 = [self.BorderTXM(nStart, nAntiRelf),
                      self.PropTXM(dAntiRelf, nAntiRelf),
                      self.BorderTXM(nAntiRelf, nSubst)]
        refl_1 = self.calcTXM(matrices_1, nPoints=len(self.lbda_0))

        # setup plot
        self.initializePlot()
        self.line_0.set_data(self.lbda_0*1e9, refl_1)
        # maxidx_1 = ha5u.findMaxIdx(refl_1)
        # ax0.set_ylim(0, refl_1[maxidx_1])

    def initTask2(self):
        """
        Displays the reflected power relative to the incident power for a
        wavelength spectrum. Ths setup uses a dielectric mirror designed for
        633 nm consisting of 20 pairs of high and low refractive indices.
        """
        self.lbdaZero = 633e-9
        self.lbdaMax = 700e-9
        self.lbdaMin = 400e-9
        self.lbda_0 = np.arange(self.lbdaMin, self.lbdaMax, self.step,
                                dtype=np.float64)
        nStart = 1
        nEnd = 1
        nHigh = 1.7
        nLow = 1.5

        # propagate
        matrices_1 = self.setupDBR(nStart=nStart, nHigh=nHigh, nLow=nLow,
                                   nEnd=nEnd, nbrPairs=20,
                                   matrices=[])
        refl_1 = self.calcTXM(matrices_1, nPoints=len(self.lbda_0))

        # setup plot
        self.txtTmpl = "Maximum reflected power: %.2f percent at %.1f nm"
        self.initializePlot()
        maxidx_1 = ha5u.findMaxIdx(refl_1)
        txts = self.txtTmpl % (refl_1[maxidx_1]*100, self.lbda_0[maxidx_1]*1e9)
        self.figTxt.set_text(txts)
        self.line_0.set_data(self.lbda_0*1e9, refl_1)

    def initTask3(self):
        """
        Displays the reflected power relative to the incident power for a
        wavelength spectrum. The plot consists of two figures. In then first
        one a dielectric mirror designed for 980 nm consisting of 30 pairs of
        high and low refractive indices is used and in then second one the same
        mirror is used but it contains a segment with double thickness.
        """
        self.lbda_0 = np.arange(self.lbdaMin, self.lbdaMax, self.step,
                                dtype=np.float64)
        nStart = 3.2
        nEnd = 1
        nHigh = 3.6
        nLow = 3.1

        # propagate w/o error
        matrices_1 = self.setupDBR(nStart=nStart, nHigh=nHigh, nLow=nLow,
                                   nEnd=nEnd, nbrPairs=30,
                                   matrices=[])
        refl_1 = self.calcTXM(matrices_1, nPoints=len(self.lbda_0))
        maxidx_1 = ha5u.findMaxIdx(refl_1,
                                   int((self.lbdaZero - self.lbdaMin) /
                                       self.step*0.9))

        # propagate w/ error
        matrices_2 = self.setupDBR(nStart=nStart, nHigh=nHigh, nLow=nLow,
                                   nEnd=nHigh, nbrPairs=15,
                                   matrices=[])
        # propagate from border between low index and high index with double
        # width to the border between low index with normal width and the next
        # pair.
        matrices_2.extend([self.PropTXM(self.lbdaZero/(2*nHigh), nHigh),
                           self.BorderTXM(nHigh, nLow),
                           self.PropTXM(self.lbdaZero/(4*nLow), nLow)])
        matrices_2 = self.setupDBR(nStart=nLow, nHigh=nHigh, nLow=nLow,
                                   nEnd=nEnd, nbrPairs=14,
                                   matrices=matrices_2)
        refl_2 = self.calcTXM(matrices_2, nPoints=len(self.lbda_0))
        maxidx_2 = ha5u.findMaxIdx(refl_2,
                                   int((self.lbdaZero - self.lbdaMin) /
                                       self.step*0.9))

        # setup plot
        self.txtTmpl = "Max refl pwr, w/o error: %.2f percent at %.1f nm; " +\
                       "Bandwidth @ %.1f nm: %.1f nm\n" +\
                       "Max refl pwr, w/ error: %.2f percent at %.1f nm; " +\
                       "Bandwidth @ %.1f nm: %.1f nm"
        self.initializePlot()
        self.initializePlot_2()
        ax0.set_ylabel("$Reflectance$ $w/o$ $error$")
        ax1.set_ylabel("$Reflectance$ $w/$ $error$")
        bw_1 = ha5u.findBandWidth(refl_1, maxidx_1, 0.99)*self.step
        bw_2 = ha5u.findBandWidth(refl_2, maxidx_2, 0.99)*self.step
        txts = self.txtTmpl % (refl_1[maxidx_1]*100, self.lbda_0[maxidx_1]*1e9,
                               self.lbda_0[maxidx_1]*1e9, bw_1*1e9,
                               refl_1[maxidx_2]*100, self.lbda_0[maxidx_2]*1e9,
                               self.lbda_0[maxidx_2]*1e9, bw_2*1e9)
        self.figTxt.set_text(txts)
        self.line_0.set_data(self.lbda_0*1e9, refl_1)
        self.line_1.set_data(self.lbda_0*1e9, refl_2)
        ax0.set_ylim(0, refl_1[maxidx_1]*1.01)
        ax1.set_ylim(0, refl_2[maxidx_2]*1.01)

    def initTask4(self):
        """
        shows how the sing of the imaginary part of the refractive index
        affects wether the material absorbs or amplifies the signal.
        """
        print("Positive imaginary parts -> material is lossy.\n\t" +
              "e^{1j*k_0(n\'+1j*n\")z} = e^{1j*k_0*n\'z}*e^{-k_0*n\"z)}")

    def initTask5(self):
        """
        Calculates the complex refractive index based on alpha=1000 cm^-1.
        """
        lbda = 980e-9
        alpha = 1e5
        print("\tE = \te^{1j*k_0*(n\'+1j*n\")*z} = " +
              "e^{1j*k_0*n\'*z}*e^{-k_0*n\"*z)}")
        print("\t <=>\t|E| = e^{-k_0*n\"*z) = e^{-alpha*z}}")
        print("\t <=>\talpha = 2*k_0*n\" = 4*pi*n\"/lambda")
        print("\t <=>\tn\" = alpha*lambda/(4*pi)\n")
        print("\talpha\t= 1e3 cm^-1 = 1e5 m^-1")
        print("\tlambda\t= 980e-9 nm")
        print("\t  =>\tn\" = %.3f" % (alpha*lbda/(4*pi)*1e3) + "e-3")

    def initTask6(self):
        """
        Displays reflected power as a fraction of incident power for a wave
        propagating through a 20 nm thick lossy plate, then some distance
        varying from 0 um to 1000 nm and then through a dielectric mirror
        designed for 980 nm with 30 layers.
        """
        self.lbda_0 = 980e-9
        self.step = 1e-9
        nStart = 3.2
        nEnd = 1
        nHigh = 3.6
        nLow = 3.1
        alpha = 1e5
        nLossy = nStart + 1j*alpha*self.lbdaZero/(4*pi)
        ds = np.arange(0, 1000e-9, self.step)

        # propagate from border between low index and high index with double
        # width to the border between low index with normal width and the next
        # pair.
        matrices = ([self.BorderTXM(nStart, nLossy),
                     self.PropTXM(20e-9, nLossy),
                     self.BorderTXM(nLossy, nStart),
                     self.PropTXM(ds, nStart)])
        matrices = self.setupDBR(nStart=nStart, nHigh=nHigh, nLow=nLow,
                                 nEnd=nEnd, nbrPairs=30,
                                 matrices=matrices)
        refl = self.calcTXM(matrices, nPoints=len(ds))

        # stup plot
        self.txtTmpl = "Max refl pwr, w/o error: %.2f percent at $s = %.1f$ nm"
        self.initializePlot()
        self.figTxt = ax0.text(0, 1.05,
                               "", color="#000000")
        ax0.set_ylabel("$Reflectance$ $@$ $980$ $nm$")
        print(len(refl))
        maxidx_1 = ha5u.findMaxIdx(refl, 0)
        bw_1 = ha5u.findBandWidth(refl, maxidx_1, 0.99)*(ds[1]-ds[0])
        txts = self.txtTmpl % (refl[maxidx_1]*100, bw_1*1e9)
        self.figTxt.set_text(txts)
        self.line_0.set_data(ds*1e9, refl)
        ax0.set_xlabel("$s$ $[nm]$")
        ax0.set_xlim(1, 1001)
        ax0.set_ylim(0, 1.01)


if __name__ == "__main__":
    """
    Called when the script runs.
    """
    fig = plt.figure()  # new figure
    hi5 = None
    plotFig = True
    if (len(sys.argv) > 1):
        args = sys.argv[1:]
        if (args[0] == "1"):
            ax0 = fig.add_subplot(111)
            hi5 = Handin5()  # the simulation class
            hi5.initTask1()  # method for task 1
        elif (args[0] == "2"):
            ax0 = fig.add_subplot(111)
            hi5 = Handin5()  # the simulation class
            hi5.initTask2()  # method for task 2
        elif (args[0] == "3"):
            ax0 = fig.add_subplot(121)
            ax1 = fig.add_subplot(122)
            hi5 = Handin5()  # the simulation class
            hi5.initTask3()  # method for task 3
        elif (args[0] == "4"):
            plotFig = False
            hi5 = Handin5()  # the simulation class
            hi5.initTask4()  # method for task 4
        elif (args[0] == "5"):
            plotFig = False
            hi5 = Handin5()  # the simulation class
            hi5.initTask5()  # method for task 5
        elif (args[0] == "6"):
            ax0 = fig.add_subplot(111)
            hi5 = Handin5()  # the simulation class
            hi5.initTask6()  # method for task 6
        else:
            sys.stdout.write("Usage: python <filename.py> <task_nbr>")
            sys.stdout.flush()
            sys.exit()
    else:
        sys.stdout.write("Usage: python <filename.py> <task_nbr>")
        sys.stdout.flush()
        sys.exit()
    if (plotFig):
        plt.show()  # plot figure
