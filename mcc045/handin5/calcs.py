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
        self.lbda_0 = None
        self.lbdaZero = 550e-9  # center wavelength
        self.lbdaMax = 3*self.lbdaZero
        self.lbdaMin = self.lbdaZero/3

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
        self.figTxt = ax1.text(0, 1.1*ymax, "", color="#000000")
        self.line, = ax1.plot([], [], linestyle="-", color="black")
        ax1.set_ylim(0, ymax)
        ax1.set_xlabel("$wavelength [nm]$", size=16)
        ax1.set_ylabel("$Reflectance$", size=16)
        ax1.set_xlim(self.lbdaMin, self.lbdaMax)

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

    def initTask1(self):
        self.taskNbr = 1
        self.lbdaZero = 550e-9  # center wavelength
        self.lbdaMax = 50*self.lbdaZero
        self.lbdaMin = self.lbdaZero/50
        nPoints = int((self.lbdaMax-self.lbdaMin)*1e9)*100
        self.lbda_0 = np.linspace(self.lbdaMin, self.lbdaMax, nPoints,
                                  dtype=np.float64)
        self.initializePlot()
        nSubst = 1.5
        nAR = np.sqrt(nSubst)
        d = self.lbdaZero/(4*nAR)  # AR thickness
        nVac = np.ones(nPoints, dtype=np.complex128)
        nsv = nSubst*nVac  # substrate refr idx
        narv = nAR*nVac  # anti-reflection coating refr idx
        matrices = [self.BorderTX(nVac, narv),
                    self.PropTX(d, narv, self.lbda_0),
                    self.BorderTX(narv, nsv)]
        dims = np.array([2, 2, nPoints], dtype=np.int32)
        matTot = np.array([[np.ones(nPoints), np.zeros(nPoints)],
                           [np.zeros(nPoints), np.ones(nPoints)]],
                          dtype=np.complex128)
        for i in range(len(matrices), 0, -1):
            matTot = ha5utils.vecMatDot(dims, matrices[i-1], matTot)
        refl = np.abs(matTot[1, 0]/matTot[1, 1])
        self.line.set_data(self.lbda_0, refl)
        ax1.set_xscale('log')


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
