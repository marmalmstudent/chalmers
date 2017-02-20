"""
This code is adapted for 80 columns
|------------------------------------------------------------------------------|
"""
import sys
try:
    import numpy as np
    import numpy.matlib as npm
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
except ImportError:
    sys.stderr.write("Could not execute program: " + str(sys.exc_info()[1])
                     + ".\nPlease install the module and try again")
    sys.stderr.flush()
    sys.exit()
pi = np.pi
c_0 = 2.99792458e8
mu_0 = pi*4e-7
eps_0 = 1/(mu_0*c_0**2)


class Handin4():
    """
    text
    """

    def __init__(self, nz=2000, nt=2000, mu_r=1,
                 eps_r=np.array([1], dtype=float), lbda_0=633e-9):
        """
        Class instantiation
        """
        self.txtTmpl = "Time: %.3f ps, velocity: %.0f Mm/s, imp: %.0f ohm"
        self.taskNbr = None
        self.figTxt = None
        self.E = None  # E-field
        self.H = None  # H-field
        self.dist = None  # Distance traveled
        if (len(np.shape(eps_r)) == 0):
            eps_r = np.array([eps_r], dtype=float)
        self.eps_r = eps_r
        self.eps = eps_0*eps_r  # spacial refractive index
        self.mu = mu_0*mu_r  # magnetic permeability
        self.EBounds = None  # Boundary condition for E-field
        self.nz = nz  # Number of z-points
        self.nt = nt  # Number of time-steps
        self.lbda_0 = lbda_0  # wavelength in vacuum
        self.period = self.lbda_0/c_0
        self.dt = self.period/30
        self.dz = self.lbda_0/30
        self.line = None
        self.line2 = None
        self.nItrPlot = 10  # nbr of iterations per plot

    def initField(self, envWidth, periodOffs,
                  peakField, q, m):
        """
        Calculates incident field in one position in time step m

        Parameters
        ----------
        envWidth : float
            width of envelope measured in periods (e.g. =20)
        periodOffs : float
            how many periods until the envelop reaches is max value (e.g. =
            1.5*envWidth)
        peakField : float
            max value of incident field (e.g. =1 V/m)
        q : int
            supergauss coefficient; the higher the flatter the envelope, q=1 is
            Gaussian (e.g. =4)
        m : int
            time step number

        Returns
        -------
        numpy.ndarray
            The incident field at time setp m.
        """
        f = 1/self.period
        envWidth = envWidth*self.period
        offset = periodOffs*self.period
        envelope = peakField*np.exp(-(m*self.dt-offset)**(2*q) /
                                    (envWidth/2)**(2*q))
        # incident field at time step m
        return np.cos(2*pi*f*m*self.dt)*envelope

    def DielMirrorNVec(self, nPairs, nLow, nHigh,
                       nSurround, zVec, zMirrorStart):
        """
        Calculates them number of vectors for the dielectric mirror.

        Parameters
        ----------
        nPairs : int
            ???
        nLow : int
            ???
        nHigh : int
            ???
        nSurround : int
            ???
        zVec : float
            ???
        zMirrorStart : float
            ???

        Returns
        -------
        int
            The number of vectors used in them dielectric mirror.
        """
        thickLow = self.lbda_0/(4*nLow)
        thickHigh = self.lbda_0/(4*nHigh)
        thickPair = thickLow + thickHigh
        zMirrorEnd = zMirrorStart + nPairs*thickPair
        return (((zVec - zMirrorStart) % thickPair) <= thickHigh) *\
            nHigh*(zVec > zMirrorStart)*(zVec < zMirrorEnd) +\
            (((zVec - zMirrorStart) % thickPair) >= thickHigh) *\
            nLow*(zVec > zMirrorStart)*(zVec < zMirrorEnd) +\
            nSurround*(zVec <= zMirrorStart | zVec >= zMirrorEnd)

    def diff(self, vec):
        """
        calculates the difference between index i and i+1 in a
        given vector

        Example

        Parameters
        ----------
        vec : array_like
            The array to perform operations on.

        Returns
        -------
        numpy.ndarray
            An array containing the difference between index i and i+1
        """
        if (len(np.shape(vec)) != 1):
            return None
        return vec[1:len(vec)] - vec[:len(vec)-1]

    def simData(self):
        """
        This method updates the data that is to be plotted onto the figure.

        The function yields different values for different tasks:
        Yields
        ------
        numpy.ndarray
            The E-field.
        int
            The current iteration number.
        """
        if (True):
            if (np.size(self.eps) != 1 and np.size(self.eps) != self.nz+1):
                print("Size mismatch for E-field and spacial refractive index")
                return
            for i in range(0, self.nt, self.nItrPlot):
                for j in range(0, self.nItrPlot):
                    self.E[0] = self.initField(envWidth=20,
                                               periodOffs=30,
                                               peakField=1, q=1, m=i+j) /\
                                               self.eps_r[0]
                    self.H += self.diff(self.E)*self.dt/(self.dz*self.mu)
                    self.E[1:self.nz] += self.diff(self.H) /\
                        self.eps*(self.dt/self.dz)
                yield self.E, i

    def updatefig(self, simData):
        """
        Handles setting the new image/graph and info text for the simulation.
        Used by the animation function.

        Parameters
        ----------
        simData : tuple
            simData[0] : The E-field
            simData[1] : The iteration index

        Returns
        -------
        matplotlib.lines.Line2D
            The plot containing the E-field.
        matplotlib.text.Text
            The text presenting how much time has elapsed.
        """
        if (self.taskNbr == 1 or self.taskNbr == 3):
            E, iter_idx = simData[0], simData[1]
            self.line.set_data(self.dist, E)
            dst = self.txtTmpl % (iter_idx*self.dt*1e12, self.dz/self.dt/1e6,
                                  np.sum(abs(self.E))/np.sum(np.abs(self.H)))
            self.figTxt.set_text(dst)
            return self.line, self.figTxt
        elif (self.taskNbr == 2):
            E, iter_idx = simData[0], simData[1]
            self.line.set_data(self.dist, E)
            # Note the minus sign
            self.line2.set_data(self.dist[:len(self.H)],
                                -self.E[:len(self.H)]*self.H)
            dst = self.txtTmpl % (iter_idx*self.dt*1e12, self.dz/self.dt/1e6,
                                  np.sum(abs(self.E))/np.sum(np.abs(self.H)))
            self.figTxt.set_text(dst)
            return self.line, self.figTxt

    def initializePlot(self):
        """
        Initializes the 2D plot of them beam diameter vs iteration number.
        """
        ymax = 2
        self.figTxt = ax1.text(0, 1.1*ymax, "", color="#000000")
        self.line, = ax1.plot([], [], linestyle="-", color="black")
        ax1.set_ylim(-ymax, ymax)
        ax1.set_xlabel("$z-position [\mu m]$", size=16)
        ax1.set_ylabel("$Amplitude$", size=16)
        ax1.set_xlim(0, self.dist[np.size(self.dist, axis=0)-1])

    def initTask1(self):
        """
        Initializes task 1.
        """
        self.taskNbr = 1
        self.nz = 2000
        self.nt = 2000
        self.E = np.zeros(self.nz+1)
        self.H = np.zeros(self.nz)
        self.dist = np.linspace(0, self.dz*self.nz,
                                np.size(self.E, axis=0))*1e6
        self.initializePlot()

    def initTask2(self):
        """
        Initializes task 2

        If self.nz > self.nt the wave will propagate to the right side
        and be reflected back. this is because them final point of the
        E-field is always zero so the second to last value will be mirrored.
        """
        self.taskNbr = 2
        self.nz = 1000
        self.nt = 2000
        self.E = np.zeros(self.nz+1)
        self.H = np.zeros(self.nz)
        self.dist = np.linspace(0, self.dz*self.nz,
                                np.size(self.E, axis=0))*1e6
        self.initializePlot()
        ymax = 0.01
        self.line2, = ax2.plot([], [], linestyle="-", color="black")
        ax2.set_ylim(-ymax, ymax)
        ax2.set_xlabel("$z-position [\mu m]$", size=16)
        ax2.set_ylabel("$Amplitude$", size=16)
        ax2.set_xlim(0, self.dist[np.size(self.dist, axis=0)-1])

    def initTask3(self):
        """
        Initializes task 3.
        """
        self.taskNbr = 3
        self.nz = 2000
        self.nt = 2000
        self.E = np.zeros(self.nz+1)
        self.H = np.zeros(self.nz)
        self.dist = np.linspace(0, self.dz*self.nz,
                                np.size(self.E, axis=0))*1e6
        self.initializePlot()
        self.eps_r = np.ones(self.nz)
        self.eps_r[int(self.nz/2):] = 1.33**2
        #self.eps = self.eps_r*eps_0
        #print(len(self.eps), self.nz)


if __name__ == "__main__":
    """
    Called when the script runs.
    """
    fig = plt.figure()  # new figure
    hi4 = None
    if (len(sys.argv) > 1):
        args = sys.argv[1:]
        if (args[0] == "1"):
            ax1 = fig.add_subplot(111)
            hi4 = Handin4()  # the simulation class
            hi4.initTask1()  # method for task 1
        elif (args[0] == "2"):
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122)
            hi4 = Handin4()  # the simulation class
            hi4.initTask2()  # method for task 2
        elif (args[0] == "3"):
            ax1 = fig.add_subplot(111)
            hi4 = Handin4()  # the simulation class
            hi4.initTask3()  # method for task 2
        else:
            sys.stdout.write("Usage: python <filename.py> <task_nbr>")
            sys.stdout.flush()
            sys.exit()
    else:
        sys.stdout.write("Usage: python <filename.py> <task_nbr>")
        sys.stdout.flush()
        sys.exit()
    # the function handling the the animation itself
    ani = animation.FuncAnimation(fig, hi4.updatefig, hi4.simData,
                                  interval=20, blit=False, repeat=False)
    plt.show()  # plot figure
