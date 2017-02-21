"""
This code is adapted for 80 columns
|------------------------------------------------------------------------------|
"""
import sys
try:
    import numpy as np
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
eta_0 = 377


class Handin4():
    """
    text
    """

    def __init__(self, nz=2000, nt=2000, mu_r=1,
                 eps_r=np.array([1], dtype=float), lbda_0=633e-9):
        """
        Class instantiation
        """
        self.txtTmpl = " "
        self.taskNbr = None
        self.figTxt = None  # the figure text object
        self.E = None  # E-field
        self.E_init = None  # initial E-field
        self.H = None  # H-field
        self.H_init = None  # initial H-field
        self.dist = None  # Distance traveled
        if (len(np.shape(eps_r)) == 0):
            eps_r = np.array([eps_r], dtype=float)
        self.eps_r = eps_r  # relative electic permittivity
        self.eps = eps_0*eps_r  # electric permittivity
        self.mu = mu_0*mu_r  # magnetic permeability
        self.nz = nz  # Number of z-points
        self.nt = nt  # Number of time-steps
        self.n = np.sqrt(self.eps_r)  # refractive index
        self.lbda_0 = lbda_0  # wavelength in vacuum
        self.period = self.lbda_0/c_0  # period of the wave [s]
        self.dt = self.period/300  # time step
        self.dz = self.lbda_0/300  # spacial step
        self.line = None  # E-field plot line
        self.line2 = None  # Poynting field plot line
        self.nItrPlot = 10  # nbr of iterations per plot
        self.pulseWidth = 10  # pulse width in periods
        self.pulseOffs = 15  # pulse offset in periods
        self.sigma = np.array([0.0])  # electric conductivity
        self.mirrBegIdx = 1  # start index (z-direction) of mirror
        self.mirrEndIdx = 1  # end index (z-direction) of mirror

    def setNT(self, newNT):
        """
        Sets the number of time steps and adjusts the number of steps
        to loop between plots to make the program run faster.

        Parameters
        ----------
        newNT : int
            The new numer of time steps
        """
        self.nt = newNT
        self.nItrPlot = int(newNT/200)

    def setNZ(self, newNZ):
        """
        Sets the number of spacial steps and adjusts the time steps
        and the spacial steps to make the program run faster.

        Parameters
        ----------
        newNZ : int
            The new numer of time steps
        """
        self.nz = newNZ
        self.dt = self.period/(0.015*newNZ)
        self.dz = self.lbda_0/(0.015*newNZ)

    def setDZ(self, newDZ):
        """
        Sets the spacial step and adjusts the number of spacial steps
        and the time step to make the program run faster.

        Parameters
        ----------
        newDZ : int
            The new numer of time steps
        """
        self.dz = newDZ
        self.nz = int(self.lbda_0/(0.015*self.dz))
        self.dt = self.period/(0.015*self.nz)

    def initField(self, peakField, q, m):
        """
        Calculates incident field in one position in time step m

        Parameters
        ----------
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
        envWidth = self.pulseWidth*self.period
        offset = self.pulseOffs*self.period
        envelope = peakField*np.exp(-(m*self.dt-offset)**(2*q) /
                                    (envWidth/2)**(2*q))
        # incident field at time step m
        return np.cos(2*pi*f*m*self.dt)*envelope

    def applyDielMirr(self, zStartIdx, nHigh, nLow, nPairs):
        """
        Applies the dielectric mirror to the relative permittivity
        vector.

        Parameters
        ----------
        zStartIdx : int
            The index in the eps_r vector where the mirror starts
        nHigh : float
            The high refractive index
        nLow : float
            The low refractive index
        nPairs : int
            The number of pairs of log and high refractive index segments
        """
        for i in range(0, nPairs):
            segmLen = int(self.lbda_0/(4*nHigh*self.dz))
            self.eps_r[zStartIdx:zStartIdx+segmLen] = nHigh**2
            zStartIdx += segmLen
            segmLen = int(self.lbda_0/(4*nLow*self.dz))
            self.eps_r[zStartIdx:zStartIdx+segmLen] = nLow**2
            zStartIdx += segmLen
        return zStartIdx

    def calcZSalmpDist(self, n1, n2, nSamplMin, lbdaFrac):
        """
        Calculates the spacial sampling distance to make fraction of lambda
        equal an integer number of sampling points for both indices.

        Parameters
        ----------
        n1 : double
            The first refractive index
        n2 : double
            The second refractive index
        nSamplMin : int
            The minimum amount of sampling points per wavelength
        lbdaFrac : float
            The fraction of the wavelength (in the material) that will be
            the width of the material.
        """
        nSamplPtsMin = np.ceil(nSamplMin*lbdaFrac)
        # n1 will have fewer points
        samplDefault = np.array([np.ceil(nSamplPtsMin/np.min((n1, n2)))],
                                dtype=np.uint32)
        nIters = 0
        while(True):
            if (((samplDefault/n1) % 1 == 0) and ((samplDefault/n2) % 1 == 0)):
                break
            samplDefault += 1
            nIters += 1
            if (nIters > 1e3):
                print("Could not find a reasonable sample distance")
                return
        self.setDZ(self.lbda_0*lbdaFrac/samplDefault)

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
        if (np.size(self.eps) != 1 and np.size(self.eps) != self.nz):
            print("Size mismatch for E-field and spacial refractive index")
            return
        # implements the Yee-algorithm
        for i in range(0, self.nt, self.nItrPlot):
            for j in range(0, self.nItrPlot):
                self.E[0] = self.initField(peakField=1, q=4, m=i+j) /\
                                           self.eps_r[0]
                self.H += self.diff(self.E)*self.dt/(self.dz*self.mu)
                self.E[1:self.nz] += (self.dt/self.eps[:self.nz-1]) *\
                                     (self.diff(self.H)/self.dz -
                                      self.sigma[:self.nz-1]*self.E[1:self.nz])
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
        if (self.taskNbr == 1):
            """
            In task 1 the electric field is plotted along with information
            about how many iterations have been run, how far the wave has
            propagated (in ps) and the velocity (in Mm/s)
            """
            E, iter_idx = simData[0], simData[1]
            self.line.set_data(self.dist[:len(self.H)], E[:len(self.H)])
            dst = self.txtTmpl % (iter_idx, iter_idx*self.dt*1e12,
                                  self.dz/self.dt/1e6)
            self.figTxt.set_text(dst)
            return self.line, self.figTxt
        elif (self.taskNbr == 2):
            """
            In task 2 the electric field and the poynting vector is plotted
            in two separete figures. Also some information about how many
            iterations that have been run, how far the wave has propagated
            (in ps) and the wave impedance.
            """
            E, iter_idx = simData[0], simData[1]
            self.line.set_data(self.dist, E)
            # Note the minus sign
            self.line2.set_data(self.dist[:len(self.H)],
                                -self.E[:len(self.H)]*self.H)
            waveImp = 0
            if (np.sum(np.abs(self.H)) != 0):
                waveImp = np.sum(abs(self.E))/np.sum(np.abs(self.H))
            dst = self.txtTmpl % (iter_idx, iter_idx*self.dt*1e12, waveImp)
            self.figTxt.set_text(dst)
            return self.line, self.figTxt
        elif (self.taskNbr == 3 or self.taskNbr == 4):
            """
            In task 3 and 4 the poynting vector is plotted along with some
            information about how many iterations that have been run, how far
            the wave has propagated (in ps), the reflected power and the
            transmitted power.
            """
            E, iter_idx = simData[0], simData[1]
            self.line2.set_data(self.dist[:len(self.H)],
                                -self.E[:len(self.H)]*self.H)
            rx = np.sum(np.abs(self.E[:self.mirrBegIdx] *
                               self.H[:self.mirrBegIdx]))*self.dz
            tx = np.sum(np.abs(self.E[self.mirrEndIdx:len(self.H)] *
                               self.H[self.mirrEndIdx:]))*self.dz
            dst = self.txtTmpl % (iter_idx, iter_idx*self.dt*1e12,
                                  rx/(rx+tx*self.n[self.mirrEndIdx])*100,
                                  tx*self.n[self.mirrEndIdx] /
                                  (rx+tx*self.n[self.mirrEndIdx])*100)
            self.figTxt.set_text(dst)
            return self.line, self.figTxt
        elif (self.taskNbr == 5):
            """
            In task 5 the poynting vector is plotted along with some
            information about how many iterations that have been run, how far
            the wave has propagated (in ps), the reflected power and the
            transmitted power.
            """
            E, iter_idx = simData[0], simData[1]
            self.line2.set_data(self.dist[:len(self.H)],
                                -self.E[:len(self.H)]*self.H)
            # these are the amplitudes
            rx = np.sum(np.abs(self.E[:self.mirrBegIdx] *
                               self.H[:self.mirrBegIdx]))*self.dz
            tx = np.sum(np.abs(self.E[self.mirrEndIdx:len(self.H)] *
                               self.H[self.mirrEndIdx:]))*self.dz
            # add 10^-40 to avoid division by zero.
            tot = np.sum(np.abs(self.E[:len(self.H)]*self.H))*self.dz + 1e-40
            dst = self.txtTmpl % (iter_idx, iter_idx*self.dt*1e12,
                                  rx/tot*100,
                                  tx/tot*100)
            self.figTxt.set_text(dst)
            return self.line, self.figTxt
        elif (self.taskNbr == 6):
            """
            In task 6 the poynting vector is plotted along with some
            information about how many iterations that have been run, how far
            the wave has propagated (in ps), the reflected power (as a
            percentage of the incident power).
            """
            E, iter_idx = simData[0], simData[1]
            self.line2.set_data(self.dist[:len(self.H)],
                                -self.E[:len(self.H)]*self.H)
            # these are the amplitudes
            rx = np.sum(np.abs(self.E[:self.mirrBegIdx] *
                               self.H[:self.mirrBegIdx]))*self.dz
            # add 10^-40 to avoid division by zero.
            tot = np.sum(np.abs(self.E[:len(self.H)]*self.H))*self.dz + 1e-40
            dst = self.txtTmpl % (iter_idx, iter_idx*self.dt*1e12,
                                  rx/tot*100)
            self.figTxt.set_text(dst)
            return self.line, self.figTxt

    def initializeEPlot(self):
        """
        Initializes the 2D plot of the E-field amplitude vs propagated
        distance.
        """
        ymax = 2
        self.figTxt = ax1.text(0, 1.1*ymax, "", color="#000000")
        self.line, = ax1.plot([], [], linestyle="-", color="black")
        ax1.set_ylim(-ymax, ymax)
        ax1.set_xlabel("$z-position [\mu m]$", size=16)
        ax1.set_ylabel("$Amplitude$", size=16)
        ax1.set_xlim(0, self.dist[np.size(self.dist, axis=0)-1])

    def initializePPlot(self, dispTxt):
        """
        Initializes the 2D plot of the poynting vecotr vs propagated distance.
        """
        ymax = 0.003
        if (dispTxt):
            self.figTxt = ax2.text(0, 1.1*ymax, "", color="#000000")
        self.line2, = ax2.plot([], [], linestyle="-", color="black")
        ax2.set_ylim(-ymax, ymax)
        ax2.set_xlabel("$z-position [\mu m]$", size=16)
        ax2.set_ylabel("$Amplitude$", size=16)
        ax2.set_xlim(0, self.dist[np.size(self.dist, axis=0)-1])

    def initTask1(self):
        """
        Initializes task 1.
        """
        self.txtTmpl = "Iterations: %.0f, Time: %.3f ps, velocity: %.0f Mm/s"
        self.taskNbr = 1
        self.setNT(2000)
        self.setNZ(2000)
        self.E = np.zeros(self.nz+1)
        self.H = np.zeros(self.nz)
        self.dist = np.linspace(0, self.dz*self.nz,
                                np.size(self.E, axis=0))*1e6
        self.initializeEPlot()

    def initTask2(self):
        """
        Initializes task 2

        If self.nz > self.nt the wave will propagate to the right side
        and be reflected back. this is because them final point of the
        E-field is always zero so the second to last value will be mirrored.
        """
        self.txtTmpl = "Iterations: %.0f, Time: %.3f ps, imp: %.0f ohm"
        self.taskNbr = 2
        self.setNT(20000)
        self.setNZ(11500)
        self.E = np.zeros(self.nz+1)
        self.H = np.zeros(self.nz)
        self.dist = np.linspace(0, self.dz*self.nz,
                                np.size(self.E, axis=0))*1e6
        self.initializeEPlot()
        self.initializePPlot(False)

    def initTask3(self):
        """
        Initializes task 3.
        """
        self.txtTmpl = "Iterations: %.0f, Time: %.3f ps, RXPWR: %.1f " +\
                       "percent, TXPWR: %.1f percent"
        self.taskNbr = 3
        self.setNT(12000)
        self.setNZ(20000)
        self.E_init = self.initField(peakField=1, q=4, m=np.arange(0, self.nt))
        self.H_init = self.diff(self.E_init)*self.dt/(self.dz*self.mu)
        self.E = np.zeros(self.nz+1)
        self.H = np.zeros(self.nz)
        self.dist = np.linspace(0, self.dz*self.nz,
                                np.size(self.E, axis=0))*1e6
        self.initializePPlot(True)
        self.eps_r = np.ones(self.nz)
        self.mirrBegIdx = int(self.nz/4)
        self.mirrEndIdx = int(self.nz/4)
        self.eps_r[self.mirrBegIdx:] = 1.5**2
        self.n = np.sqrt(self.eps_r)
        self.eps = self.eps_r*eps_0

    def initTask4(self):
        """
        Initializes task 4
        More samples makes the enegery coservation -> 100 %
        """
        self.txtTmpl = "Iterations: %.0f, Time: %.3f ps, RXPWR: %.1f " +\
                       "percent, TXPWR: %.1f percent"
        self.taskNbr = 4
        self.setNT(12000)
        self.setNZ(20000)
        self.E_init = self.initField(peakField=1, q=4, m=np.arange(0, self.nt))
        self.H_init = self.diff(self.E_init)*self.dt/(self.dz*self.mu)
        self.E = np.zeros(self.nz+1)
        self.H = np.zeros(self.nz)
        self.dist = np.linspace(0, self.dz*self.nz,
                                np.size(self.E, axis=0))*1e6
        self.initializePPlot(True)
        n_sub = 1.5
        n_coat = np.sqrt(n_sub)
        nARSampl = self.lbda_0/(4*n_coat*self.dz)
        print(nARSampl)
        self.mirrBegIdx = int(self.nz/4)
        self.mirrEndIdx = self.mirrBegIdx + int(nARSampl)
        self.eps_r = np.ones(self.nz)
        self.eps_r[self.mirrBegIdx:self.mirrEndIdx] = n_coat**2
        self.eps_r[self.mirrEndIdx:] = n_sub**2
        self.n = np.sqrt(self.eps_r)
        self.eps = self.eps_r*eps_0

    def initTask5(self):
        """
        Initializes task 5
        longer pulses yields values closer to them theoretical.
        """
        self.txtTmpl = "Iterations: %.0f, Time: %.3f ps, RXPWR: %.1f " +\
                       "percent, TXPWR: %.1f percent"
        self.taskNbr = 5
        nHigh = 1.7
        nLow = 1.5
        self.pulseWidth = 20
        self.pulseOffs = 30
        self.calcZSalmpDist(nLow, nHigh, 200, 0.25)
        self.nz = 2*self.nz
        self.setNT(int(self.nz*0.75))
        self.E_init = self.initField(peakField=1, q=4, m=np.arange(0, self.nt))
        self.H_init = self.diff(self.E_init)*self.dt/(self.dz*self.mu)
        self.E = np.zeros(self.nz+1)
        self.H = np.zeros(self.nz)
        self.dist = np.linspace(0, self.dz*self.nz,
                                np.size(self.E, axis=0))*1e6
        self.initializePPlot(True)
        self.eps_r = np.ones(self.nz)
        self.mirrBegIdx = int(self.nz/3)
        self.mirrEndIdx = self.applyDielMirr(self.mirrBegIdx, nHigh, nLow, 20)
        self.n = np.sqrt(self.eps_r)
        self.eps = self.eps_r*eps_0

    def initTask6(self):
        """
        Initializes task 6
        """
        self.txtTmpl = "Iterations: %.0f, Time: %.3f ps, RXPWR: %.1f percent"
        self.taskNbr = 6
        self.pulseWidth = 30
        self.pulseOffs = 40
        self.setNT(13000)
        self.setNZ(8000)
        self.nz = 2*self.nz
        self.E_init = self.initField(peakField=1, q=4, m=np.arange(0, self.nt))
        self.H_init = self.diff(self.E_init)*self.dt/(self.dz*self.mu)
        self.E = np.zeros(self.nz+1)
        self.H = np.zeros(self.nz)
        self.dist = np.linspace(0, self.dz*self.nz,
                                np.size(self.E, axis=0))*1e6
        self.initializePPlot(True)
        self.eps_r = np.ones(self.nz)
        self.sigma = np.zeros(self.nz)
        n_sub = 1.5
        self.mirrBegIdx = int(self.nz/3)
        self.mirrEndIdx = int(self.nz/3)
        self.eps_r[self.mirrBegIdx:] = n_sub**2
        self.sigma[self.mirrBegIdx:] = -1e2
        self.n = np.sqrt(self.eps_r)
        self.eps = self.eps_r*eps_0


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
            ax2 = fig.add_subplot(111)
            hi4 = Handin4()  # the simulation class
            hi4.initTask3()  # method for task 3
        elif (args[0] == "4"):
            ax2 = fig.add_subplot(111)
            hi4 = Handin4()  # the simulation class
            hi4.initTask4()  # method for task 4
        elif (args[0] == "5"):
            ax2 = fig.add_subplot(111)
            hi4 = Handin4()
            hi4.initTask5()
        elif (args[0] == "6"):
            ax2 = fig.add_subplot(111)
            hi4 = Handin4()  # the simulation class
            hi4.initTask6()  # method for task 4
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
