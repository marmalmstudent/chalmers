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

    def __init__(self, nz=30, nt=30, mu_r=1, eps_r=1):
        """
        Class instantiation
        """
        self.E = None  # E-field
        self.H = None  # H-field
        self.eps = eps_0*eps_r  # spacial refractive index
        self.mu = mu_0*mu_r  # magnetic permeability
        self.EBounds = None  # Boundary condition for E-field
        self.nz = nz  # Number of z-points
        self.nt = nt  # Number of time-steps

    def incField(self, period, envWidth, periodOffs,
                 peakField, q, dt, m):
        """
        Calculates incident field in one position in time step m

        Parameters
        ----------
        period : float
            some period i guess
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
        dt : float
            time step in FDTD
        m : int
            time step number

        Returns
        -------
        numpy.ndarray
            The incident field at time setp m.
        """
        f = 1/period
        envWidth = envWidth*period
        offset = periodOffs*period
        envelope = peakField*np.exp(-(m*dt-offset)**(2*q)/(envWidth/2)**(2*q))
        # incident field at time step m
        return np.cos(2*pi*f*m*dt)*envelope

    def DielMirrorNVec(self, nPairs, nLow, nHigh, nSurround,
                       zVec, zMirrorStart, lbda_0):
        """
        Calculates them number of vectors for the dielectric mirror.

        Parameters
        ----------
        nPairs : int
        nLow : int
        nHigh : int
        nSurround : int
        zVec : float
        zMirrorStart : float
        lbda_0 : float
            Wavelength in vacuum

        Returns
        -------
        int
            The number of vectors used in them dielectric mirror.
        """
        thickLow = lbda_0/(4*nLow)
        thickHigh = lbda_0/(4*nHigh)
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

    def yeeAlg(self):
        """
        Implements the Yee algorithm
        """
        if (np.size(self.eps) != 1 or np.size(self.eps) != self.nz+1):
            print("Size mismatch for E-field and spacial refractive index")
            return
        # set_array boundary condition for all time steps
        np.transpose(self.E)[0] = self.EBounds
        for m in range(0, self.nt):
            self.H[m] += self.diff(self.E[m])*self.dt/(self.dz*self.mu)
            self.E[m, 1:self.nz] += self.diff(self.H[m]) /\
                self.eps*(self.dt/self.dz)

    def initTask1(self):
        # e-field in space (note extra point!)
        self.E = np.zeros((self.nt, self.nz+1))
        # h-field in space
        self.H = np.zeros((self.nt, self.nz))


if __name__ == "__main__":
    """
    Called when the script runs.
    """
