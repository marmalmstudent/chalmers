"""
This code is adapted for 80 columns
|------------------------------------------------------------------------------|
"""
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys
pi = np.pi


def loadHolo():
    """
    loads the hologram provided for the task123.
    """
    mat = scipy.io.loadmat("Hologram_transmission_function.mat")
    return mat["T_hologram"]


def fft2c(x):
    """
    shifted fourier transform,
    """
    return np.fft.fftshift(np.fft.fft2(np.fft.fftshift(x)))


def ifft2c(x):
    """
    shifted inverse fourier transform.
    """
    return np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(x)))


def getCoordVect(nPoints, samplDist):
    """
    A static method that returns an array of nPoints with a spacing of
    sampleDist starting from -nPoints/2 to nPoints/2.
    """
    return np.arange(-nPoints/2*samplDist, (nPoints/2)*samplDist, samplDist)


def TSM(E1, a, b, L, lbda_0, n):
    """
    This function implements the conventional Huygens-Fresnel method (TSM) for
    free space propagation.
    """
    N = np.size(E1, axis=0)  # Matrix size (one side of the square matrix E1)
    lambda_medium = lbda_0/n
    k = 2*pi/lambda_medium

    # Plane 1 coordinates:
    xvect = getCoordVect(N, a)
    yvect = getCoordVect(N, a)
    [xmat, ymat] = np.meshgrid(xvect, yvect)

    # Plane 2 coordinates
    uvect = getCoordVect(N, b)
    vvect = getCoordVect(N, b)
    [umat, vmat] = np.meshgrid(uvect, vvect)

    # Distance to dummy plane
    L1 = L*a/(a-b)  # plane 1 to dummy plane
    L2 = L*b/(a-b)  # plane 2 to dummy plane
    # Sampling distance in Dummy plane
    c = lambda_medium/(N*a)*L1

    # Dummy plane coordinates
    xiVect = getCoordVect(N, c)
    etaVect = getCoordVect(N, c)
    [xiMat, etaMat] = np.meshgrid(xiVect, etaVect)

    # Calculating the field in Plane 2
    fctr0 = a**2 / b**2 * L2/L1 * np.exp(1j*k*(L1-L2))
    fctr1 = np.exp(-1j*k *
                   (umat**2 + vmat**2) / (2*L2))
    prefactor = np.exp(
        1j*k * (xiMat**2 + etaMat**2) / 2 * (1/L1 - 1/L2))
    fft_part = fft2c(
        E1*np.exp(1j*k * (xmat**2 + ymat**2) / (2*L1)))
    fctr2 = ifft2c(prefactor * fft_part)
    E2 = fctr0*fctr1*fctr2
    return E2, L1, L2, uvect, vvect, umat, vmat


class handin1():
    """
    This class contains the necessary functions to run the simulations.
    For simplicity and for the sake of not having to pass too many paramters
    arround most variables are class variables.
    """
    
    def __init__(self):
        """
        class instantiation
        """
        self.txt_tmpl = "Distance: %.4f m"
        self.img_txt = plt.text(0.0, -5.0, "")
        self.N = 0x1 << 8  # number fo sample points
        self.lbda_0 = 635e-9  # wavelength in vacuum
        self.n = 1.0  # refractive index
        self.k = 2*pi*self.n/self.lbda_0
        self.a_anim = 100e-7  # sample distance at animation plane 1
        self.b_anim = 99e-7  # sample distance at animation plane 2
        self.E1 = None  # the starting field for the animation
        self.img = None  # the animation image
        self.doe = None  # difractive optical element
        self.lense = None
        self.f_lense = None  # focal length of lense
        self.L_bounds = [-1, 1]
        self.task_nbr = 0  # If a circle to block the bright spot

    def mainFun(self, L):
        """
        The purpose of this function is to calculcate the next field
        during the animation process. If necessary it also blocks
        the bright spot in the center for the hologram task123.
        """
        E2, L1, L2, uvect, vvect, umat, vmat = TSM(self.E1, self.a_anim,
                                                   self.b_anim, L,
                                                   self.lbda_0, self.n)
        if (self.task_nbr == 4):
            [xmat, ymat] = np.meshgrid(getCoordVect(self.N, self.b_anim),
                                       getCoordVect(self.N, self.b_anim))
            rmat = np.sqrt(xmat**2 + ymat**2)
            maskRadius = 40e-5
            E2 = E2*(rmat > maskRadius)
        return E2

    def simData(self):
        """
        generates the data needed in each simulation step. This includes
        the distance and the matrix.
        """
        if (self.task_nbr > 0 and self.task_nbr < 4):
            # L = [2e-1]
            L = np.logspace(-8.0, 0.0, 100)
        elif (self.task_nbr == 4):
            # L = [1.5]
            L = np.linspace(1.4, 1.6, 100)
        for i in L:
            yield self.mainFun(i), i
        return

    def updatefig(self, simData):
        """
        Handles setting the new image and distance text for the simulation.
        """
        E2, L = simData[0], simData[1]
        E2_plt = np.array(np.abs(E2)/np.max(np.abs(E2)), dtype=np.float64)**2
        dst = self.txt_tmpl % L
        self.img_txt.set_text(dst)
        self.img.set_array(E2_plt)
        return self.img, self.img_txt

    def task123(self):
        """
        Initializes the task123 computations
        """
        self.task_nbr = 1
        [xmat, ymat] = np.meshgrid(getCoordVect(self.N, self.a_anim),
                                   getCoordVect(self.N, self.a_anim))
        rmat = np.sqrt(xmat**2 + ymat**2)
        maskRadius = 400e-6
        bool_mat = (rmat > maskRadius)
        self.E1 = np.ones((self.N, self.N))*bool_mat
        L_start = 1e-8
        E2 = self.mainFun(L_start)
        E2_plt = np.array(np.abs(E2)/np.max(np.abs(E2)), dtype=np.float64)**2
        self.img = plt.imshow(E2_plt, animated=True, cmap=plt.get_cmap('gray'),
                              interpolation='quadric', vmin=0, vmax=1)
        return

    def task4(self):
        """
        Initializes the task4 computations
        """
        """ The hologram i created myself
        self.doe = np.loadtxt('DOE.txt', comments='#',
                              delimiter=None, converters=None, skiprows=5,
                              usecols=None, unpack=False, ndmin=0)
        """
        self.task_nbr = 4
        self.doe = loadHolo()  # the hologram provided by for the task
        self.E1 = self.doe  # starting field
        E2_tmp = self.mainFun(1e-1)  # field just before the lens
        self.a_anim = 80e-6
        self.b_anim = 78e-6
        [xmat, ymat] = np.meshgrid(getCoordVect(self.N, self.b_anim),
                                   getCoordVect(self.N, self.b_anim))
        # lens transfer function.
        self.f_lense = 1.5
        self.lense = np.exp(-1j*self.k*(xmat**2 + ymat**2)/(2*self.f_lense))
        L_start = 1.4
        self.E1 = E2_tmp*self.lense  # starting field for TSM (after the lens)
        # new sampling distance after the lens
        self.a_anim = self.b_anim
        # new sampling distance in the plane we look at
        self.b_anim = self.lbda_0/(self.n*self.N*self.a_anim)*self.f_lense
        E2 = self.mainFun(L_start)  # field in the plane we look at
        E2_plt = np.array(np.abs(E2)/np.max(np.abs(E2)), dtype=np.float64)**2
        self.img = plt.imshow(E2_plt, animated=True, cmap=plt.get_cmap('gray'),
                              interpolation='quadric', vmin=0, vmax=1)
        return


if (__name__ == "__main__"):
    """
    This functuin is called when the program executes
    """
    fig = plt.figure()  # new figure
    hi1 = None  # the simulation class
    if (len(sys.argv) > 1):
        hi1 = handin1()  # the simulation class
        args = sys.argv[1:]
        if (args[0] == "1" or args[0] == "2" or args[0] == "3"):
            hi1.task123()  # method for task 1, 2 and 3
        elif (args[0] == "4"):
            hi1.task4()  # method for task 4
        else:
            sys.stdout.write("Usage: python <filename.py> <task_nbr>")
            sys.stdout.flush()
            sys.exit()
    else:
        sys.stdout.write("Usage: python <filename.py> <task_nbr>")
        sys.stdout.flush()
        sys.exit()
    # the function handling the the animation itself
    ani = animation.FuncAnimation(fig, hi1.updatefig, hi1.simData,
                                  interval=50, blit=False, repeat=True)
    plt.show()  # plot figure
