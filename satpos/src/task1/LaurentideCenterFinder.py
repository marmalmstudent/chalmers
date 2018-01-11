from cmath import phase
from numpy import array, average, real, imag, exp, sqrt
from scipy.optimize import fmin


def euler(r, v): return r*exp(1j*v)


class LaurentideCenterFinder(object):
    def __init__(self, station_models):
        self.theta = list()
        self.phi = list()
        self.rad = list()
        self.dr = list()
        self.w = list()
        for mdl in station_models:
            c = mdl.coords
            self.theta.append(c.theta().ref)
            self.phi.append(c.phi().ref)
            self.rad.append(c.rad().ref)
            self.dr.append(mdl.total_rad_move())
            self.w.append(average(c.rad().error, weights=1/c.rad().error))
        self.R = average(self.rad)

    def solver(self):
        _dr = array(self.dr)
        w = sqrt(abs(_dr))*_dr/abs(_dr)
        sd = _SolveData(phi=average(self.phi, weights=w), theta=average(self.theta, weights=w))
        sd.stp(fmin(lambda X: self.theta_phi(sd, X[0], X[1]), [sd.phi, sd.theta], disp=True))
        return (sd.phi, sd.theta, sd.a, sd.b, sd.c1, sd.c2)

    def theta_phi(self, sd, phi, theta):
        d = self.R*((self.phi-phi) + 1j*(self.theta-theta))
        Y0 = [abs(sd.c1), phase(sd.c1), abs(sd.c2), phase(sd.c2)]
        Y = (fmin(lambda Y: self.c(sd, d, euler(Y[0], Y[1]), euler(Y[2], Y[3])), Y0, disp=True))
        sd.c1, sd.c2 = euler(Y[0], Y[1]), euler(Y[2], Y[3])
        return sd.f

    def c(self, sd, d, c1, c2):
        r = exp(-(real(c1*d)**2 + imag(c2*d)**2))
        Y0 = [sd.a, sd.b]
        Y = (fmin(lambda Y: self.a_b(sd, Y[0] + Y[1]*r), Y0, disp=False))
        sd.a, sd.b = Y[0], Y[1]
        return sd.f

    def a_b(self, sd, f):
        sd.f = average(abs(self.dr - f))
        return sd.f


class _StationData(object):
    def __init__(self, x, y, z, dr):
        self.x = x
        self.y = y
        self.z = z
        self.dr = dr


class _SolveData(object):
    def __init__(self, phi, theta, c1=1, c2=1, a=0, b=1):
        self.theta = theta
        self.phi = phi
        self.c1 = c1
        self.c2 = c2
        self.a = a
        self.b = b

    def stp(self, X):
        self.phi, self.theta = X

    def sab(self, X):
        self.a, self.b = X

    def sc(self, X):
        self.c1 = euler(X[0], X[1])
        self.c2 = euler(X[2], X[3])
