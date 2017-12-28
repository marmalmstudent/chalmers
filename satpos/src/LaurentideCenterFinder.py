import numpy as np
import scipy.optimize as spo
import copy


class LaurentideCenterFinder(object):
    def __init__(self, station_models):
        self.theta = list()
        self.phi = list()
        self.rad = list()
        self.dr = list()
        for mdl in station_models:
            c = mdl.coords
            self.theta.append(c.theta().ref)
            self.phi.append(c.phi().ref)
            self.rad.append(c.rad().ref)
            self.dr.append(mdl.total_rad_move())

    def solver(self):
        elev = copy.copy(self.dr)
        elev.sort(reverse=True)
        (i0, i1) = (self.dr.index(elev[0]), self.dr.index(elev[1]))
        t_start = np.average((self.theta[i0], self.theta[i1]))
        p_start = np.average((self.phi[i0], self.phi[i1]))
        bgn = [t_start, p_start, 0, 1, 1, 3*np.pi/4, 1, 3*np.pi/4]
        return spo.fmin(lambda X: self.target(X), bgn)

    def target(self, X):
        (theta, phi) = X[0], X[1]
        (a, b, c1_r, c1_v, c2_r, c2_v) = X[2], X[3], X[4], X[5], X[6], X[7]
        c1 = c1_r*np.exp(1j*c1_v)
        c2 = c2_r*np.exp(1j*c2_v)
        R = np.average(self.rad)
        (dxi, deta) = (self.theta-theta, self.phi-phi)
        d = R*(dxi + 1j*deta)
        elevation_model = a + b*np.exp(-(np.real(c1*d)**2 + np.imag(c2*d)**2))
        w = np.sqrt(abs(np.array(self.dr)/d))
        f = np.average(abs(self.dr - elevation_model), weights=w)
        return f


class _StationData(object):
    def __init__(self, x, y, z, dr):
        self.x = x
        self.y = y
        self.z = z
        self.dr = dr
