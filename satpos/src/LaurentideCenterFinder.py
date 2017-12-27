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
            self.theta.append(c.xi.ref)
            self.phi.append(c.eta.ref)
            self.rad.append(c.zeta.ref)
            self.dr.append(mdl.total_rad_move())

    def solver(self):
        elev = copy.copy(self.dr)
        elev.sort(reverse=True)
        (i0, i1) = (self.dr.index(elev[0]), self.dr.index(elev[1]))
        t_start = np.average((self.theta[i0], self.theta[i1]))
        p_start = np.average((self.phi[i0], self.phi[i1]))
        return spo.fmin(lambda X: self.target(X), [t_start, p_start, 1, 1])

    def target(self, X):
        theta, phi = X[0], X[1]
        R = np.average(self.rad)
        (dxi, deta) = (self.theta-theta, self.phi-phi)
        d = R*np.sqrt(dxi**2 + deta**2)
        elevation_model = X[2]/np.exp(X[3]*d**2)
        w = abs(100*np.array(self.dr))**(1/16)
        f = np.average(abs(self.dr - elevation_model), weights=w)
        return f


class _StationData(object):
    def __init__(self, x, y, z, dr):
        self.x = x
        self.y = y
        self.z = z
        self.dr = dr
