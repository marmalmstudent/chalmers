import numpy as np
import scipy.optimize as spo


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
        start_idx = self.dr.index(max(self.dr))
        return spo.fmin(lambda X: self.target(X),
                        [self.theta[start_idx], self.phi[start_idx], 1, 1])

    def target(self, X):
        theta = X[0]
        phi = X[1]
        R = np.average(self.rad)
        dxi = self.theta - theta
        deta = self.phi - phi
        d = R*np.sqrt(dxi**2 + deta**2)
        f = np.average(abs(self.dr - X[2]/np.cosh(X[3]*d)),
                       weights=abs(100*np.array(self.dr))**(1/16))
        return f


class _StationData(object):
    def __init__(self, x, y, z, dr):
        self.x = x
        self.y = y
        self.z = z
        self.dr = dr
