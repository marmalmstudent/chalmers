import math
import csv
from Station import StationModel
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d
import numpy as np


pi = np.pi
factor = 1e-6
R = 6371000 * factor


def csvl(fpath):
    out = dict()
    with open(fpath, "r") as f:
        for l in csv.reader(f, delimiter=","):
            out[l[0]] = coord_dict_from_csv_list(l[1:])
    return out


def coord_dict_from_csv_list(l):
    out = dict()
    for e in l:
        k, v = e.split(":")
        out[k] = float(v)
    return out


def tpr_to_xyz(t, p, r):
    x = r*np.cos(t)*np.cos(p)
    y = r*np.cos(t)*np.sin(p)
    z = r*np.sin(t)
    return (x, y, z)


class StationManager(object):
    def __init__(self, respath, stations, years):
        self.sm = dict()
        for stn in stations:
            self.sm[stn] = StationModel(stn, respath, years)

    def print_station_data(self):
        for mdl in self.sm.values():
            tom = mdl.total_lon_move()*100
            tam = mdl.total_lat_move()*100
            trm = mdl.total_rad_move()*100
            print(("Total move for > %5s: Lon:%+6.2f cm, "
                   + "Lat:%+6.2f cm, Rad:%+6.2f cm, "
                   + "Lat-Lon:%+6.2f cm, Lat-Lon-Rad:%+6.2f cm")
                  % (mdl.acro, tom, tam, trm,
                     math.sqrt(tom**2 + tam**2),
                     math.sqrt(tom**2 + tam**2 + trm**2)))

    def load_station_ref_coords(self, respath, fname):
        coords_data = csvl("%s/%s" % (respath, fname))
        for acro, mdl in self.sm.items():
            cdict = coords_data[acro.lower()]
            coords = mdl.coords
            coords.set_ref_pos(cdict["x"] * factor,
                               cdict["y"] * factor,
                               cdict["z"] * factor)

    def make_spherical(self):
        for mdl in self.sm.values():
            mdl.coords.cart_to_sphere()


class StationManagerView(object):
    def __init__(self, sm_mgr):
        self.sm_mgr = sm_mgr
        self.fig = plt.figure()
        # self.ax1 = self.fig.add_subplot(111, projection='3d')
        self.ax1 = self.fig.add_subplot(111)
        self.ax1.set_xlabel("Latitude (North-South)")
        self.ax1.set_ylabel("Longitude (West-East)")
        # self.ax1.view_init(elev=-90, azim=0)
        self.settings = {"linestyle": "None"}

    def present(self):
        plt.show()

    def plot_lon_lat(self):
        for mdl in self.sm_mgr.sm.values():
            (xi, eta, zeta) = (mdl.coords.xi, mdl.coords.eta, mdl.coords.zeta)
            msize = max(1, 200*(zeta.val[-1]-zeta.val[0]))
            self.ax1.plot([eta.ref], [xi.ref], marker="o", markersize=msize)
        self.ax1.set_aspect("equal")

    def plot_laurentide_center(self, theta, phi, a, b):
        rad = list()
        for mdl in self.sm_mgr.sm.values():
            rad.append(mdl.coords.zeta.ref)
        R = np.average(rad)
        (U, V) = np.meshgrid(np.linspace(phi - pi/4, phi + pi/4, 101),
                             np.linspace(theta - pi/4, theta + pi/4, 101))
        dphi = U - phi
        dtheta = V - theta
        d = R*np.sqrt(dphi**2 + dtheta**2)
        f = a/np.cosh(b*d)
        self.ax1.plot([phi], [theta], marker="x", markersize=10)
        # self.ax1.plot_surface(U, V, f)
        self.ax1.contour(U, V, f, 10)
