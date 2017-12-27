import math
import csv
from Station import StationModel
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d
import numpy as np


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


"""
container for all active stations.
handles loading position for all stations into staiton xi,eta,zeta

fin station with highest zeta movement
subtract eta-xi movement of station with highest zeta movement from
all stations.
create line between station with highest zeta movement and all other
stations and project those stations eta-xi movement to tha
"""


class StationManager(object):
    def __init__(self, respath, stations, years):
        self.sm = dict()
        for stn in stations:
            self.sm[stn] = StationModel(stn, respath, years)

    def print_station_data(self):
        for mdl in self.sm.values():
            print(("Total move for > %5s: Lon:%+6.3f m, "
                   + "Lat:%+6.3f m, Rad:%+6.3f m, "
                   + "Lat-Lon:%+6.3f m, Lat-Lon-Rad:%+6.3f m")
                  % (mdl.acro, mdl.total_lon_move(),
                     mdl.total_lat_move(), mdl.total_rad_move(),
                     math.sqrt(mdl.total_lon_move()**2
                               + mdl.total_lat_move()**2),
                     math.sqrt(mdl.total_lon_move()**2
                               + mdl.total_lat_move()**2
                               + mdl.total_rad_move()**2)))

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
        self.ax1 = self.fig.add_subplot(111, projection='3d')
        self.ax1.set_xlabel("Latitude (North-South)")
        self.ax1.set_ylabel("Longitude (West-East)")
        self.ax1.view_init(elev=-90, azim=0)
        self.settings = {"linestyle": "None"}

    def present(self):
        plt.show()

    def _plot_lon_lat(self):
        """
        U, V = np.meshgrid(np.linspace(7*np.pi/6, 11*np.pi/6, 33),
                           np.linspace(0, np.pi/2, 21))
        """
        U, V = np.meshgrid(np.linspace(0, 2*np.pi, 33),
                           np.linspace(0, np.pi, 21))
        x = R*np.cos(U)*np.sin(V)
        y = R*np.sin(U)*np.sin(V)
        z = R*np.cos(V)
        # self.ax1.plot_surface(x, y, z)
        for mdl in self.sm_mgr.sm.values():
            xi = mdl.coords.xi
            eta = mdl.coords.eta
            zeta = mdl.coords.zeta
            (_x, _y, _z) = (xi.val, eta.val, zeta.val)
            if mdl.coords._spherical:
                (_x, _y, _z) = tpr_to_xyz(xi.val, eta.val, zeta.val)
            self.ax1.plot([_x[0], _x[0]], [_y[0], _y[0]], [_z[0], _z[0]],
                          marker="x")
            self.ax1.plot(_x, _y, _z)
        self.ax1.set_aspect("equal")

    def plot_lon_lat(self):
        for mdl in self.sm_mgr.sm.values():
            xi = mdl.coords.xi
            eta = mdl.coords.eta
            zeta = mdl.coords.zeta
            (_x, _y, _z) = (xi.val, eta.val, zeta.val - zeta.val[0])
            self.ax1.plot([_x[0], _x[0]], [_y[0], _y[0]], [_z[0], _z[0]],
                          marker="x")
            self.ax1.plot(_x, _y, _z)
        self.ax1.set_aspect("equal")

    def plot_laurentide_center(self, theta, phi, a, b):
        rad = list()
        for mdl in self.sm_mgr.sm.values():
            rad.append(mdl.coords.zeta.ref)
        R = np.average(rad)
        U, V = np.meshgrid(np.linspace(theta - np.pi/4, theta + np.pi/4, 101),
                           np.linspace(phi - np.pi/4, phi + np.pi/4, 101))
        dtheta = U - theta
        dphi = V - phi
        d = R*np.sqrt(dtheta**2 + dphi**2)
        f = a/np.cosh(b*d)
        self.ax1.contour(U, V, f)
        # self.ax1.plot_surface(U, V, f)
        """
        (_x, _y, _z) = (theta, phi, 0)
        self.ax1.plot([_x, _x], [_y, _y], [_z, _z], marker="o")
        """
