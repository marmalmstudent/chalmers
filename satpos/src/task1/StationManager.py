import math
import csv
from .Station import StationModel
from .MapPlot import MapPlot
import matplotlib.pyplot as plt
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
    def __init__(self, respath, stations, years, yearly_average, correct_errors):
        self.sm = dict()
        for stn in stations:
            self.sm[stn] = StationModel(stn, respath, years, yearly_average, correct_errors)

    def print_station_data(self):
        for mdl in self.sm.values():
            tom = mdl.total_lon_move()*1000
            tam = mdl.total_lat_move()*1000
            trm = mdl.total_rad_move()*1000
            print(("Total move for    > %5s: Lon:%+6.2f mm, "
                   + "Lat:%+6.2f mm, Rad:%+6.2f mm, "
                   + "Lat-Lon:%+6.2f mm, Lat-Lon-Rad:%+6.2f mm")
                  % (mdl.acro, tom, tam, trm,
                     math.sqrt(tom**2 + tam**2),
                     math.sqrt(tom**2 + tam**2 + trm**2)))
            aoe = np.abs(mdl.average_lon_error())*1000
            aae = np.abs(mdl.average_lat_error())*1000
            are = np.abs(mdl.average_rad_error())*1000
            print(("Average error for > %5s: Lon:%+6.2f mm, "
                   + "Lat:%+6.2f mm, Rad:%+6.2f mm, "
                   + "Lat-Lon:%+6.2f mm, Lat-Lon-Rad:%+6.2f mm")
                  % (mdl.acro, aoe, aae, are,
                     math.sqrt(aoe**2 + aae**2),
                     math.sqrt(aoe**2 + aae**2 + are**2)))
            print((" Error/Move       > %5s: Lon:%+6.2f %%, "
                   + "Lat:%+6.2f %%, Rad:%+6.2f %%, "
                   + "Lat-Lon:%+6.2f %%, Lat-Lon-Rad:%+6.2f %%")
                  % (mdl.acro, 100*aoe/tom, 100*aae/tam, 100*are/trm,
                     100*math.sqrt(aoe**2 + aae**2)/math.sqrt(tom**2 + tam**2),
                     100*math.sqrt(aoe**2 + aae**2 + are**2)/math.sqrt(tom**2 + tam**2 + trm**2)))

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
        self.ax1 = self.fig.add_subplot(111)
        self.plt_stn = MapPlot(self.ax1)
        self.settings = {"linestyle": "None"}
        
    def draw_map_background(self):
        self.plt_stn.plot_bg()

    def present(self):
        plt.show()

    def plot_lon_lat(self):
        for acro, mdl in self.sm_mgr.sm.items():
            coords = mdl.coords
            (xi, eta) = (coords.theta(), coords.phi())
            self.plt_stn.add_station(eta.ref, xi.ref, "%s\n(%.0f mm)" % (acro, 1000*mdl.total_rad_move()))

    def plot_laurentide_center(self, phi, theta, a, b, c1, c2):
        rad = list()
        for mdl in self.sm_mgr.sm.values():
            rad.append(mdl.coords.rad().ref)
        R = np.average(rad)
        (U, V) = np.meshgrid(np.linspace(phi - pi/2, phi + pi/2, 501),
                             np.linspace(theta - pi/2, theta + pi/2, 501))
        dphi = U - phi
        dtheta = V - theta
        d = R*(dphi + 1j*dtheta)
        f = a + b*np.exp(-(np.real(c1*d)**2 + np.imag(c2*d)**2))
        self.plt_stn.add_station(phi, theta, "$C_{L}$\n", dotcolor="blue")
        xpt, ypt = self.plt_stn.lat_lon_to_map_xy(U, V)
        self.plt_stn.plot_contour(xpt, ypt, f)
