import numpy as np
from Coordinate import Coordinate
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d


class StationModel(object):
    def __init__(self, acronyme, path_to_csv, years):
        self.acro = acronyme
        self.coords = Coordinate(years)
        self.coords.from_file("%s/%s.lon.csv" % (path_to_csv, self.acro),
                              "%s/%s.lat.csv" % (path_to_csv, self.acro),
                              "%s/%s.rad.csv" % (path_to_csv, self.acro))

    def total_lon_move(self):
        xi = self.coords.xi
        eta = self.coords.eta
        zeta = self.coords.zeta
        if self.coords._spherical:
            bgn = zeta.val[0]*np.cos(xi.val[0])*np.cos(eta.val[0])
            end = zeta.val[-1]*np.cos(xi.val[-1])*np.cos(eta.val[-1])
            return end - bgn
        else:
            return xi.val[-1] - xi.val[0]

    def total_lat_move(self):
        xi = self.coords.xi
        eta = self.coords.eta
        zeta = self.coords.zeta
        if self.coords._spherical:
            bgn = zeta.val[0]*np.cos(xi.val[0])*np.sin(eta.val[0])
            end = zeta.val[-1]*np.cos(xi.val[-1])*np.sin(eta.val[-1])
            return end - bgn
        else:
            return eta.val[-1] - eta.val[0]

    def total_rad_move(self):
        xi = self.coords.xi
        eta = self.coords.eta
        zeta = self.coords.zeta
        if self.coords._spherical:
            return zeta.val[-1] - zeta.val[0]
        else:
            bgn = np.sqrt(xi.val[0]**2 + eta.val[0]**2 + zeta.val[0]**2)
            end = np.sqrt(xi.val[-1]**2 + eta.val[-1]**2 + zeta.val[-1]**2)
            return end - bgn


class StationViewer(object):
    class Errorbar3D(object):
        def __init__(self, ax3d, xval, xerr, yval, yerr, zval, zerr):
            self._ax = ax3d
            self._xv = xval
            self._xe = xerr
            self._yv = yval
            self._ye = yerr
            self._zv = zval
            self._ze = zerr

        def plot_xyz(self):
            self._ax.plot(self._xv, self._yv, self._zv,
                          linestyle="None", marker=".", color="#000000")
            self._plot_errorbars_3d()

        def _plot_errorbars_3d(self):
            itr = zip(self._xv, self._yv, self._zv,
                      self._xe, self._ye, self._ze)
            for (_xv, _yv, _zv, _xe, _ye, _ze) in itr:
                self._ax.plot([_xv+_xe, _xv-_xe], [_yv, _yv], [_zv, _zv],
                              marker="None", color="#ff0000")
                self._ax.plot([_xv, _xv], [_yv+_ye, _yv-_ye], [_zv, _zv],
                              marker="None", color="#00ff00")
                self._ax.plot([_xv, _xv], [_yv, _yv], [_zv+_ze, _zv-_ze],
                              marker="None", color="#0000ff")

    def __init__(self, smodel):
        self.smodel = smodel
        self._fig = plt.figure()
        self._ax1 = self._fig.add_subplot(321)
        self._ax2 = self._fig.add_subplot(323)
        self._ax3 = self._fig.add_subplot(325)
        self._ax4 = self._fig.add_subplot(122, projection='3d')
        self._ax4.set_xlabel("Longitude (West-East)")
        self._ax4.set_ylabel("Latitude (North-South)")
        self._ax4.set_zlabel("Radial (Up-Down)")
        odict = self.smodel.coords.xi
        adict = self.smodel.coords.eta
        rdict = self.smodel.coords.zeta
        self._errbar = self.Errorbar3D(self._ax4,
                                       odict.val, odict.error,
                                       adict.val, adict.error,
                                       rdict.val, rdict.error)
        self._fig.suptitle("Data for %s" % self.smodel.acro)
        self.settings = {"linestyle": "None"}

    def present(self):
        plt.show()

    def plot_lat_lon_rad(self):
        self._errbar.plot_xyz()
        self.plot_linreg_3d(self._ax4, "#000000")
        # self._ax4.set_aspect("equal", adjustable="box")

    def plot_data(self, axis, point_dict, plot_title,
                  xlabel, ylabel, line_color):
        axis.set_title(plot_title)
        axis.plot()
        axis.errorbar(self.smodel.coords.time,
                      point_dict.val, yerr=point_dict.error,
                      ecolor=line_color, **self.settings)
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)

    def plot_linreg_3d(self, axis, line_color):
        m_o, c_o = self.smodel.coords.linreg_polyfit(self.smodel.coords.xi)
        m_a, c_a = self.smodel.coords.linreg_polyfit(self.smodel.coords.eta)
        m_r, c_r = self.smodel.coords.linreg_polyfit(self.smodel.coords.zeta)
        t = self.smodel.coords.time
        o = m_o + c_o*t
        a = m_a + c_a*t
        r = m_r + c_r*t
        axis.plot(o, a, r, color=line_color, linestyle="--", linewidth=1)

    def plot_linreg(self, axis, point_list, line_color):
        m, c = self.smodel.coords.linreg_polyfit(point_list)
        t = self.smodel.coords.time
        y = m + c*np.array(t)
        axis.plot(t, y, color=line_color, linestyle="--", linewidth=1)

    def plot_lon(self):
        self.plot_data(self._ax2, self.smodel.coords.xi,
                       "", "Year", "Longitude (East-West)", "#ff0000")
        self.plot_linreg(self._ax2, self.smodel.coords.xi, "#000000")

    def plot_lat(self):
        self.plot_data(self._ax1, self.smodel.coords.eta,
                       "", "Year", "Latitude (North-South)", "#00ff00")
        self.plot_linreg(self._ax1, self.smodel.coords.eta, "#000000")

    def plot_rad(self):
        self.plot_data(self._ax3, self.smodel.coords.zeta,
                       "", "Year", "Height (Up-Down)", "#0000ff")
        self.plot_linreg(self._ax3, self.smodel.coords.zeta, "#000000")
