import numpy as np
from coordinate import coordinate
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d


class station_model(object):
    def __init__(self, acronyme, path_to_csv, years):
        self.acro = acronyme
        self.coords = coordinate(years)
        self.coords.from_file("%s/%s.lon.csv" % (path_to_csv, self.acro),
                              "%s/%s.lat.csv" % (path_to_csv, self.acro),
                              "%s/%s.rad.csv" % (path_to_csv, self.acro))

    def total_lon_move(self):
        return self.total_move(self.coords.lon.val,
                               self.coords.lon.error)

    def total_lat_move(self):
        return self.total_move(self.coords.lat.val,
                               self.coords.lat.error)

    def total_rad_move(self):
        return self.total_move(self.coords.rad.val,
                               self.coords.rad.error)

    def total_move(self, axis_val, axis_err):
        return max(axis_val) - min(axis_val)

    def time_data(self):
        return self.coords.time

    def __str__(self):
        return str(self.coords)


class station_viewer(object):
    class errorbar3d(object):
        def __init__(self, ax3d, xval, xerr, yval, yerr, zval, zerr):
            self.ax = ax3d
            self.xv = xval
            self.xe = xerr
            self.yv = yval
            self.ye = yerr
            self.zv = zval
            self.ze = zerr

        def plot_xyz(self):
            self.ax.plot(self.xv, self.yv, self.zv,
                         linestyle="None", marker=".", color="#000000")
            self.plot_errorbars_3d()

        def plot_errorbars_3d(self):
            for (xv, yv, zv, xe, ye, ze) in zip(
                    self.xv, self.yv, self.zv, self.xe, self.ye, self.ze):
                self.ax.plot([xv+xe, xv-xe], [yv, yv], [zv, zv],
                             marker="_", color="#ff0000")
                self.ax.plot([xv, xv], [yv+ye, yv-ye], [zv, zv],
                             marker="_", color="#00ff00")
                self.ax.plot([xv, xv], [yv, yv], [zv+ze, zv-ze],
                             marker="_", color="#0000ff")

    def __init__(self, smodel):
        self.smodel = smodel
        self.fig = plt.figure()
        self.ax1 = self.fig.add_subplot(321)
        self.ax2 = self.fig.add_subplot(323)
        self.ax3 = self.fig.add_subplot(325)
        self.ax4 = self.fig.add_subplot(122, projection='3d')
        self.ax4.set_xlabel("Longitude (West-East)")
        self.ax4.set_ylabel("Latitude (North-South)")
        self.ax4.set_zlabel("Radial (Up-Down)")
        odict = self.smodel.coords.lon
        adict = self.smodel.coords.lat
        rdict = self.smodel.coords.rad
        self.errbar = self.errorbar3d(self.ax4,
                                      odict.val, odict.error,
                                      adict.val, adict.error,
                                      rdict.val, rdict.error)
        self.fig.suptitle("Data for %s" % self.smodel.acro)
        self.settings = {"linestyle": "None"}

    def present(self):
        plt.show()

    def plot_lat_lon_rad(self):
        self.errbar.plot_xyz()
        self.plot_linreg_3d(self.ax4, "#000000")
        # self.ax4.set_aspect("equal", adjustable="box")

    def plot_data(self, axis, point_dict, plot_title,
                  xlabel, ylabel, line_color):
        axis.set_title(plot_title)
        axis.plot()
        axis.errorbar(self.smodel.time_data(),
                      point_dict.val,
                      yerr=point_dict.error,
                      ecolor=line_color,
                      **self.settings)
        axis.set_xlabel(xlabel)
        axis.set_ylabel(ylabel)

    def plot_linreg_3d(self, axis, line_color):
        m_o, c_o = self.smodel.coords.linreg_polyfit(self.smodel.coords.lon)
        m_a, c_a = self.smodel.coords.linreg_polyfit(self.smodel.coords.lat)
        m_r, c_r = self.smodel.coords.linreg_polyfit(self.smodel.coords.rad)
        t = np.array(self.smodel.coords.time)
        o = m_o + c_o*t
        a = m_a + c_a*t
        r = m_r + c_r*t
        axis.plot(o, a, r, color=line_color, linestyle="--", linewidth=1)

    def plot_linreg(self, axis, point_list, line_color):
        m, c = self.smodel.coords.linreg_polyfit(point_list)
        x = np.array(self.smodel.coords.time)
        y = m + c*x
        axis.plot(x, y, color=line_color, linestyle="--", linewidth=1)

    def plot_lon(self):
        self.plot_data(self.ax2, self.smodel.coords.lon,
                       "", "Year", "Longitude (East-West)", "#ff0000")
        self.plot_linreg(self.ax2, self.smodel.coords.lon, "#000000")

    def plot_lat(self):
        self.plot_data(self.ax1, self.smodel.coords.lat,
                       "", "Year", "Latitude (North-South)", "#00ff00")
        self.plot_linreg(self.ax1, self.smodel.coords.lat, "#000000")

    def plot_rad(self):
        self.plot_data(self.ax3, self.smodel.coords.rad,
                       "", "Year", "Height (Up-Down)", "#0000ff")
        self.plot_linreg(self.ax3, self.smodel.coords.rad, "#000000")
