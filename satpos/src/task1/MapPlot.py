from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np


def rtod(r):
    return r*180/np.pi


class MapPlot(object):
    def __init__(self, ax):
        self.ax = ax
        # setup Lambert Conformal basemap.
        self.m = Basemap(width=12e6, height=9e6,
                         projection='lcc', resolution=None,
                         lat_1=45., lat_2=55, lat_0=50, lon_0=-107.)

    def plot_bg(self):
        self.m.shadedrelief()

    def add_station(self, lon, lat, text, mark='o', dotcolor='red'):
        # convert to map projection coords. 
        # Note that lon,lat can be scalars, lists or numpy arrays.
        xpt,ypt = self.m(rtod(lon), rtod(lat))

        # convert back to lat/lon
        lonpt, latpt = self.m(xpt, ypt, inverse=True)
        self.m.plot(xpt, ypt, color=dotcolor, marker=mark)  # plot a blue dot there

        # put some text next to the dot, offset a little bit
        # (the offset is in map projection coordinates)
        #self.ax.text(xpt+100000, ypt+100000, text, horizontalalignment='center')
        self.ax.text(xpt, ypt-250000+100000, text, horizontalalignment='center')

    def lat_lon_to_map_xy(self, lon, lat):
        return self.m(rtod(lon), rtod(lat))

    def plot_contour(self, U, V, f):
        cs = self.ax.contour(U, V, f, cmap=plt.get_cmap("jet"))
        try:
            self.ax.clabel(cs, inline=True, colors="black")
        except UnboundLocalError:
            print("UnboundLocalError encountered when creating contour labels")


if __name__ == '__main__':
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt_stn = MapPlot(ax)
    plt_stn.plot_bg()
    # plot blue dot on Boulder, colorado and label it as such.
    lon, lat = -104.237, 40.125  # Location of Boulder
    plt_stn.add_station(lon, lat, "test: (%.1f cm)" % (14.9))
    plt.show()
