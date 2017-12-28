import traceback
from Station import StationViewer
from StationManager import StationManager, StationManagerView
from LaurentideCenterFinder import LaurentideCenterFinder
import os
import numpy as np


class ProgramConfiguration(object):
    def __init__(self):
        self.respath = "res/csv"
        self.fstr = ".rad.csv"
        self._years = [i for i in range(2005, 2018)]
        self.valid_plot_types = {"lon": StationViewer.plot_lat,
                                 "lat": StationViewer.plot_lon,
                                 "rad": StationViewer.plot_rad,
                                 "3d": StationViewer.plot_lat_lon_rad}
        self.valid_stations = [f.replace(self.fstr, "")
                               for f in os.listdir(self.respath)
                               if f.endswith(self.fstr)]
        self.stations = list()
        self.plot_types = list()
        self.years = list()

    def find_settings(self, argv, look_for_str):
        s = sum([p.replace(look_for_str, "").split(",")
                 for p in argv if p.find(look_for_str) >= 0], [])
        entries = list()
        for e in s:
            if ":" in e:
                try:
                    vals = e.split(":")
                    start = int(vals[0]) if vals[0] else min(self._years)
                    if len(vals) == 2:
                        step = 1
                        stop = int(vals[1]) if vals[1] else max(self._years)
                    else:
                        step = int(vals[1]) if vals[1] else 1
                        stop = int(vals[2]) if vals[2] else max(self._years)
                    for i in range(start, stop+1, step):
                        entries.append(i)
                except ValueError:
                    entries.append(e)
            else:
                entries.append(e)
        return entries

    def filter_list(self, raw_list, list_of_valids):
        if "all" in raw_list:
            return list_of_valids
        vals = list()
        for p in raw_list:
            if p in list_of_valids and p not in vals:
                vals.append(p)
        return vals

    def apply_settings(self, argv):
        for s in self.filter_list(self.find_settings(argv, "--station="),
                                  self.valid_stations):
            self.stations.append(s)
        for t in self.filter_list(self.find_settings(argv, "--plot="),
                                  self.valid_plot_types.keys()):
            self.plot_types.append(t)
        for y in self.filter_list(self.find_settings(argv, "--year="),
                                  self._years):
            self.years.append(y)


def plot_station(sm, cfg):
    if (len(cfg.plot_types) > 0):
        if (len(sm) > 1):
            print("Plotting data for first station (%s)"
                  % cfg.stations[0])
        sv = StationViewer(sm[cfg.stations[0]])
        for key in cfg.plot_types:
            cfg.valid_plot_types[key](sv)
        sv.present()


if __name__ == "__main__":
    try:
        cfg = ProgramConfiguration()
        cfg.apply_settings(os.sys.argv[1:])
        if len(cfg.stations) == 0:
            raise Exception("Unknown station")
        smgr = StationManager(cfg.respath, cfg.stations, cfg.years)
        plot_station(smgr.sm, cfg)
        smgr.load_station_ref_coords("res", "coords.txt")
        smgr.print_station_data()
        smgr.make_spherical()

        lcf = LaurentideCenterFinder(smgr.sm.values())
        ctheta, cphi, a, b, c1_r, c1_v, c2_r, c2_v = lcf.solver()
        c1 = c1_r*np.exp(1j*c1_v)
        c2 = c2_r*np.exp(1j*c2_v)

        sm_mgr_view = StationManagerView(smgr)
        sm_mgr_view.plot_laurentide_center(ctheta, cphi, a, b, c1, c2)
        sm_mgr_view.plot_lon_lat()
        sm_mgr_view.present()
    except ValueError:
        print("ValueError")
        traceback.print_exc()
