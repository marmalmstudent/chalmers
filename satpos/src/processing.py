import station
import os
import math


class prog_settings(object):
    def __init__(self):
        self.respath = "res/csv"
        self.fstr = ".rad.csv"
        self.valid_years = [i for i in range(2005, 2018)]
        self.valid_plot_types = {"lon": station.station_viewer.plot_lat,
                                 "lat": station.station_viewer.plot_lon,
                                 "rad": station.station_viewer.plot_rad,
                                 "3d": station.station_viewer.plot_lat_lon_rad}
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
                    start = int(vals[0]) if vals[0] else min(self.valid_years)
                    if len(vals) == 2:
                        step = 1
                        stop = int(vals[1]) if vals[1] else max(self.valid_years)
                    else:
                        step = int(vals[1]) if vals[1] else 1
                        stop = int(vals[2]) if vals[2] else max(self.valid_years)
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
                                  self.valid_years):
            self.years.append(y)


if __name__ == "__main__":
    try:
        settings = prog_settings()
        settings.apply_settings(os.sys.argv[1:])
        if len(settings.stations) == 0:
            raise Exception("Unknown station")
        sm = [station.station_model(stn, settings.respath, settings.years)
              for stn in settings.stations]
        for stn_mdl in sm:
            print(("Total move for > %5s: Lon:%+6.2f cm, "
                   + "Lat:%+6.2f cm, Rad:%+6.2f cm, "
                   + "Lat-Lon:%+6.2f cm, Lat-Lon-Rad:%+6.2f cm")
                  % (stn_mdl.acro, stn_mdl.total_lon_move(),
                     stn_mdl.total_lat_move(), stn_mdl.total_rad_move(),
                     math.sqrt(stn_mdl.total_lon_move()**2
                               + stn_mdl.total_lat_move()**2),
                     math.sqrt(stn_mdl.total_lon_move()**2
                               + stn_mdl.total_lat_move()**2
                               + stn_mdl.total_rad_move()**2)))
        if (len(settings.plot_types) > 0):
            if (len(sm) > 1):
                print("Plotting data for first station (%s)"
                      % settings.stations[0])
            sv = station.station_viewer(sm[0])
            for key in settings.plot_types:
                settings.valid_plot_types[key](sv)
            sv.present()
    except ValueError:
        print("ValueError")
    """
    except IndexError:
        print("Usage: python %s <station-acronyme>" % os.sys.argv[0])
    except FileNotFoundError:
        print(("Could not find data for station with acronyme '%s'\n"
              + "Try one of the following: '%s'")
              % (os.sys.argv[1], "', '".join(settings.valid_stations)))
    """
