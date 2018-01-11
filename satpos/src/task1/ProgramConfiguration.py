from .Station import StationViewer
from os import listdir


class ProgramConfiguration(object):
    def __init__(self, respath):
        self.respath = respath
        fstr = ".rad.csv"
        self._years = [i for i in range(2005, 2018)]
        self.valid_plot_types = {"lon": StationViewer.plot_lat,
                                 "lat": StationViewer.plot_lon,
                                 "rad": StationViewer.plot_rad,
                                 "3d": StationViewer.plot_lat_lon_rad}
        self.valid_stations = [f.replace(fstr, "") for f in listdir(respath) if f.endswith(fstr)]
        self.stations = list()
        self.plot_types = list()
        self.years = list()
        
    def _year_range(self, vals):
        start = int(vals[0]) if vals[0] else min(self._years)
        if len(vals) == 2:
            step = 1
            stop = int(vals[1]) if vals[1] else max(self._years)
        else:
            step = int(vals[1]) if vals[1] else 1
            stop = int(vals[2]) if vals[2] else max(self._years)
        return range(start, stop+1, step)

    def find_settings(self, argv, look_for_str):
        s = sum([p.replace(look_for_str, "").split(",") for p in argv if p.find(look_for_str) >= 0], [])
        entries = list()
        for e in s:
            if ":" in e:
                try:
                    entries.extend(i for i in self._year_range(e.split(":")))
                except ValueError:
                    print("%s is not a range!" % e)
            else:
                entries.append(e)
        return entries

    def filter_list(self, raw_list, list_of_valids):
        return list_of_valids if "all" in raw_list else list(set([p for p in raw_list if p in list_of_valids]))

    def apply_settings(self, argv):
        self.stations.extend(s for s in self.filter_list(self.find_settings(argv, "--station="), self.valid_stations))
        self.plot_types.extend(t for t in self.filter_list(self.find_settings(argv, "--plot="), self.valid_plot_types.keys()))
        self.years.extend(y for y in self.filter_list(self.find_settings(argv, "--year="), self._years))
        self.yearly_average = True if self.find_settings(argv, "--yearly-average") else False
        self.correct_errors = True if self.find_settings(argv, "--correct-errors") else False
        self.find_center = True if self.find_settings(argv, "--find-center") else False
        self.plot_map = True if self.find_settings(argv, "--plot-map") else False
