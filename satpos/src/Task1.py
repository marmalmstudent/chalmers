import traceback
from task1.Station import StationViewer
from task1.StationManager import StationManager, StationManagerView
from task1.LaurentideCenterFinder import LaurentideCenterFinder
from task1.ProgramConfiguration import ProgramConfiguration
from os import sys
from numpy import pi


def plot_station(sm, cfg):
    if (len(cfg.plot_types) > 0):
        if (len(sm) > 1):
            print("Plotting data for first station (%s)" % cfg.stations[0])
        sv = StationViewer(sm[cfg.stations[0]])
        for key in cfg.plot_types:
            cfg.valid_plot_types[key](sv)
        sv.present()


if __name__ == "__main__":
    try:
        respath = sys.argv[1] + "/task1"
        cfg = ProgramConfiguration(respath + "/csv")
        cfg.apply_settings(sys.argv[2:])
        if len(cfg.stations) == 0:
            raise Exception("Unknown station")
        
        smgr = StationManager(cfg.respath, cfg.stations, cfg.years, cfg.yearly_average, cfg.correct_errors)
        plot_station(smgr.sm, cfg)
        smgr.load_station_ref_coords(respath, "coords.txt")
        smgr.make_spherical()
        smgr.print_station_data()
        
        if cfg.plot_map:
            sm_mgr_view = StationManagerView(smgr)
            sm_mgr_view.draw_map_background()
            if cfg.find_center:
                lcf = LaurentideCenterFinder(smgr.sm.values())
                cphi, ctheta, a, b, c1, c2 = lcf.solver()
                print(ctheta*180/pi, cphi*180/pi, a, b, c1, c2)
                sm_mgr_view.plot_laurentide_center(cphi, ctheta, a, b, c1, c2)
            sm_mgr_view.plot_lon_lat()
            sm_mgr_view.present()
    except ValueError:
        print("ValueError")
        traceback.print_exc()
