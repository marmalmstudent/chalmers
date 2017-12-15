import numpy as np
import csv
import copy


def csvl(fpath, years):
    p = point()
    with open(fpath, "r") as f:
        for l in csv.reader(f, delimiter=","):
            if int(float(l[0])) in years:
                p.add(float(l[0]), float(l[1]), float(l[2]))
    return p


def within_tol(point1, point2, tol):
    return abs(point1.val - point2.val) < tol


class coordinate(object):
    def __init__(self, years, diff_tol=5.0):
        self.years = years
        self.diff_tol = diff_tol
        self.time = list()
        self.lat = point()
        self.lon = point()
        self.rad = point()

    def from_file(self, file_lon, file_lat, file_rad):
        self.fill_coordinates(csvl(file_lon, self.years),
                              csvl(file_lat, self.years),
                              csvl(file_rad, self.years))

    def fill_coordinates(self, lon_points, lat_points, rad_points):
        self.time = copy.copy(rad_points.time)
        self.lon = lon_points
        self.lat = lat_points
        self.rad = rad_points

    def condition(self, i, p, mean):
        cond = False
        if p.val[i] > p.val[i-1]:
            cond = cond or ((p.val[i-1]+p.error[i-1])
                            - (p.val[i]-p.error[i])) < self.diff_tol
        else:
            cond = cond or ((p.val[i-1]-p.error[i-1])
                            - (p.val[i]+p.error[i])) < self.diff_tol
        if p.val[i] > p.val[i+1]:
            cond = cond or ((p.val[i+1]+p.error[i+1])
                            - (p.val[i]-p.error[i])) < self.diff_tol
        else:
            cond = cond or ((p.val[i+1]-p.error[i+1])
                            - (p.val[i]+p.error[i])) < self.diff_tol
        return cond

    def linreg_polyfit(self, point_list):
        return np.polynomial.polynomial.polyfit(
            np.array(self.time), np.array(point_list.val), 1,
            w=1.0/np.array(point_list.error), full=False)


class point(object):
    def __init__(self):
        self.time = list()
        self.val = list()
        self.error = list()

    def add(self, time, val, error):
        self.time.append(time)
        self.val.append(val)
        self.error.append(error)

    def __str__(self):
        return "[%8.4f, %6.2f, %6.2f]" % (self.time, self.val, self.error)
