import numpy as np
from collections import OrderedDict
from numpy.polynomial.polynomial import polyfit
import csv
import copy


def csvl(fpath, years):
    p = _Point()
    with open(fpath, "r") as f:
        for l in csv.reader(f, delimiter=","):
            if int(float(l[0])) in years:
                p.add(float(l[0]), float(l[1]), float(l[2]))
    return p


def list_norm(l):
    return np.sqrt(sum(np.array(l)**2))


def zrot_mat(v):
    return np.array([[np.cos(v), -np.sin(v), 0],
                     [np.sin(v), np.cos(v), 0],
                     [0, 0, 1]])


def yrot_mat(v):
    return np.array([[np.cos(v), 0, np.sin(v)],
                     [0, 1, 0],
                     [-np.sin(v), 0, np.cos(v)]])


class Coordinate(object):
    def __init__(self, years, diff_tol=10.0):
        self.years = years
        self.diff_tol = diff_tol
        self.time = list()

        self._spherical = False
        self.xi = _Point()
        self.eta = _Point()
        self.zeta = _Point()

    def set_ref_pos(self, x, y, z):
        theta = np.arctan2(z, list_norm((x, y)))
        phi = np.arctan2(y, x)
        _trot = yrot_mat(-theta)
        _prot = zrot_mat(phi)
        R = np.sqrt(x**2 + y**2 + z**2)
        # reference point is as x=1, y=0, z=0
        _p_val = np.array([R + self.zeta.val, self.xi.val, self.eta.val])
        _p_error = np.array([self.zeta.error, self.xi.error, self.eta.error])
        _p_rot_val = np.dot(_prot, np.dot(_trot, _p_val))
        _p_rot_error = np.dot(_prot, np.dot(_trot, _p_error))
        _x_val = _p_rot_val[0]
        _x_error = _p_rot_error[0]
        _y_val = _p_rot_val[1]
        _y_error = _p_rot_error[1]
        _z_val = _p_rot_val[2]
        _z_error = _p_rot_error[2]
        self.xi.ref = x
        self.xi.val = _x_val
        self.xi.error = _x_error
        self.eta.ref = y
        self.eta.val = _y_val
        self.eta.error = _y_error
        self.zeta.ref = z
        self.zeta.val = _z_val
        self.zeta.error = _z_error

    def cart_to_sphere(self):
        if self._spherical:
            return
        self._spherical = True
        (x, y, z) = (self.xi, self.eta, self.zeta)
        _theta_ref = np.arctan2(z.ref, list_norm((x.ref, y.ref)))
        _theta_val = np.arctan2(z.val, list_norm((x.val, y.val)))
        _theta_error = np.arctan2(z.error, list_norm((x.error, y.error)))
        _phi_ref = np.arctan2(y.ref, x.ref)
        _phi_val = np.arctan2(y.val, x.val)
        _phi_error = np.arctan2(y.error, x.error)
        _rad_ref = list_norm((x.ref, y.ref, z.ref))
        _rad_val = list_norm((x.val, y.val, z.val))
        _rad_error = list_norm((x.error, y.error, z.error))
        self.xi.transform(_theta_ref, _theta_val, _theta_error)
        self.eta.transform(_phi_ref, _phi_val, _phi_error)
        self.zeta.transform(_rad_ref, _rad_val, _rad_error)

    def sphere_to_cart(self):
        if not self._spherical:
            return
        self._spherical = False
        (t, p, r) = (self.xi, self.eta, self.zeta)
        _x_ref = r.ref*np.cos(t.ref)*np.cos(p.ref)
        _x_val = r.val*np.cos(t.val)*np.cos(p.val)
        _x_error = r.error*np.cos(t.error)*np.cos(p.error)
        _y_ref = r.ref*np.cos(t.ref)*np.sin(p.ref)
        _y_val = r.val*np.cos(t.val)*np.sin(p.val)
        _y_error = r.error*np.cos(t.error)*np.sin(p.error)
        _z_ref = r.ref*np.sin(t.ref)
        _z_val = r.val*np.sin(t.val)
        _z_error = r.error*np.sin(t.error)
        self.xi.transform(_x_ref, _x_val, _x_error)
        self.eta.transform(_y_ref, _y_val, _y_error)
        self.zeta.transform(_z_ref, _z_val, _z_error)

    def linreg_polyfit(self, point_list):
        return polyfit(self.time, point_list.val, 1, w=1.0/point_list.error)

    def from_file(self, lon_file, lat_file, rad_file):
        (self.time, self.xi, self.eta, self.zeta) = self._fill_coordinates(
            csvl(lon_file, self.years),
            csvl(lat_file, self.years),
            csvl(rad_file, self.years))
        self.time = np.array(self.time, dtype=np.float64)
        self.xi.lock()
        self.xi.val *= 1e-2
        self.xi.error *= 1e-2
        self.eta.lock()
        self.eta.val *= 1e-2
        self.eta.error *= 1e-2
        self.zeta.lock()
        self.zeta.val *= 1e-2
        self.zeta.error *= 1e-2

    def _fill_coordinates(self, xi, eta, zeta):
        xi = self._form_yearly_average(self._sort_by_year(xi))
        eta = self._form_yearly_average(self._sort_by_year(eta))
        zeta = self._form_yearly_average(self._sort_by_year(zeta))

        (fxi, feta, fzeta) = self._fwd_fill(xi, eta, zeta)
        time = copy.copy(fzeta.time)

        (bxi, beta, bzeta) = self._bwd_fill(xi, eta, zeta)
        avg = np.average(fzeta.val, weights=fzeta.error)
        for (i, t) in enumerate(bzeta.time):
            bzeta.val[i] -= bzeta.val[i] - avg

        # add points from backward iteration, if not already added in forward
        for (i, t) in enumerate(bzeta.time):
            if t not in fzeta.time:
                time.append(t)
                fxi.add(beta.time[i], bxi.val[i], bxi.error[i])
                feta.add(bxi.time[i], beta.val[i], beta.error[i])
                fzeta.add(bzeta.time[i], bzeta.val[i], bzeta.error[i])

        # average each year
        return (time, fxi, feta, fzeta)

    def _sort_by_year(self, p):
        _data = dict()
        for (t, v, e) in zip(p.time, p.val, p.error):
            k = int(t)
            if k not in _data:
                _data[k] = dict()
                _data[k] = _Point()
            _data[k].add(t, v, e)
        return _data

    def _form_yearly_average(self, _data):
        _point = _Point()
        for p in OrderedDict(sorted(_data.items())).values():
            w = 1/np.array(p.error)
            _point.add(np.average(p.time),
                       np.average(p.val, weights=w),
                       np.average(p.error, weights=w))
        return _point

    def _fwd_fill(self, sxi, seta, szeta):
        dxi = _Point()
        deta = _Point()
        dzeta = _Point()
        j = 0
        dxi.add(sxi.time[j], sxi.val[j], sxi.error[j])
        deta.add(seta.time[j], seta.val[j], seta.error[j])
        dzeta.add(szeta.time[j], szeta.val[j], szeta.error[j])
        for (i, t) in enumerate(szeta.time):
            if self._condition(dzeta, j, szeta, i):
                dxi.add(sxi.time[i], sxi.val[i], sxi.error[i])
                deta.add(seta.time[i], seta.val[i], seta.error[i])
                dzeta.add(szeta.time[i], szeta.val[i], szeta.error[i])
                j += 1
        return (deta, dxi, dzeta)

    def _bwd_fill(self, sxi, seta, szeta):
        sxi.reverse()
        seta.reverse()
        szeta.reverse()
        (dxi, deta, dzeta) = self._fwd_fill(sxi, seta, szeta)
        sxi.reverse()
        seta.reverse()
        szeta.reverse()

        dxi.reverse()
        deta.reverse()
        dzeta.reverse()
        return (dxi, deta, dzeta)

    def _condition(self, p_tst, i_tst, p_cmp, i_cmp):
        return p_tst.within_tol(i_tst, p_cmp, i_cmp, self.diff_tol)


class _Point(object):
    def __init__(self):
        self.ref = 0
        self.time = list()
        self.val = list()
        self.error = list()

    def lock(self):
        self.time = np.array(self.time, dtype=np.float64)
        self.val = np.array(self.val, dtype=np.float64)
        self.error = np.array(self.error, dtype=np.float64)

    def add(self, time, val, error):
        self.time.append(time)
        self.val.append(val)
        self.error.append(error)

    def reverse(self):
        self.time.reverse()
        self.val.reverse()
        self.error.reverse()

    def within_tol(self, i, p, j, tol):
        if self.val[i] > p.val[j]:
            return (self.val[i]-self.error[i]) - (p.val[j]+p.error[j]) < tol
        else:
            return (p.val[j]-p.error[j]) - (self.val[i]+self.error[i]) < tol

    def transform(self, ref, val, error):
        self.ref = ref
        self.val = val
        self.error = error

    def __str__(self):
        out = np.array([self.time, self.val, self.error])
        return str(out.transpose().tolist())
