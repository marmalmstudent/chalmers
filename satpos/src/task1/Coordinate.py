from numpy import array, cos, sin, arctan2, sqrt, dot, average, float64
from collections import OrderedDict
from numpy.polynomial.polynomial import polyfit
from csv import reader
from copy import copy


def csvl(fpath, years):
    p = _Point()
    with open(fpath, "r") as f:
        for l in reader(f, delimiter=","):
            if int(float(l[0])) in years:
                p.add(float(l[0]), float(l[1]), float(l[2]))
    return p


def list_norm(l): return sqrt(sum(array(l)**2))


def zrot_mat(v): return array([[ cos(v), -sin(v),      0],
                               [ sin(v),  cos(v),      0],
                               [      0,       0,      1]])


def yrot_mat(v): return array([[  cos(v),      0, sin(v)],
                               [       0,      1,      0],
                               [ -sin(v),      0, cos(v)]])


def cart_to_theta(x, y, z): return arctan2(z, list_norm((x, y)))


def cart_to_phi(x, y, z): return arctan2(y, x)


def cart_to_rad(x, y, z): return list_norm((x, y, z))


def sphere_to_x(theta, phi, r): return r*cos(theta)*cos(phi)


def sphere_to_y(theta, phi, r): return r*cos(theta)*sin(phi)


def sphere_to_z(theta, phi, r): return r*sin(theta)


class Coordinate(object):
    def __init__(self, years, yearly_average=False, correct_errors=False, diff_tol=10.0):
        self.years = years
        self.diff_tol = diff_tol
        self.yearly_average = yearly_average
        self.correct_errors = correct_errors
        self._t = list()

        self._spherical = False
        self._xi = _Point()
        self._eta = _Point()
        self._zeta = _Point()

    def time(self): return self._t

    def theta(self): return self._eta

    def phi(self): return self._xi

    def rad(self): return self._zeta

    def x(self): return self._xi

    def y(self): return self._eta

    def z(self): return self._zeta

    def set_ref_pos(self, x, y, z):
        y_rot_matrix = yrot_mat(-arctan2(z, list_norm((x, y))))
        z_rot_matrix = zrot_mat(arctan2(y, x))
        # x-axis is radial, y-axis is longitude, z-axis is latitude
        R = list_norm((x, y, z))
        p_val = array([R + self._zeta.val, self._xi.val, self._eta.val])
        p_error = array([self._zeta.error, self._xi.error, self._eta.error])
        p_rot_val = dot(z_rot_matrix, dot(y_rot_matrix, p_val))
        p_rot_error = dot(z_rot_matrix, dot(y_rot_matrix, p_error))
        self._xi.transform(x, p_rot_val[0], p_rot_error[0])
        self._eta.transform(y, p_rot_val[1], p_rot_error[1])
        self._zeta.transform(z, p_rot_val[2], p_rot_error[2])

    def cart_to_sphere(self):
        if self._spherical:
            return
        self._spherical = True
        (x, y, z) = (self._xi, self._eta, self._zeta)
        phi_ref = cart_to_phi(x.ref, y.ref, z.ref)
        phi_val = cart_to_phi(x.val, y.val, z.val)
        phi_error = cart_to_phi(x.error, y.error, z.error)
        theta_ref = cart_to_theta(x.ref, y.ref, z.ref)
        theta_val = cart_to_theta(x.val, y.val, z.val)
        theta_error = cart_to_theta(x.error, y.error, z.error)
        rad_ref = cart_to_rad(x.ref, y.ref, z.ref)
        rad_val = cart_to_rad(x.val, y.val, z.val)
        rad_error = cart_to_rad(x.error, y.error, z.error)
        self._xi.transform(phi_ref, phi_val, phi_error)
        self._eta.transform(theta_ref, theta_val, theta_error)
        self._zeta.transform(rad_ref, rad_val, rad_error)

    def linreg_polyfit(self, point_list):
        return polyfit(self._t, point_list.val, 1, w=1.0/point_list.error)

    def from_file(self, lon_file, lat_file, rad_file):
        (self._t, self._xi, self._eta, self._zeta) = self._fill_coordinates(
            csvl(lon_file, self.years),
            csvl(lat_file, self.years),
            csvl(rad_file, self.years))
        self._t = array(self._t, dtype=float64)
        self._xi.lock()
        self._eta.lock()
        self._zeta.lock()
        # values are in cm; transform to m
        self._xi.val *= 1e-2
        self._xi.error *= 1e-2
        self._eta.val *= 1e-2
        self._eta.error *= 1e-2
        self._zeta.val *= 1e-2
        self._zeta.error *= 1e-2

    def _fill_coordinates(self, xi, eta, zeta):
        if self.yearly_average:
            xi = self._form_yearly_average(self._sort_by_year(xi))
            eta = self._form_yearly_average(self._sort_by_year(eta))
            zeta = self._form_yearly_average(self._sort_by_year(zeta))

        if self.correct_errors:
            (fxi, feta, fzeta) = self._fwd_fill(xi, eta, zeta)
            time = copy(fzeta.time)

            (bxi, beta, bzeta) = self._bwd_fill(xi, eta, zeta)
            m_r, c_r = polyfit(fzeta.time, fzeta.val, 1, w=1.0/array(fzeta.error))
            bzeta.val -= bzeta.val[0] - (m_r + c_r*bzeta.time[0])

            # add points from backward iteration, if not already added in forward
            for (i, t) in enumerate(bzeta.time):
                if t not in fzeta.time:
                    time.append(t)
                    fxi.add(beta.time[i], bxi.val[i], bxi.error[i])
                    feta.add(bxi.time[i], beta.val[i], beta.error[i])
                    fzeta.add(bzeta.time[i], bzeta.val[i], bzeta.error[i])
            return (time, fxi, feta, fzeta)
        return (copy(zeta.time), xi, eta, zeta)

    def _sort_by_year(self, p):
        yearly_points = dict()
        for (t, v, e) in zip(p.time, p.val, p.error):
            k = int(t)
            if k not in yearly_points:
                yearly_points[k] = _Point()
            yearly_points[k].add(t, v, e)
        return yearly_points

    def _form_yearly_average(self, yearly_points):
        _point = _Point()
        for p in OrderedDict(sorted(yearly_points.items())).values():
            w = 1/array(p.error)
            _point.add(average(p.time),
                       average(p.val, weights=w),
                       average(p.error, weights=w))
        return _point

    def _fwd_fill(self, sxi, seta, szeta):
        j = 0
        dxi = _Point()
        dxi.add(sxi.time[j], sxi.val[j], sxi.error[j])
        deta = _Point()
        deta.add(seta.time[j], seta.val[j], seta.error[j])
        dzeta = _Point()
        dzeta.add(szeta.time[j], szeta.val[j], szeta.error[j])
        for (i, _) in enumerate(szeta.time):
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
        self.time = array(self.time, dtype=float64)
        self.val = array(self.val, dtype=float64)
        self.error = array(self.error, dtype=float64)

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
        out = array([self.time, self.val, self.error])
        return str(out.transpose().tolist())
