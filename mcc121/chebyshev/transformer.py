import numpy as np
import numpy.polynomial.chebyshev as cheb
import numpy.matlib as npm
import scipy.optimize as spo
import matplotlib.pyplot as plt
pi = np.pi
c_0 = 2.99792458e8
mu_0 = pi*4e-7
epsilon_0 = 1/(mu_0*np.power(c_0, 2))
z_0 = mu_0*c_0


class transformer(object):
    """
    transformer(object)

    This class is a transformer that supports variable dielectric or variable
    dimentsions.
    """
    gamma_m = 0.07  # Group 6
    load_ratio = 3.4  # Group 6
    bw_frac = 1  # Group 6

    def __init__(self, a, b, f_0):
        """
        __init__(self, a, b, f_0)

        Class initialzation.

        Parameters
        ----------
        a : number
            length along x-axis (a>b).
        b : number
            length along y-axis (b<a).
        f_0 : number
            Operating frequency of the transformer.
        """
        self.a = a  # waveguide size along x-axis
        self.b = b  # waveguide size along y-axis
        self.f_0 = f_0  # frequency
        self.e_r_wg = 1  # relative permittivity in waveguide
        self.mu_r_wg = 1  # relative permeability in waveguide
        self.e_r = None  # relative permittivity in transformer sections
        self.mu_r = 1  # relative permeability in transformer sections
        self.m = None  # modes along x-axis that can propagate in waveguide
        self.n = None  # modes along y-axis that can propagate in waveguide
        self.bool_mn = None  # comb of m and n that can propagate in waveguide
        self.nbr_sctn = None  # number of sections in them transformer
        self.sctn_len = None  # length of sections in transformer
        # normalized impedances of sections in transformer
        self.z_norm = np.array([1, 0, 0, 0, 0], dtype=np.float64)
        self.targ_theta_m = pi/4*(2-transformer.bw_frac)  # target bw
        self.targ_gamma = None  # target reflection coefficients
        self.gamma = None  # final reflection coefficients
        self.theta_m = None  # final bw
        self.bw_frac = None  # final fractional bandwidth
        self.g0_tmp = None  # tmp gamma_0 used in calculations

    def find_modes(self):
        """
        find_modes(self)

        Calculates which modes that can exist in a given rectangular waveguide.
        It does so by calculating the maximum m-mode when n is zero and n-mode
        when m is zero. Then it calculates looks at every combination of nodes
        that could exist. The result is a boolean matrix containing information
        about whether a particular mode propagates or not. For example, if
        bool_mn[m][n] returns true then it means that the mn-mode can
        propagate in this waveguide at this frequency.
        """
        m = np.arange(np.floor(
            2*np.sqrt(self.e_r_wg*self.mu_r_wg)*self.f_0/c_0*self.a)+1)
        n = np.arange(np.floor(
            2*np.sqrt(self.e_r_wg*self.mu_r_wg)*self.f_0/c_0*self.b)+1)
        self.m = npm.repmat(m, np.size(n), 1).T
        self.n = npm.repmat(n, np.size(m), 1)
        self.bool_mn = (
            c_0/(2*np.sqrt(self.e_r_wg*self.mu_r_wg))
            * np.sqrt((self.m/self.a)**2 + (self.n/self.b)**2)) <= self.f_0
        return

    def calc_norm_imp(self):
        """
        calc_norm_imp(self)

        Calculates the normalized impedances, i.e. relative to the impedance of
        the waveguide. In this version of the program it only works with three
        segments.

        Returns
        -------
        z_norm: ndarray
            the normalized impedances for the segments of the transformer.
        """
        gamma_0 = -transformer.gamma_m/(2*np.cos(self.targ_theta_m)**3)
        gamma_1 = -3/2*transformer.gamma_m*(
            1/np.cos(self.targ_theta_m)**3 -
            1/np.cos(self.targ_theta_m))
        self.targ_gamma = np.array([gamma_0, gamma_1, gamma_1, gamma_0])
        self.correct_props()  # correct provided values to solve the problem
        for i in range(0, np.size(self.z_norm)-1):
            self.z_norm[i+1] = self.z_norm[i]*((1+self.gamma[i]) /
                                               (1-self.gamma[i]))
        return

    def correct_props(self):
        """
        correct_props(self)

        Calculates the bandwidth of them transformer, the reflection
        coefficients at the interfaces between the sections in the transformer,
        and calculates gamma_m for these values

        Returns
        -------
        """
        # inital value of g0_tmp
        self.g0_tmp = -transformer.gamma_m/(2*np.cos(self.targ_theta_m)**3)
        self.theta_m = spo.fmin(
            lambda targ_theta_m: self.main_cond(targ_theta_m, self.g0_tmp,
                                                3*np.sin(targ_theta_m)**2),
            self.targ_theta_m)[0]
        self.bw_frac = (2-4*self.theta_m/pi)
        self.gamma = np.array(
            [self.g0_tmp, 3*np.sin(self.theta_m)**2*self.g0_tmp,
             3*np.sin(self.theta_m)**2*self.g0_tmp, self.g0_tmp])
        self.gamma_m = np.abs(self.gamma[0]*(2*np.cos(self.theta_m)**3))
        return

    def main_cond(self, targ_theta_m, g, ratio):
        """
        main_cond(self, targ_theta_m, g, ratio)

        This is them main condition that will be used to calculated them values
        of targ_theta_m, gamma_0 and gamma_1. It will calculate these values
        so that the overall norm of the deviation from the inital values is as
        small as possible.

        Parameters
        ----------
        targ_theta_m : number
            The bandwidth expressed in electrical length.
        g : number
            The reflection coefficient for the first interface, gamma_0
        ratio : number
            The ratio between gamma_1 and gamma_0

        Returns
        -------
        new_theta_m : number
            The new bandwidth that minimizes the overall norm of the
            deviation
        new_g : number
            The new reflection coefficient at them first interface that
            minimizes the overall norm of the deviation.
        """
        print(np.linalg.norm(np.array([
            np.abs(targ_theta_m - self.targ_theta_m) / self.targ_theta_m,
            np.abs(self.g0_tmp - self.targ_gamma[0]) / self.targ_gamma[0],
            np.abs(3*np.sin(targ_theta_m)**2*self.g0_tmp-self.targ_gamma[1]) /
            self.targ_gamma[1]])))
        # optimization of g
        self.g0_tmp = spo.fsolve(
            lambda g: self.gamma_cond(g, 3*np.sin(targ_theta_m)**2), g)[0]
        # calculate residual (optimization of targ_theta_m)
        return np.linalg.norm(np.array([
            np.abs(targ_theta_m - self.targ_theta_m) / self.targ_theta_m,
            np.abs(self.g0_tmp - self.targ_gamma[0]) / self.targ_gamma[0],
            np.abs(3*np.sin(targ_theta_m)**2*self.g0_tmp-self.targ_gamma[1]) /
            self.targ_gamma[1]]))

    def gamma_cond(self, g, ratio):
        """
        gamma_cond(self, g, ratio)

        finds the difference between z_0/z_L calculated using the reflection
        coefficients at each interfaces and load_ratio that was given in the
        problem description.

        Parameters
        ----------
        g : number
            The value of gamma_0 where the function will be evaluated at.
        ratio : number
            The ratio between gamma_1 and gamma_0 that satisfies the chebyshev
            polynomial

        Returns
        -------
        residual : number
            The difference between z_0/z_L calculated using the reflection
            coefficients at each interfaces and load_ratio that was given in
            the problem description.
        """
        return ((1+g)/(1-g))**2*((1+ratio*g)/(1-ratio*g))**2 - \
            1/transformer.load_ratio

    def calc_n_sctn(self):
        """
        calc_n_sctn(self)

        Calculates the number of segments required to meet the requirements of
        the transformer. Note that this method uses an approximation that is
        only somewhat accurate when z_L is close to z_0.

        Returns
        -------
        N : int
            The number of segment needed in order to achieve a targ_theta_m
            that is as close to the supplied desired value as possible along
            with the acquired targ_theta_m for the calculated number of
            segments.
        """
        """
        This is an approximation of N because the total reflection coefficient
        is not (z_L-z_0)/(z_L+z_0) when frequency is zero.
        """
        self.nbr_sctn = int(np.round(spo.fsolve(lambda N: (
            np.cosh(np.arccosh(
                np.abs((1-transformer.load_ratio) /
                       (1+transformer.load_ratio))/transformer.gamma_m)/N)
            - 1/np.cos(self.targ_theta_m)), 4)))
        return

    def calc_charr_imp(self, m, n, e_r, mu_r):
        """
        calc_charr_imp(self, m, n, e_r, mu_r)

        Calculates the characteristic impedance of a rectangular wavegauide
        with crossection x=a and y=b.

        Parameters
        ----------
        m : int, matrix_like
            m-mode value. Dimensions must be the same as n.
        n : int, matrix_like
            n-mode value. Dimensions must be the same as m.
        e_r : number
            relative permittivity of waveguide/section
        mu_r : number
            relative permeability of waveguide/section

        Returns
        -------
        z_tm : number, matrix_like
            (characteristic) wave impedance of the waveguide.
        """
        w = 2*pi*self.f_0
        return w*mu_0/np.sqrt(w**2*mu_0*epsilon_0*e_r-(pi/self.a)**2)

    def calc_e_r(self, m, n, z_c):
        """
        calc_e_r(self, m, n, z_c)

        This methos assumes that mu_r is zeros, since it really calculates the
        product of e_r and mu_r.

        Parameters
        ----------
        m : int, matrix_like
            m-mode value. Dimensions must be the same as n.
        n : int, matrix-like
            n-mode value. Dimensions must be the same as m.
        z_c : ndarray
            (characteristic) wave impedances
        """
        w = 2*pi*self.f_0
        self.e_r = (c_0/w)**2*(w**2*mu_0**2/z_c**2+(pi/self.a)**2)
        return

    def calc_trsf_dim(self, z_c):
        """
        calc_trsf_dim(self, z_c)

        calculates dimensions for segments with variable size.

        Parameters
        ----------
        z_c : ndarray
            (characteristic) wave impedances.
        """
        self.a = self.a*np.ones(np.size(z_c))
        self.b = self.b*z_c/z_c[0]
        return

    def calc_trsf_prop(self, trsf_type=0):
        """
        calc_trsf_prop(self, trsf_type)

        calculates the transformer properties for a transformer with fixed
        dimensions and variable dielectric constant

        Parameters
        ----------
        trsf_type : int
            Then type of transformer. 0 for variable dielectric, 1 for variable
            dimensions.
        """
        self.calc_n_sctn()  # number of sections
        self.calc_norm_imp()  # normalized impedance
        self.find_modes()  # mode info
        # characteristic impedance for all posible modes
        z_c = self.calc_charr_imp(
            self.m*self.bool_mn, self.n*self.bool_mn,
            self.e_r_wg, self.mu_r_wg)*self.bool_mn
        z_10 = z_c[1][0]*self.z_norm  # characteristic impedance for 10 mode
        if (trsf_type == 0):
            # variable dielectric type
            self.calc_e_r(1, 0, z_10)
        elif (trsf_type == 1):
            # variable dimensions type
            self.calc_trsf_dim(z_10)
            self.e_r = 1
        # phase constant
        beta = np.sqrt((2*pi*self.f_0)**2*self.e_r/c_0**2 -
                       ((pi*self.m[1][0]/self.a)**2 +
                        (pi*self.n[0][0]/self.b)**2))
        lbda_g = 2*pi/beta  # guided wavelength
        self.sctn_len = lbda_g/4
        return


trsf_1 = transformer(80e-3, 10e-3, 6e9)
trsf_1.calc_trsf_prop(trsf_type=0)
trsf_2 = transformer(80e-3, 10e-3, 6e9)
trsf_2.calc_trsf_prop(trsf_type=1)
print("---------------Filter properties------------------")
print("Design ID:\t\t\t6")
print("Z_0/Z_L:\t\t\t"+str(transformer.load_ratio))
print("Target Gamma_m:\t\t\t"+str(transformer.gamma_m))
print("Target fractional bandwidth:\t" + str(transformer.bw_frac))
print("Target theta_m:\t\t\t"+str((pi/4*(2-transformer.bw_frac))) +
      "\n`-Number of segments:\t\t" + str(trsf_1.nbr_sctn))
print("Final Gamma_m:\t\t\t"+str((trsf_1.gamma_m)))
print("Final fractional bandwidth:\t" + str(np.round(trsf_1.bw_frac, 3)))
print("Final theta_m:\t\t\t"+str((trsf_1.theta_m)))
print("---------------Variable dielectric----------------")
print("Epsilon_r relative to first segment: \n\t", (trsf_1.e_r))
print("Segment length (lambda/4) [mm]: \n\t", (trsf_1.sctn_len*1e3))
print("---------------Variable dimensions----------------")
print("Waveguide height [mm]: \n\t", (trsf_2.b*1e3))
print("Segment length (lambda/4) [mm]: \n\t", (trsf_2.sctn_len*1e3))
"""
|----How-to-simulate-in-HFSS----|
Create variables: enter length as a variable (e.g. XSize = x1). An option will
popup that lets you save the variable.
Edit variables: HFSS -> Design properties
Join together: Assign dielectrics to each segment and put them together. Do not
 join just pot them together and it will work.
"""

# thet = np.linspace(trsf_1.theta_m, pi - trsf_1.theta_m, 1000)
thet = np.linspace(0, pi, 1000)
gmma1 = 2*np.exp(1j*trsf_1.nbr_sctn*thet) * \
      (trsf_1.gamma[0]*np.cos(trsf_1.nbr_sctn*thet) +
       trsf_1.gamma[1]*np.cos((trsf_1.nbr_sctn-2)*thet))
f_0 = 6e9
bw = (2-4*trsf_1.theta_m/pi)*f_0
f_span = 2*f_0
f = np.linspace(f_0-f_span/2, f_0+f_span/2, 1000)
gmma_m = trsf_1.gamma_m*np.ones(1000)
y = np.linspace(0, trsf_1.gamma_m)
f1 = (f_0-bw/2)*np.ones(np.size(y))
f2 = (f_0+bw/2)*np.ones(np.size(y))
plt.plot(f, np.abs(gmma1), linestyle="-", color="#000000")
plt.plot(f, np.abs(gmma_m), linestyle="--", color="#000000")
plt.plot(f1, y, linestyle="--", color="#000000")
plt.plot(f2, y, linestyle="--", color="#000000")
plt.xlabel("Frequency, [Hz]")
plt.ylabel("$|\Gamma|$")
print(2-4*trsf_1.targ_theta_m/pi)
plt.show()
