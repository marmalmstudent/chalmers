import numpy as np
import scipy.optimize as spo


__stability_codes = {0: "Both input and output is unstable.",
                     1: "Input is unstable, output is conditionally stable.",
                     2: "Input is conditionally stable, output is unstable.",
                     3: "Both input and output is conditionally stable.",
                     4: "Input is unstable, output is stable.",
                     6: "Input is conditionally stable, output is stable.",
                     8: "Input is stable, output is unstable.",
                     9: "Input is stable, output is conditionally stable.",
                     12: "Both input and output are stable."}


def p_to_c(ampl, phase, radians=False):
    """
    Converts a complex number on polar form to cartesian complex form.

    Parameters
    ----------
    ampl : float
        The amplitude of the number.
    phase : float
        The phase of the number (in radians or degrees).
    radians : bool
        Whether or not the phase is given in radians.
        Default is False.

    Returns
    -------
    complex
        The complex value in catresian form.
    """
    if (ampl < 0):
        print("Warning: Amplitude is negative. Using absolute value")
        ampl = abs(ampl)
    if (not radians):
        phase = phase*np.pi/180
    return ampl*np.exp(1j*phase)


def fromdb(db_val, use_10=True):
    """
    Converts a value in db to a normal value.

    Parameters
    ----------
    db_val : float
        The value in dB to be converted to normal value.
    use_10 : bool
        True if the db value is defined as 10*log. False if it is defined as
        20*log.
        Default is True.
    """
    if (type(db_val) == complex):
        print("Values in dB can not be complex!")
        return None
    if (use_10):
        return 10**(db_val/10)
    else:
        return 10**(db_val/20)


def todb(val, use_10=True):
    """
    Converts a value to dB.

    Parameters
    ----------
    val : float
        The value to be converted to dB.
    use_10 : bool
        True if the db value should be defined as 10*log. False if should be
        defined as 20*log.
        Default is True.
    """
    if (val < 0):
        print("Values < 0 can not be converted to dB!")
        return None
    if (type(val) == complex):
        print("Complex values can not be converted to dB!")
        return None
    if (use_10):
        return 10*np.log10(val)
    else:
        return 20*np.log10(val)


def input_stability_circle(s_mat, delta):
    """
    Calculates the location of the center of the input stability circle
    as well as the radius if the input stability circle.

    Parameters
    ----------
    s_mat : matrix_like
        A 2x2 matrx containing the S-parameters.
    delta : float
        the delta value

    Returns
    -------
    tuple
        complex
            The center of the input stability circle.
        float
            The radius of the input stability circle.
    """
    r_l = (abs(s_mat[0, 1]*s_mat[1, 0] /
               (abs(s_mat[1, 1])**2 - abs(delta)**2)))
    c_l = ((s_mat[1, 1] - delta*s_mat[0, 0].conjugate()).conjugate() /
           (abs(s_mat[1, 1])**2 - abs(delta)**2))
    return c_l, r_l


def output_stability_circle(s_mat, delta):
    """
    Calculates the location of the center of the output stability circle
    as well as the radius if the output stability circle.

    Parameters
    ----------
    s_mat : matrix_like
        A 2x2 matrx containing the S-parameters.
    delta : float
        the delta value

    Returns
    -------
    tuple
        complex
            The center of the output stability circle.
        float
            The radius of the output stability circle.
    """
    r_s = (abs(s_mat[0, 1]*s_mat[1, 0] /
               (abs(s_mat[0, 0])**2 - abs(delta)**2)))
    c_s = ((s_mat[0, 0] - delta*s_mat[1, 1].conjugate()).conjugate() /
           (abs(s_mat[0, 0])**2 - abs(delta)**2))
    return c_s, r_s


def noise_figure_circle(f, f_min, r_n, z_0, gma_opt):
    """
    Calculates the center and radius of the noise figure circle.

    Parameters
    ----------
    f : float
        Noise figure of the transistor
    f_min : float
        Minimum noise figure of transistor, attained when Y_S = Y_opt.
    r_n : float
        Equivalent noise resistance of transistor.
    z_0 : float
        System impedance.
    gma_opt : complex
        Reflection coefficient resulting from an optimum source impedance that
        results in minimum noise figure.

    Returns
    -------
    tuple
        complex
            The center of the noise figure circle.
        float
            The radius of the noise figure circle.
    """
    noise_fig_param = (f - f_min)*abs(1+gma_opt)**2/(4*r_n/z_0)
    center = gma_opt/(noise_fig_param + 1)
    radius = (np.sqrt(noise_fig_param*(noise_fig_param + 1 - abs(gma_opt)**2))
              / (noise_fig_param + 1))
    return center, radius


def input_const_gain_circle(gma_s, gma_in):
    """
    Calculates the constant gain circle for the input.

    Parameters
    ----------
    gma_s : complex
        Reflection coefficient at the input looking at the input matching
        network.
    gma_in : complex
        Input reflection coefficient. The reflection coefficient seen from
        the input matching network looking at the transistor.
        For the unilateral case this should be equal to S_11.

    Returns
    -------
    tuple
        complex
            The center of the input constant gain circle.
        float
            The radius of the input constant gain circle.
    """
    g_s = (1-abs(gma_s)**2) * (1-abs(gma_in)**2) / abs(1-gma_in*gma_s)**2
    center = g_s*gma_in.conjugate() / (1 - (1-g_s)*abs(gma_in)**2)
    radius = np.sqrt(1-g_s) * (1-abs(gma_in)**2) / (1 - (1-g_s)*abs(gma_in)**2)
    return center, radius


def stability_check(cl, rl, cs, rs):
    """
    Calculates if the values fulfill the stability conditions.

    stable:
    if the real part of Z_in and Z_out is greater than zero for all passive
    load and source impedances.

    conditionally stable:
    if any passive load or source termination can produce input and output
    impedances having a negative real part.

    Parameters
    ----------
    cl : complex
        The location of the center of the input stability circle in the
        complex plane.
    rl : float
        The radius of the input stability circle.
    cs : complex
        The location of the center of the output stability circle in the
        complex plane.
    rs : float
        The radius of the output stability circle.

    Returns
    -------
    int
        The stability mask
        0:  0b0000, both input and output is unstable.
        1:  0b0001, input is unstable, output is conditionally stable.
        2:  0b0010, input is conditionally stable, output is unstable.
        3:  0b0011, both input and output is conditionally stable.
        4:  0b0100, input is unstable, output is stable.
        6:  0b0110, input is conditionally stable, output is stable.
        8:  0b1000, input is stable, output is unstable.
        9:  0b1001, input is stable, output is conditionally stable.
        12: 0b1100, both input and output are stable.
    """
    closest_from_origin_l = abs(cl)-rl
    closest_from_origin_s = abs(cs)-rs
    stable = 0  # unstable
    if (closest_from_origin_l > 1.0):
        # input is unconditionally stable
        stable |= 0b1000
    elif (closest_from_origin_l > -1.0):
        # input is conditionally stable
        stable |= 0x0010
    # else input is unstable
    if (closest_from_origin_s > 1.0):
        # output is unconditionally stable
        stable |= 0b0100
    elif (closest_from_origin_s > -1.0):
        # output is conditionally stable
        stable |= 0b0001
    return stable


def create_circle(center_x, center_y, radius, n_points=101):
    """
    Draws a circle

    Parameters
    ----------
    center_x : float
        The x-coordinate of the center of the circle.
    center_y : float
        The y-coordinate of the center of the circle.
    radius : float
        The radius of the circle.
    n_points : int
        The numer of data points that will be used in the plot.
        Default is 101.

    Returns
    -------
    tuple
        numpy.ndarray
            A numpy array containing the x-coordinates of the circle.
        numpy.ndarray
            A numpy array containing the y-coordinates of the circle.
    """
    theta = np.linspace(0, 2*np.pi, n_points)
    x = center_x + radius*np.cos(theta)
    y = center_y + radius*np.sin(theta)
    return x, y


def plot_circle(ax, center, radius, color=None):
    """
    Plots a circle

    Parameters
    ----------
    ax : matplotlib.axes._subplots.AxesSubplot
        The subplot where the circle will be drawn
    center : complex
        The coordinate of the center of the circle (real, imaginary).
    radius : float
        The radius of the circle.
    color : str
        A hexadecimal colorstring i.e. "#rrggbb".
    """
    x, y = create_circle(np.real(center),
                         np.imag(center),
                         radius,
                         n_points=1000)
    if (color is not None):
        ax.plot(x, y, color=color)
    else:
        ax.plot(x, y)


def get_delta(s_mat):
    """
    Calculates delta (stability)

    Parameters
    ----------
    s_mat : matrix_like
        A 2x2 matrx containing the S-parameters.

    Returns
    -------
    float
        The delta-value
    """
    return s_mat[0, 0]*s_mat[1, 1] - s_mat[0, 1]*s_mat[1, 0]


def get_k(s_mat, delta):
    """
    Calculates K (stability)

    Parameters
    ----------
    s_mat : matrix_like
        A 2x2 matrx containing the S-parameters.
    delta : float
        the delta value

    Returns
    -------
    float
        The K-value
    """
    return (1 - abs(s_mat[0, 0])**2 - abs(s_mat[1, 1])**2 + abs(delta)**2) /\
        (2*abs(s_mat[0, 1]*s_mat[1, 0]))


def get_gma_in(s_mat, gma_ld):
    """
    Calculates gamma_in.

    Parameters
    ----------
    s_mat : matrix_like
        A 2x2 matrx containing the S-parameters.
    gma_ld : complex
        Reflection coefficient at the output looking at the output matching
        network.

    Returns
    -------
    float
        gamma_in
    """
    return s_mat[0, 0] + ((s_mat[0, 1]*s_mat[1, 0]*gma_ld) /
                          (1 - s_mat[1, 1]*gma_ld))


def get_gma_out(s_mat, gma_src):
    """
    Calculates gamma_out.

    Parameters
    ----------
    s_mat : matrix_like
        A 2x2 matrx containing the S-parameters.
    gma_src : complex
        Reflection coefficient at the input looking at the input matching
        network.

    Returns
    -------
    float
        gamma_out
    """
    return s_mat[1, 1] + ((s_mat[0, 1]*s_mat[1, 0]*gma_src) /
                          (1 - s_mat[0, 0]*gma_src))


def get_abs_gma_ports(gma_transistor, gma_matching):
    """
    Calculates the magnitude of the reflection coefficients measured at the
    ports, i.e. before input matching circuit or after output matching circuit.

    Parameters
    ----------
    gma_transistor : complex
        The in/out reflection coefficient for the transistor itself.
        E.g. Gamma_in or Gamma_out.
    gma_matching : complex
        The input/output matching circuit reflection coefficient.
        E.g. Gamma_S or Gamma_L.

    Returns
    -------
    float
        The magnitude of the reflection coefficient seen from the load/source
        looking at the output/input circuit.
    """
    return abs((gma_transistor - gma_matching.conjugate()) /
               (1 - gma_transistor*gma_matching))


def get_noise_figure(f_min, r_n, gma_opt, gma_src, z_0):
    """
    Calculates the noise figure of the transistor.
    See Pozar chapter 12.3 Low amplifier design section.

    Parameters
    ----------
    f_min : float
        Minimum noise figure of transistor, attained when Y_S = Y_opt.
    r_n : float
        Equivalent noise resistance of transistor.
    gma_opt : complex
        Reflection coefficient resulting from an optimum source impedance that
        results in minimum noise figure.
    gma_src : complex
        The source reflection coefficient. The reflection coefficient seen
        from the transistor looking at the input matching network.
    z_0 : float
        System impedance.

    Returns
    -------
    float
        The noise figure of the transistor.
    """
    return f_min + (4*r_n/z_0*abs(gma_src-gma_opt)**2
                    / ((1-abs(gma_src)**2) * abs(1+gma_opt)**2))


def get_vswr(gma):
    """
    Calculates the VSWR (voltage standing wave ratio) for a given reflection
    coefficient.

    Parameters
    ----------
    gma : complex
        The reflection coefficient seen from the point of interest looking
        at the system with this reflection coefficient.

    Returns
    -------
    float
        The VSWR.
        If gma < 1e-16 then this function will return "inf" (as a float).
    """
    gma_abs = abs(gma)
    if (gma_abs < 1e-16):
        return float("inf")
    return (1+gma_abs)/(1-gma_abs)


def get_g_s(gma_in, gma_src):
    """
    Calculates the gain of the input matching circuit.

    Parameters
    ----------
    gma_in : complex
        Input reflection coefficient. The reflection coefficient seen from
        the input matching network looking at the transistor.
    gma_src : complex
        The source reflection coefficient. The reflection coefficient seen
        from the transistor looking at the input matching network.

    Returns
    -------
    float
        The gain of the input matching circuit.
    """
    return (1-abs(gma_src)**2)/abs(1-gma_in*gma_src)**2


def get_g_0(s_mat):
    """
    Calculates the gain of the transistor.

    Parameters
    ----------
    s_mat : matrix_like
        A 2x2 matrx containing the S-parameters.

    Returns
    -------
    float
        The gain of the transistor.
    """
    return abs(s_mat[1, 0])**2


def get_g_l(gma_ld, s_mat):
    """
    Calculates the gain of the output matching circuit.

    Parameters
    ----------
    s_mat : matrix_like
        A 2x2 matrx containing the S-parameters.
    gma_ld the load reflection coefficient. The reflection coefficient seen
        from the transistor looking at the output matching network.

    Returns
    -------
    float
        The gain of the output matching circuit.
    """
    return (1-abs(gma_ld)**2)/abs(1-s_mat[1, 1]*gma_ld)**2


def get_g_t(s_mat, gma_in, gma_src, gma_ld):
    """
    Calculates the transducer power gain.

    Parameters
    ----------
    s_mat : matrix_like
        A 2x2 matrx containing the S-parameters.
    gma_in : complex
        Input reflection coefficient. The reflection coefficient seen from
        the input matching network looking at the transistor.
    gma_src : complex
        The source reflection coefficient. The reflection coefficient seen
        from the transistor looking at the input matching network.
    gma_ld the load reflection coefficient. The reflection coefficient seen
        from the transistor looking at the output matching network.

    Returns
    -------
    float
        The transducer power gain.
    """
    return get_g_s(gma_in, gma_src)*get_g_0(s_mat)*get_g_l(gma_ld, s_mat)


def get_angle(cs, cf):
    """
    Calculates the angle between the points [0, 0], [real(cs), imag(cs)],
    [real(cf), imag(cf)] using expression:

    cs = cf + abs(cf-cs)*exp(1j*theta)

    Parameters
    ----------
    cs : complex
        Center of gain circle.
    cf : complex
        Center of noise circle.

    Returns
    -------
    float
        The angle between the points.
    """
    return np.real(-1j*np.log((cf-cs)/abs(cf-cs)))


def optim_gain_noise_helpfun(
        gma_s, s_mat, f_targ, f_min, r_n, z_0, gma_opt, abs_gma_out_targ):
    """
    This function performs the calculations for the optimize_gma_s function
    and returns the target function evaluated with the supplied parameters.

    Parameters
    ----------
    gmma_s: complex
        The source reflection coefficient. The reflection coefficient seen
        from the transistor looking at the input matching network.
    s_mat : matrix_like
        A 2x2 matrx containing the S-parameters.
    f_targ : dict
        The target noise figure of the transistor.
    f_min : float
        Minimum noise figure of transistor, attained when Y_S = Y_opt.
    r_n : float
        Equivalent noise resistance of transistor.
    z_0 : float
        System impedance.
    gma_opt : complex
        Reflection coefficient resulting from an optimum source impedance that
        results in minimum noise figure.
    abs_gma_out_targ : float
        The maximum absolute value of gma_out allowed by the VSWRout
        constraints.

    Returns
    -------
    tuple
        float
            The real part of the target function. This function should be
            linked to the real part of gamma_s.
        float
            The imaginary part of then target function. This function should
            be liked to the imaginary part of gamma_s.
    """
    gma_out = get_gma_out(s_mat, gma_s)
    gma_l = gma_out.conjugate()
    gma_in = get_gma_in(s_mat, gma_l)

    cf, rf = noise_figure_circle(
        f=f_targ, f_min=f_min, r_n=r_n, z_0=50, gma_opt=gma_opt)
    cs, rs = input_const_gain_circle(gma_s, gma_in)
    v = get_angle(cs, cf)
    out = ((cf-cs) - (rf+rs)*np.exp(1j*v))/abs(cf-cs)
    return np.real(out), np.imag(out)


def optimize_gma_s(
        gma_s_start, s_mat, vswr_out_targ, f_targ, f_min, r_n, z_0, gma_opt):
    """
    This function runs the equation solving code for gamma_s. It attempts to
    optimize gamma_s such that the constant noise circle and the constant gain
    circles are separated by the sum of the radius'. This means that they
    intersect in one point which will be the final value of gamma_s.

    Parameters
    ----------
    gmma_s_start: complex
        Initial value of gamma_s, the source reflection coefficient. The
        reflection coefficient seen from the transistor looking at the input
        matching network.
    s_mat : matrix_like
        A 2x2 matrx containing the S-parameters.
    vswr_out_targ : float
        The VSWRout target value.
    f_targ : float
        The noise figure target value.
    f_min : float
        Minimum noise figure of transistor, attained when Y_S = Y_opt.
    r_n : float
        Equivalent noise resistance of transistor.
    z_0 : float
        System impedance.
    gma_opt : complex
        Reflection coefficient resulting from an optimum source impedance that
        results in minimum noise figure.

    Returns
    -------
    complex
        The optimized gamma_s.
    """
    abs_gma_out_targ = (vswr_out_targ-1)/(vswr_out_targ+1)
    gma_s = spo.fsolve(lambda gma_s: optim_gain_noise_helpfun(
        gma_s[0]+1j*gma_s[1], s_mat, f_targ, f_min, r_n, z_0, gma_opt,
        abs_gma_out_targ),
                       [np.real(gma_s_start), np.imag(gma_s_start)])
    return gma_s[0]+1j*gma_s[1]


def optimize_vswr_helpfun(abs_gma_sys_targ, gma_transistor, gma_matching):
    """
    Calculates the residual of the magnitude of the reflection coefficient of
    the system supplied to then function and calclated using the remaining
    parameters supplied to the function.

    Parameters
    ----------
    abs_gma_sys_targ : float
        Target system reflection coefficient amplitude seen from the
        source/load looking at the input/output matching network.
    gma_transistor : complex
        The reflection coefficient of the transistor, e.g. gamma_in or
        gamma_out.
    gma_matching_start : complex
        Starting value of the reflection coefficient of the matching
        network.

    Returns
    -------
    tuple
        float
            The magnitude of the residual.
        float
            The magnitude of the residual.
    """
    out = (abs_gma_sys_targ
           - get_abs_gma_ports(gma_transistor, gma_matching))
    return out, out


def optimize_vswr_matching(vswr_targ, gma_transistor, gma_matching_start):
    """
    Optimizes the reflection coefficient of the matching circuit to make the
    reflection coefficient of the system equal to the target reflection
    coefficient.

    Parameters
    ----------
    vswr_targ : float
        Target voltage standing wave ratio.
    gma_transistor : complex
        The reflection coefficient of the transistor, e.g. gamma_in or
        gamma_out.
    gma_matching_start : complex
        Starting value of the reflection coefficient of the matching
        network.

    Returns
    -------
    complex
        The new reflection coefficient for the matching circuit.
    """
    abs_gma_sys_targ = (vswr_targ - 1)/(vswr_targ + 1)
    gma_matching = spo.fsolve(
        lambda gma_matching: optimize_vswr_helpfun(
            abs_gma_sys_targ, gma_transistor,
            gma_matching[0]+1j*gma_matching[1]),
        [np.real(gma_matching_start), np.imag(gma_matching_start)])
    return gma_matching[0]+1j*gma_matching[1]


"""
def optimize_gain_helpfun(s_mat, gma_src, gma_ld, g_t_targ):
    gma_in = get_gma_in(s_mat, gma_ld)
    out = g_t_targ - get_g_t(s_mat, gma_in, gma_src, gma_ld)
    print(out)
    return (out, out, out, out)
"""


def optimize_gain(s_mat, gma_src, gma_ld, g_t_targ):
    """
    Optimizes then transducer gain.

    Parameters
    ----------
    s_mat : matrix_like
        A 2x2 matrx containing the S-parameters.
    gma_in : complex
        Input reflection coefficient. The reflection coefficient seen from
        the input matching network looking at the transistor.
    gma_src : complex
        The source reflection coefficient. The reflection coefficient seen
        from the transistor looking at the input matching network.
    gma_ld the load reflection coefficient. The reflection coefficient seen
        from the transistor looking at the output matching network.

    Returns
    -------
    float
        The optimizes transducer power gain.
    """
    gma_matching = spo.fsolve(
        lambda gma_matching: optimize_gain_helpfun(
            s_mat, gma_matching[0]+1j*gma_matching[1],
            gma_matching[2]+1j*gma_matching[3], g_t_targ),
        [np.real(gma_ld), np.imag(gma_ld), np.real(gma_src), np.imag(gma_src)])
    return gma_matching[0]+1j*gma_matching[1], gma_matching[2]+1j*gma_matching[3]
