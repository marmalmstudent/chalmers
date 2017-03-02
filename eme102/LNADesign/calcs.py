import numpy as np
import scipy.optimize as spo
import matplotlib.pyplot as plt


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
    return s_mat[0, 0] + (s_mat[0, 1]*s_mat[1, 0]*gma_ld) /\
        (1 - s_mat[1, 1]*gma_ld)


def func_in(s_mat, gma_vec):
    """
    Calculates the conjugate matching error for the input.

    Parameters
    ----------
    s_mat : matrix_like
        A 2x2 matrx containing the S-parameters.
    gma_vec : array_like
        An array containing the real and imaginary part of gamma_in.

    Returns
    -------
    tuple
        float
            The real part of the residual.
        float
            The imaginary part of then residual.
    """
    gma_in = gma_vec[0] + 1j*gma_vec[1]
    ret = gma_in - get_gma_in(s_mat, gma_in)
    return np.real(ret), np.imag(ret)


def get_gma_in_init(s_mat):
    """
    Calculates gamma_in under the assumption that gamma_load and gamma_in
    are conjugate matched (gamma_in = gamma_load*).

    Parameters
    ----------
    s_mat : matrix_like
        A 2x2 matrx containing the S-parameters.

    Returns
    -------
    float
        gamma_in under the conjugate matching assumption.
    """
    gma_in = spo.fsolve(lambda gma_in: func_in(s_mat, gma_in), [0, 0])
    return (gma_in[0] + 1j*gma_in[1])


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


def func_out(s_mat, gma_vec):
    """
    Calculates the conjugate matching error for the output.

    Parameters
    ----------
    s_mat : matrix_like
        A 2x2 matrx containing the S-parameters.
    gma_vec : array_like
        An array containing the real and imaginary part of gamma_out.

    Returns
    -------
    tuple
        float
            The real part of the residual.
        float
            The imaginary part of then residual.
    """
    gma_out = gma_vec[0] + 1j*gma_vec[1]
    ret = gma_out - get_gma_out(s_mat, gma_out)
    return np.real(ret), np.imag(ret)


def get_gma_out_init(s_mat):
    """
    Calculates gamma_out under the assumption that gamma_load and gamma_out
    are conjugate matched (gamma_out = gamma_load*).

    Parameters
    ----------
    s_mat : matrix_like
        A 2x2 matrx containing the S-parameters.

    Returns
    -------
    float
        gamma_out under the conjugate matching assumption.
    """
    gma_out = spo.fsolve(lambda gma_out: func_out(s_mat, gma_out), [0, 0])
    return (gma_out[0] + 1j*gma_out[1])


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


def get_g_t(s_mat, gma_in, gma_src, gma_ld):
    """
    Calculates then transducer power gain.

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
    return (1-abs(gma_src)**2)/abs(1-gma_in*gma_src)**2 *\
        abs(s_mat[1, 0])**2 *\
        (1-abs(gma_ld)**2)/abs(1-s_mat[1, 1]*gma_ld)**2


def func_optimize(conds, s_mat, gma_in, gma_out, gma_src, gma_ld):
    """
    This function connects all of the functions and the variables that
    controls the amplifier specifications.

    Idea: perform an fmin on the function:
        abs(conds["F"]-f, conds["VSWRin"]-vswr_in,
            conds["VSWRout"]-vswr_out, conds["GT"]-g_t).
    This function must also meed the requirements:
        f <= conds["F"]
        vswr_in <= conds["VSWRin"]
        vswr_out <= conds["VSWRout"]
        g_t <= conds["GT"]
    This could be done by for example:
        if (vswr_in > conds["VSWRin"]):
            vswr_in = conds["VSWRin"] - abs(conds["VSWRin"]-vswr_in)
            gma_in = get_gma(vswr_in)
    However, this might be difficult for parameters with multiple input
    parameters.
    """
    # F here

    # VSWRin here
    vswr_in = get_vswr(gma_in)
    # VSWRout
    vswr_out = get_vswr(gma_out)
    # GT
    g_t = get_g_t(s_mat, gma_in, gma_src, gma_ld)
    return vswr_in, vswr_out, g_t


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


def plot_stab_circle(ax, center, radius, color=None):
    """
    Plots stability circles

    Parameters
    ----------
    ax : matplotlib.axes._subplots.AxesSubplot
        The subplot where the circle will be drawn
    center : complex
        The coordinate of the center of the circle (real, imaginary).
    radius : float
        The radius of then circle.
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


if (__name__ == '__main__'):
    """
    LNA design -> available gain: G_a is important.

    -------------------------- Stability circles ------------------------------
    - Stable if |gma_in|<1 and |gma_out|<1;
    - How can we check this?
    - Calculate gma_s and gma_l for |gma_out|=1 and |gma_in|=1.
    - Ex. gma_l values resulting in |gma_in|=1 => stability circles for the
      output matching network.

    ------------------------ Potentially unstable -----------------------------
    - Calculated circle provides the boundary between stable and unstable
      regions.
    - Which part is stable?
    - Simple check: assume gma_l = 0 => gma_l=s_11. If |s_11|<1, the centre of
      the Smith chart represents a stable operating point.

    ----------------------- Unconditionally stable ----------------------------
    - |gma_in|<1; |gma_out|<1; |gma_s|<1; and |gma_l|<1;
    - From a graphical view - we would like to see the stability circles to
      fall completely outside the SC.

    -------------------------- Stability circles ------------------------------
    - Input:
        The intersection between |s_11| < 1 and ||c_l| - r_l | > 1:
            If gma_in < 1 for |gma_l| = 0 then then intersection represents
            that we can have |gma_in| > 1 for some |s_11| < 1.
            Conditionally unstable/Conditionally stable.
            Typically larger stable region.

            If gma_in > 1 for |gma_l| = 0 then then intersection represents
            that we can have |gma_in| < 1 for some |s_11| < 1.
            Conditionally stable/Conditionally unstable.
            Typically smaller stable region.
    - Output:
        The intersection between |s_22| < 1 and ||c_s| - r_s | > 1:
            If gma_out < 1 for |gma_s| = 0 then then intersection represents
            that we can have |gma_out| > 1 for some |s_22| < 1.
            Conditionally unstable/Conditionally stable.
            Typically larger stable region.

            If gma_out > 1 for |gma_s| = 0 then then intersection represents
            that we can have |gma_out| < 1 for some |s_22| < 1.
            Conditionally stable/Conditionally unstable.
            Typically smaller stable region.

    """
    lna_cond = {"F": fromdb(0.9), "VSWRin": 2,
                "VSWRout": 1.9, "GT": fromdb(14)}
    f = 1.8e9
    s = np.array([[p_to_c(0.93, -60.1), p_to_c(0.028, 57.7)],
                  [p_to_c(1.61, 128.7), p_to_c(0.43, -172.3)]],
                 dtype=np.complex128)
    gma_opt = p_to_c(0.78, 58)
    f_min = fromdb(0.7)
    r_n = 22

    # stability check
    delta = get_delta(s_mat=s)
    cl, rl = input_stability_circle(s_mat=s, delta=delta)
    cs, rs = output_stability_circle(s_mat=s, delta=delta)
    stab = stability_check(cl, rl, cs, rs)
    print(__stability_codes[stab])

    # plot stability circles
    fig = plt.figure()
    ax0 = fig.add_subplot(121)
    ax0.set_aspect('equal')
    plot_stab_circle(ax0, 0, 1, color="#0000ff")
    plot_stab_circle(ax0, cl, rl, color="#ff0000")
    ax1 = fig.add_subplot(122)
    ax1.set_aspect('equal')
    plot_stab_circle(ax1, 0, 1, color="#0000ff")
    plot_stab_circle(ax1, cs, rs, color="#ff0000")

    # initial values form gma_s and gma_l, seen lecture 8 slide 6.
    gma_s = gma_opt
    gma_out = get_gma_out(s, gma_s)
    gma_l = gma_out.conjugate()
    gma_in = get_gma_out(s, gma_l)
    """
    print(gma_s, gma_out, gma_l, gma_in)
    print(get_g_t(s, gma_in, gma_s, gma_l))
    print(get_vswr(gma_in), get_vswr(gma_out))
    """



    # plt.show()
