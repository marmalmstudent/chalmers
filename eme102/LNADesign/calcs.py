import numpy as np


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
        The value in dB.
    use_10 : bool
        True if the db value is defined as 10*log. False if it is defined as
        20*log.
        Default is True.
    """
    if (use_10):
        return 10**(db_val/10)
    else:
        return 10**(db_val/20)


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
        gamma_in
    """
    return s_mat[1, 1] + ((s_mat[0, 1]*s_mat[1, 0]*gma_src) /
                          (1 - s_mat[0, 0]*gma_src))


def input_stability_circle(s_mat, delta):
    """
    --------------------------------------------------------------------
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
    --------------------------------------------------------------------
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


def stability_check(k, delta):
    """
    Calculates if the values fulfill the stability conditions.

    Unconditionally stable:
    if the real part of Z_in and Z_out is greater than zero for all passive
    load and source impedances.

    Potentially unstable:
    if any passive load or source termination can produce input and output
    impedances having a negative real part.

    Parameters
    ----------
    k : float
        The K-value
    delta : float
        The delta-value

    Returns
    -------
    int
        The stability mask
    """
    stable = 0  # unstable
    if (False):
        stable = 1  # conditionally stable
    if (k > 1 and abs(delta) < 1):
        stable = 2  # unconditionally stable
    return stable


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
    s = np.array([[p_to_c(0.93, -60.1), p_to_c(0.028, 57.7)],
                  [p_to_c(1.61, 128.7), p_to_c(0.43, -172.3)]],
                 dtype=np.complex128)
    gma_opt = p_to_c(0.78, 58)
    f_min = fromdb(0.7)
    r_n = 22

    lna_cond = {"f": fromdb(0.9), "vswr_in": 2,
                "vswr_out": 1.9, "g_t": fromdb(14)}
