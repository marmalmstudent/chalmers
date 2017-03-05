import numpy as np
import matplotlib.pyplot as plt
import calcs


if (__name__ == '__main__'):
    lna_cond = {"F": calcs.fromdb(0.9), "VSWRin": 2,
                "VSWRout": 1.9, "GT": calcs.fromdb(14.0)}
    f = 1.8e9
    s = np.array([[calcs.p_to_c(0.93, -60.1), calcs.p_to_c(0.028, 57.7)],
                  [calcs.p_to_c(1.61, 128.7), calcs.p_to_c(0.43, -172.3)]],
                 dtype=np.complex128)
    gma_opt = calcs.p_to_c(0.78, 58)
    f_min = calcs.fromdb(0.7)
    r_n = 22

    # stability check
    delta = calcs.get_delta(s_mat=s)
    cl, rl = calcs.input_stability_circle(s_mat=s, delta=delta)
    cs, rs = calcs.output_stability_circle(s_mat=s, delta=delta)
    stab = calcs.stability_check(cl, rl, cs, rs)
    print(calcs.__stability_codes[stab])

    fig = plt.figure()
    ax0 = fig.add_subplot(111)
    ax0.set_aspect('equal')
    cf, rf = calcs.noise_figure_circle(
        f=lna_cond["F"], f_min=f_min, r_n=r_n, z_0=50, gma_opt=gma_opt)
    calcs.plot_circle(ax0, 0, 1, color="#000000")
    calcs.plot_circle(ax0, cf, rf, color="#ff0000")
    """
    # plot stability circles
    fig = plt.figure()
    ax0 = fig.add_subplot(121)
    ax0.set_aspect('equal')
    calcs.plot_circle(ax0, 0, 1, color="#0000ff")
    calcs.plot_circle(ax0, cl, rl, color="#ff0000")
    ax1 = fig.add_subplot(122)
    ax1.set_aspect('equal')
    calcs.plot_circle(ax1, 0, 1, color="#0000ff")
    calcs.plot_circle(ax1, cs, rs, color="#ff0000")
    """

    # find gamma_s such that constant gain circle and constant noise circle
    # intersect.
    gma_s = calcs.optimize_gma_s(
        gma_opt, s, lna_cond["VSWRout"],
        lna_cond["F"], f_min, r_n, 50, gma_opt)
    gma_out = calcs.get_gma_out(s, gma_s)
    gma_l = gma_out.conjugate()
    gma_in = calcs.get_gma_in(s, gma_l)

    # calculate gamma_s and gamma_l that satisfies vswr constraints.
    count = 0
    while (((calcs.get_vswr(calcs.get_abs_gma_ports(gma_out, gma_l))
             > lna_cond["VSWRout"])
            or (calcs.get_vswr(calcs.get_abs_gma_ports(gma_in, gma_s))
                > lna_cond["VSWRin"]))
           and (count < 100)):
        if (calcs.get_vswr(calcs.get_abs_gma_ports(gma_out, gma_l))
                > lna_cond["VSWRout"]):
            gma_l = calcs.optimize_vswr_matching(
                lna_cond["VSWRout"], gma_out, gma_l)
            gma_in = calcs.get_gma_in(s, gma_l)
        if (calcs.get_vswr(calcs.get_abs_gma_ports(gma_in, gma_s))
                > lna_cond["VSWRin"]):
            gma_s = calcs.optimize_vswr_matching(
                lna_cond["VSWRin"], gma_in, gma_s)
            gma_out = calcs.get_gma_out(s, gma_s)
        if (calcs.get_g_t(s, gma_in, gma_s, gma_l) < lna_cond["GT"]):
            gma_s, gma_l = calcs.optimize_gain(s, gma_s, gma_l, lna_cond["GT"])
            gma_in = calcs.get_gma_in(s, gma_l)
        count += 1

    """
    gma_s = calcs.optimize_gma_s(gma_opt, s, lna_cond, f_min, r_n, 50, gma_opt)
    gma_out = calcs.get_gma_out(s, gma_s)
    """
    """
    gma_l = calcs.optimize_vswr_matching(
        lna_cond["VSWRout"], gma_out, gma_l)
    gma_in = calcs.get_gma_in(s, gma_l)
    """
    print(calcs.todb(calcs.get_noise_figure(f_min, r_n, gma_opt, gma_s, 50)))
    print(calcs.get_abs_gma_ports(gma_out, gma_l),
          calcs.get_vswr(calcs.get_abs_gma_ports(gma_out, gma_l)))
    print(calcs.get_abs_gma_ports(gma_in, gma_s),
          calcs.get_vswr(calcs.get_abs_gma_ports(gma_in, gma_s)))
    print(calcs.todb(calcs.get_g_t(s, gma_in, gma_s, gma_l)))

    cs, rs = calcs.input_const_gain_circle(gma_s, gma_in)
    calcs.plot_circle(ax0, cs, rs, color="#00ff00")
    """
    gma_s = gma_opt*1.156*0.9
    gma_out = calcs.get_gma_out(s, gma_s)
    gma_l = gma_out.conjugate()
    gma_in = calcs.get_gma_in(s, gma_l)

    cs, rs = calcs.input_const_gain_circle(gma_s, gma_in)
    calcs.plot_circle(ax0, cs, rs, color="#00ff00")
    print(calcs.todb(calcs.get_g_s(gma_in, cs+rs)))
    print(calcs.todb(calcs.get_noise_figure(f_min, r_n, gma_opt, cf+rf, 50)))
    vec1 = np.array([np.real(cf) - np.real(cs),
                     np.imag(cf) - np.imag(cs)])
    vec2 = np.array([np.real(cs)+1 - np.real(cs),
                     np.imag(cs) - np.imag(cs)])
    v = -np.arccos(np.dot(vec1, vec2)
                   / (np.linalg.norm(vec1)*np.linalg.norm(vec2)))
    print(360 + v*180/np.pi)
    print(calcs.todb(calcs.get_noise_figure(
        f_min, r_n, gma_opt, cs+rs*(np.cos(v)+1j*np.sin(v)), 50)))
    ax0.plot(np.real(cs+rs*(np.cos(v)+1j*np.sin(v))),
             np.imag(cs+rs*(np.cos(v)+1j*np.sin(v))), marker='o')
    ax0.plot(np.real(cs),
             np.imag(cs), marker='o')
    """
    """
    print(noise_figure_to_gma_s(f=lna_cond["F"], f_min=f_min, r_n=r_n,
                                gma_opt=gma_opt, z_0=50))
    """
    """
    print(func_optimize(conds=lna_cond, f_min=f_min, r_n=r_n, s_mat=s,
                        gma_in=gma_in, gma_out=gma_out, gma_src=gma_s,
                        gma_ld=gma_l, z_0=50))
    """
    """
    print(gma_s, gma_out, gma_l, gma_in)
    print(get_g_t(s, gma_in, gma_s, gma_l))
    print(get_vswr(gma_in), get_vswr(gma_out))
    """



    plt.show()
