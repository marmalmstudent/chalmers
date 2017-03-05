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
    z_0 = 50

    # stability check
    delta = calcs.get_delta(s_mat=s)
    cl, rl = calcs.input_stability_circle(s_mat=s, delta=delta)
    cs, rs = calcs.output_stability_circle(s_mat=s, delta=delta)
    stab = calcs.stability_check(cl, rl, cs, rs)
    print(calcs.__stability_codes[stab])

    fig = plt.figure()
    ax0 = fig.add_subplot(111)
    ax0.set_aspect('equal')
    ax0.axis([-1.1, 1.1, -1.1, 1.1])

    calcs.plot_circle(ax0, 0, 1)
    calcs.plot_circle(ax0, cl, rl, linestyle='--')
    calcs.plot_circle(ax0, cs, rs, linestyle='--')

    # find gamma_s such that constant gain circle and constant noise circle
    # intersect.
    gma_s = gma_opt
    gma_out = calcs.get_gma_out(s, gma_s)
    gma_l = gma_out.conjugate()
    gma_in = calcs.get_gma_in(s, gma_l)

    # calculate gamma_s and gamma_l that satisfies the constraints.
    count = 0
    _max_iter = 1000
    while (((calcs.get_vswr(calcs.get_abs_gma_ports(gma_out, gma_l))
             > lna_cond["VSWRout"])
            or (calcs.get_vswr(calcs.get_abs_gma_ports(gma_in, gma_s))
                > lna_cond["VSWRin"])
            or (calcs.get_noise_figure(f_min, r_n, gma_opt, gma_s, z_0)
                > lna_cond["F"])
            or (calcs.get_g_t(s, gma_in, gma_s, gma_l)
                < lna_cond["GT"]))
           and (count < _max_iter)):
        gma_l = calcs.optimize_vswr_matching(
            lna_cond["VSWRout"], gma_out, gma_l)
        gma_in = calcs.get_gma_in(s, gma_l)
        gma_s = calcs.optimize_vswr_matching(
            lna_cond["VSWRin"], gma_in, gma_s)
        gma_out = calcs.get_gma_out(s, gma_s)
        gma_s = calcs.optimize_gma_s(
            gma_s, lna_cond["F"], f_min, r_n, z_0, gma_opt, gma_in)
        gma_out = calcs.get_gma_out(s, gma_s)
        """ This needs fixing
        if (calcs.get_noise_figure(f_min, r_n, gma_opt, gma_s, z_0)
                > lna_cond["F"]):
            print("Fixing F!")
            gma_s = calcs.optimize_gma_s(
                gma_opt, s, lna_cond["F"], f_min, r_n, z_0, gma_opt, gma_in)
            gma_out = calcs.get_gma_out(s, gma_s)
        if (calcs.get_vswr(calcs.get_abs_gma_ports(gma_out, gma_l))
                > lna_cond["VSWRout"]):
            print("Fixing VSWRout!")
            gma_l = calcs.optimize_vswr_matching(
                lna_cond["VSWRout"], gma_out, gma_l)
            gma_in = calcs.get_gma_in(s, gma_l)
        if (calcs.get_vswr(calcs.get_abs_gma_ports(gma_in, gma_s))
                > lna_cond["VSWRin"]):
            print("Fixing VSWRin!")
            gma_s = calcs.optimize_vswr_matching(
                lna_cond["VSWRin"], gma_in, gma_s)
            gma_out = calcs.get_gma_out(s, gma_s)
        # Gain calculations do not work yet
        if (calcs.get_g_t(s, gma_in, gma_s, gma_l) < lna_cond["GT"]):
            print("Fixing GT!")
            gma_s, gma_l = calcs.optimize_gain(s, gma_s, gma_l, lna_cond["GT"])
            gma_in = calcs.get_gma_in(s, gma_l)
        """
        print("%d:\t%.3f\t%.3f\t%.3f\t%.3f" % (
            count, calcs.todb(calcs.get_noise_figure(
                f_min, r_n, gma_opt, gma_s, z_0)),
            calcs.get_vswr(calcs.get_abs_gma_ports(gma_out, gma_l)),
            calcs.get_vswr(calcs.get_abs_gma_ports(gma_in, gma_s)),
            calcs.todb(calcs.get_g_t(s, gma_in, gma_s, gma_l))))
        count += 1
    if (count < _max_iter):
        F = calcs.get_noise_figure(f_min, r_n, gma_opt, gma_s, z_0)
        abs_gma_b = calcs.get_abs_gma_ports(gma_out, gma_l)
        VSWRout = calcs.get_vswr(abs_gma_b)
        abs_gma_a = calcs.get_abs_gma_ports(gma_in, gma_s)
        VSWRin = calcs.get_vswr(abs_gma_a)
        GT = calcs.get_g_t(s, gma_in, gma_s, gma_l)
        print("\n\n\033[92m", end="")
        print("########################################"
              + "########################################")
        print("#                       Program terminated succe"
              + "ssfully                        #")
        print("########################################"
              + "########################################")
        print("# Number of iterations:% 3d              " % count
              + "                                       #")
        print("# Noise Figure:% 15.3f dB                " % calcs.todb(F)
              + "                              #")
        print("# VSWRout:% 20.3f                        " % VSWRout
              + "                         #")
        print("# VSWRin:% 21.3f                         " % VSWRin
              + "                        #")
        print("# GT:% 25.3f dB                          " % calcs.todb(GT)
              + "                    #")
        print("#                                       "
              + "                                       #")
        print("# Gamma_S:% 20.3f + % 4.3fj              " % (
            np.real(gma_s), np.imag(gma_s))
              + "                         #")
        print("# Gamma_out:% 18.3f + % 4.3fj            " % (
            np.real(gma_out), np.imag(gma_out))
              + "                           #")
        print("# Gamma_L:% 20.3f + % 4.3fj              " % (
            np.real(gma_l), np.imag(gma_l))
              + "                         #")
        print("# Gamma_in:% 19.3f + % 4.3fj             " % (
            np.real(gma_in), np.imag(gma_in))
              + "                          #")
        print("########################################"
              + "########################################")

        cvi = (gma_in.conjugate()*(1-abs_gma_a**2)
               / (1 - (abs_gma_a*abs(gma_in))**2))
        rvi = (abs_gma_a*(1-abs(gma_in)**2)
               / (1 - (abs_gma_a*abs(gma_in))**2))
        calcs.plot_circle(ax0, cvi, rvi, color="#00ff00")

        cvo = (gma_out.conjugate()*(1-abs_gma_b**2)
               / (1 - (abs_gma_b*abs(gma_out))**2))
        rvo = (abs_gma_b*(1-abs(gma_out)**2)
               / (1 - (abs_gma_b*abs(gma_out))**2))
        calcs.plot_circle(ax0, cvo, rvo, color="#00ff00")

        cf, rf = calcs.noise_figure_circle(
            f=calcs.get_noise_figure(f_min, r_n, gma_opt, gma_s, z_0),
            f_min=f_min, r_n=r_n, z_0=z_0, gma_opt=gma_opt)
        calcs.plot_circle(ax0, cf, rf, color="#0000ff")
        cs, rs = calcs.input_const_gain_circle(gma_s, gma_in)
        calcs.plot_circle(ax0, cs, rs, color="#ff0000", linestyle=':')
    else:
        print("\n\n\033[91m", end="")
        print("########################################"
              + "########################################")
        print("#                      Program terminated unsucce"
              + "ssfully                       #")
        print("########################################"
              + "########################################")
    print("\033[0m")  # normal text color
    if (count < 1000):
        plt.show()
