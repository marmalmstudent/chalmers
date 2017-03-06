import numpy as np
import scipy.optimize as spo


z_0 = 50
r_l = 50
f = 1.8e9
c = 3e8
lbda = c/f


def funcsolve(z_targ, bd, bl, short_term=True):
    if (short_term):
        z_l = z_0*np.tanh(1j*bl)
    else:
        z_l = z_0/np.tanh(1j*bl)
    z_im = r_l*z_l/(r_l+z_l)
    z_d = z_0*((z_im + z_0*np.tanh(1j*bd))
               / (z_0 + z_im*np.tanh(1j*bd)))
    out = abs(z_targ - z_d)
    return out, out


gamma_s = 0.462+0.762j
gamma_l = -0.454 + 0.065j
z_s = z_0*(1+gamma_s)/(1-gamma_s)
z_l = z_0*(1+gamma_l)/(1-gamma_l)
lengths_l_o = abs(spo.fsolve(lambda lengths_l_o: funcsolve(
    z_l, abs(lengths_l_o[0]), abs(lengths_l_o[1]), short_term=False),
                     [0.1, 0.1])) % np.pi
print((lengths_l_o % np.pi)*180/np.pi)
"""
lengths_l_s = abs(spo.fsolve(lambda lengths_l_s: funcsolve(
    z_l, abs(lengths_l_s[0]), abs(lengths_l_s[1]), short_term=True),
                     [0.1, 0.1])) % np.pi
l_o_diff = sum(abs(lengths_l_o - np.average(lengths_l_o)))
l_s_diff = sum(abs(lengths_l_s - np.average(lengths_l_s)))
if (l_o_diff > l_s_diff):
    print("Output matching network should use short stub with electrical length %.3f degrees closest to the load and a transmission line with electrical length %.3f degrees" % (lengths_l_s[1]*180/np.pi, lengths_l_s[0]*180/np.pi))
else:
    print("Output matching network should use open stub with electrical length %.3f degrees closest to the load and a transmission line with electrical length %.3f degrees" % (lengths_l_o[1]*180/np.pi, lengths_l_o[0]*180/np.pi))
"""

lengths_s_s = abs(spo.fsolve(lambda lengths_s_s: funcsolve(
    z_s, abs(lengths_s_s[0]), abs(lengths_s_s[1]), short_term=True),
                     [0.1, 0.1])) % np.pi
print((lengths_s_s % np.pi)*180/np.pi)
"""
lengths_s_o = abs(spo.fsolve(lambda lengths_s_o: funcsolve(
    z_s, abs(lengths_s_o[0]), abs(lengths_s_o[1]), short_term=False),
                     [0.1, 0.1])) % np.pi
s_o_diff = sum(abs(lengths_s_o - np.average(lengths_s_o)))
s_s_diff = sum(abs(lengths_s_s - np.average(lengths_s_s)))
if (s_o_diff > s_s_diff):
    print("Input matching network should use short stub with electrical length %.3f degrees closest to the load and a transmission line with electrical length %.3f degrees" % (lengths_s_s[1]*180/np.pi, lengths_s_s[0]*180/np.pi))
else:
    print("Input matching network should use open stub with electrical length %.3f degrees closest to the load and a transmission line with electrical length %.3f degrees" % (lengths_s_o[1]*180/np.pi, lengths_s_o[0]*180/np.pi))
"""
