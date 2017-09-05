import numpy as np


def task2():
    m = 0.010  # [kg] mass of silver plate
    phi = 1e10  # [neutrons/(m^2 s)] neutron flux density
    N_a = 6.022e23  # Avogadro's constant
    M_Ag = 0.10787  # [kg/mol] molar mass of Ag (average over all isotopes)
    rho_Ag = 10.49e3  # [kg/m^3] density of Ag
    prcnt_iso = np.array([0.518, 0.482])  # [Ag107, Ag109]
    s_iso = np.array([35e-28, 89e-28])  # thrml absrp cross sctn, Ag iso [m^2]
    v = m/rho_Ag  # volume of Ag plate
    # N = rho_Ag/M_Ag*N_a  # nbr scatt Ag nuclei/m^3
    N_iso = prcnt_iso*rho_Ag/M_Ag*N_a  # nbr scatt Ag isotope nuclei/m^3
    molar_new_iso_rate = phi*s_iso*v*N_iso/N_a
    return molar_new_iso_rate


def task3():
    """
    Chemical formula for paraffin: C_j H_2j+2
    """
    E_0 = 5e6  # [eV] initial neutron energy
    E_n = 25e-3  # [eV] final neutron energy
    n = np.ceil(np.log2(E_0/E_n))  # nbr of collisions to reach E_n
    
    N_A = 6.022e23  # Avogadro's constant
    s = 20e-28  # [m^2]
    j = 58  # The paraffin used
    n_H = 2*j+2  # number of hydrogen atom per paraffin molecule
    rho_paraffin = 900  # [kg/m^3] density of paraffin
    M_paraffin = (58*12.0107+118*1.00794)/1000  # [kg/mol] of paraffin
    # N_parraf = rho_paraffin/M_paraffin*N_A  # nbr of scattering paraff/m^3
    N_H = n_H*rho_paraffin/M_paraffin*N_A  # nbr of scattering hydrogen/m^3
    sigma = N_H*s  # macroscopic cross section
    l = 1/sigma  # mean free path

    dx_ub = n*l  # lower bound of paraffin thickness
    return np.array([n, l, dx_ub])


def task4():
    # activation time of 90 % saturation
    max_hl = 180  # [s]
    lbda = np.log(2)/max_hl
    t_90 = -np.log(0.1)/lbda
    return t_90

print task2()
print task3()
print task4()
