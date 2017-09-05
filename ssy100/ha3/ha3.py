import scipy as sp
import scipy.io as spi
import numpy as np


s_par = sp.io.loadmat("S_par.mat")
header = s_par['__header__']
version = s_par['__version__']
global_val = s_par['__globals__']
s_parameter = s_par['S_parameter'][0][0]
s_par_data_0 = s_parameter[0]
s_par_data_1 = s_parameter[1]
s_par_data_2 = s_parameter[2]
s_par_data_3 = s_parameter[3]
s_par_data_4 = s_parameter[4]
s_par_data_5 = s_parameter[5]
s_par_data_6 = s_parameter[6]
s_par_data_7 = s_parameter[7]
s_par_freq = s_parameter[8]
s_par_cplx = s_parameter[9][0]
s_par_cplx_0 = s_par_cplx[0][0][0][0]
s_par_cplx_1 = s_par_cplx[1][0]
print(s_par_cplx_0)
# print(s_par_cplx_1)
