f = 4e4;
w = 2*pi*f;
sigma = 4;
mu_0 = 4*pi*1e-7;
e_0 = 8.85e-12;
e_r = 80;
e_d = e_r*e_0;
c_0 = 3e8;

alpha = -w*sqrt(mu_0*e_d)*power(1+(sigma/(w*e_d))^22, 1/4)*sin(0.5*atan(-sigma/(w*e_d)));
beta = w*sqrt(mu_0*e_d)*power(1+(sigma/(w*e_d))^2, 1/4)*cos(0.5*atan(-sigma/(w*e_d)));
lbda = 2*pi/beta;

k_c = beta - 1j*beta;

z_10db = log(10)/alpha;
disp(z_10db)

f = logspace(-1, 11);
w = 2*pi*f;
k_c = w.*sqrt(mu_0*e_d).*power(1-1j*sigma./(w*e_d), 1/2);
delta = -1./imag(k_c);
loglog(f, delta)
