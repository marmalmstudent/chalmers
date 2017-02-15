N = 256;  # number fo sample points
lbda_0 = 635e-9;  # wavelength in vacuum
n = 1.33;  # refractive index
k = 2*pi*n/lbda_0;
a_anim = 100e-7;  # sample distance at animation plane 1
b_anim = 99e-7;  # sample distance at animation plane 2
E1 = 0;  # the starting field for the animation
img = 0;  # the animation image
doe = 0;  # difractive optical element
lense = 0;
f_lense = 0;  # focal length of lense


[xmat, ymat] = meshgrid(getCoordVect(N, a_anim), getCoordVect(N, a_anim));
rmat = sqrt(xmat.^2 + ymat.^2);
maskRadius = 400e-6;
bool_mat = (rmat > maskRadius);
E1 = ones(N, N).*bool_mat;
L = 1e-2;
[E2, L1, L2, uvect, vvect, umat, vmat] = HFM(E1, a_anim, b_anim, L, lbda_0, n);
E2_plt = abs(E2)./max(max(abs(E2)));
imshow(E2_plt)
# [xmat, ymat] = np.meshgrid(getCoordVect(self.N, self.b_anim), getCoordVect(self.N, self.b_anim))
# rmat = np.sqrt(xmat**2 + ymat**2)
# maskRadius = 70e-5
# E2 = E2*(rmat > maskRadius)

# E2 = E2*(np.abs(E2)/np.max(np.abs(E2)) > 0.2)
