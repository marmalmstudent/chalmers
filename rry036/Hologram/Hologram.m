clc, clear

[img, map, alpha] = imread('Arch1.png');
[x, map] = rgb2ind(img);
I_2t = double(x);  # target intensity, white motive, black background

grid_size = power(2,8);  # no. samplses
holo_size = 0.02;  # [m], hologram size hsould be 2 cm
lbda = 635e-9;  # ~ red laser
k = 2*pi/lbda;  # wave number
L = 1e3;  # Distance to p2 from p1 m
sigma = 1e-2;  #std div, roughly equal to grid size

p_1_sample = holo_size/(grid_size-1);  # distance btwn each sample ping p1
x_vect = -grid_size/2*p_1_sample:p_1_sample:(grid_size/2-1)*p_1_sample;
[X, Y] = meshgrid(x_vect, x_vect);
r_mat_p_1 = sqrt(X.^2+Y.^2); # matrix containing dist to origin p1
p_2_sample = lbda*L/holo_size;  # distance between each sample point p2
u_vect = -grid_size/2*p_2_sample:p_2_sample:(grid_size/2-1)*p_2_sample;
[U, V] = meshgrid(u_vect, u_vect);
r_mat_p_2 = sqrt(U.^2+V.^2);  # matrix containing dist to origin p2
gauss_beam = exp(-r_mat_p_1.^2/sigma^2);

e_1 = double(rand(grid_size, grid_size) > 0.5);
for i=0:1
  e_2_tmp = fftshift(fft2(fftshift(e_1)));
  e_2 = sqrt(I_2t).*exp(j*angle(e_2_tmp));  # set amplitude
  e_1_tmp = fftshift(ifft2(fftshift(e_2)));
  e_1 = double(real(e_1_tmp) > 0);  # real parts>0 1, neg parts<0 0
end
DOE = e_1;
f_1 = DOE;#gauss_beam.*exp(j*k/(2*L)*(X.^2+Y.^2)).*DOE;
f_2 = fftshift(fft2(fftshift(f_1)));

black_out_radius = 8*p_2_sample;
I_2 = abs(f_2).^2;
I_2_tmp = I_2.*double(r_mat_p_2 > black_out_radius);
I_2_plot = I_2_tmp./max(max(I_2_tmp));
fig1 = figure();
imshow(I_2_plot);