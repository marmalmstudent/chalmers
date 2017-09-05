clear all;
load("constants.mat")  % the constants from the table
load("sarlab.mat")

hold("on");
len_n = 1442;  % synthetic aperture in terms of pixels
len_ol = 4096;  % length of oneline
n = -len_n/2:len_n/2-1;
PRF_v = 2.33;

r = R0_near + dR*(0:size(signal, 1)-1);
[N, R] = meshgrid(n, r);
h = exp(2j*pi ./ (lbda*R*PRF_v^2) .* N.^2);

% h = createHfilter(len_n, PRF_v);
image = sarfilter(h, signal);
imageplot(image)

ns = 9;
Ns = ones(1,ns)/ns;
s = ceil(ns/2); e = (s-1) + size(image, 2);
nlook_image = sqrt(conv2(abs(image).^2, Ns)(:, s:e));
figure()
imageplot(nlook_image)

figure()
plot(abs(nlook_image(37, :)).^2)
