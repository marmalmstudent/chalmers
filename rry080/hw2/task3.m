clear all;
load("constants.mat")  % the constants from the table
load("sarlab.mat")

hold("on");
% plot(abs(oneline).^2)
% print "task3.pdf" -dpdflatex
len_n = 1442;  % synthetic aperture in terms of pixels
len_ol = 4096;  % length of oneline
n = -len_n/2:len_n/2-1;
R0 = R0_near + 36*dR;
h = exp(2j*pi/(lbda*R0*PRF_v^2) * n.^2);

% FFT
printf("FFT:\n");
subplot(121)
tic
sig_fft = ifft(fft(oneline, 2*len_ol).*fft(h, 2*len_ol));
toc
plot(abs(sig_fft(len_n/2:len_n+len_ol/2)).^2)
title("FFT");

% Convolve
printf("Convolve:\n");
subplot(122)
tic
sig_conv = conv(oneline,h);
toc
plot(abs(sig_conv(len_n/2:len_n+len_ol/2)).^2)
title("Convolve");