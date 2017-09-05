clear all;
load("constants.mat")  % the constants from the table
load("sarlab.mat")

% print "task3.pdf" -dpdflatex
len_n = 1442;  % synthetic aperture in terms of pixels
len_ol = 4096;  % length of oneline
n = -len_n/2:len_n/2-1;
R0 = R0_near + 36*dR;

figure();
hold("on");
vals = cell(1,4)
PRF_v = 2.31:0.01:2.34;
for i = 1:size(PRF_v, 2)
	h = exp(2j*pi/(lbda*R0*PRF_v(i)^2) * n.^2);
	sig_fft = ifft(fft(oneline, 2*len_ol).*fft(h, 2*len_ol));
	plot(abs(sig_fft(len_n/2:len_n+len_ol/2)).^2,
		"LineWidth", 3, "Color", rand(1,3))
	vals{i} = sprintf("PRF_v: %.2f", PRF_v(i));
end
title("PRF_v")
legend(vals)