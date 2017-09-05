clear all;
load("constants.mat")  % the constants from the table
load("sarlab.mat")

% spectrumplot(signal, size(signal, 1))
% print "task1.pdf" -dpdflatex
B = 200;
azres1 = v/B;  % best azimuth resolution
azres2 = D_a/2;
printf("Best azimuth resolution using\nBandwidth: %f m\nAntenna length: %f m\n",
	azres1, azres2);