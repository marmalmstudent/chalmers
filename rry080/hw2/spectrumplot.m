%
%	M-file for easy access to estimates of 
%	a signal's spectrum.
%	
%	spectrumplot(signal,nlines)
%
%	Signal is the raw data and nlines is the number 
%	of lines you want in the estimate of the
%	spectrum. 

%	Version 1.0 Patrik Dammert 950523
%   Version 1.1 Anders Berg    100827 clean-up



function spectrumplot(signal,nlines)

powof2 = 256;
signal_reshaped = reshape(signal(1:nlines,:)',powof2,[])';

%figure;
plot(linspace(-155,155,powof2),mean(abs(fftshift(fft(signal_reshaped,[],2))).^2));
grid;
xlabel('Frequency (Hz)')
ylabel('Power')
title('Power Spectral Density of azimuth spectrum')


