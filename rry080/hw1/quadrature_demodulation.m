% p=quadrature_demodulation(s,t) 
% This function performs quadrature demodulation on the signal s(t). The
% output p is the complex phasor, p=I-iQ, where I and Q are the in-phase
% and quadrature components of the signal.
% Note: It is assumed that the apperent carrier frequency is Fs/4, thus the
% program lacks generality.
function p=quadrature_demodulation(s,t)
% Find the sampling frequency
Fs = 1/(t(2)-t(1));
% Assume thet the apperent carrier frequency is Fs/4
fc_app=Fs/4;
% Find the number of samples
N = length(s);
% Define the frequency vector for zero-centered FFT
f   = (-N/2:N/2-1)*(Fs/N);
% Calculate the FFT of the I-channel, before LP-filtering
tmp1=fftshift(fft(s.*cos(2*pi*t*fc_app)));   
% Calculate the FFT of the Q-channel, before LP-filtering
tmp2=fftshift(fft(s.*sin(2*pi*t*fc_app)));   
% Do LP-filtering for I
tmp1(abs(f)>fc_app)=0;
% Do LP-filtering for Q
tmp2(abs(f)>fc_app)=0;
% Combine the I and Q channels in the frequency domain
tmp3=tmp1-1i*tmp2;
% Do an inverse fftshift and inverse FFT to get P
p = 2*ifft(ifftshift(tmp3)); % Complex phasor P

