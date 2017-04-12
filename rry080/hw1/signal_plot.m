% [S,f]=signal_plot(s,t) 
% Plots the real part of the signal s as a
% function of time t, as well as the magnitude of the Fourier transform of
% s (denoted S) as a function of frequency f. S and f are also given as
% output variables. 
% The time vector t (in seconds) is assumed to be evenly spaced.
function [S,f]=signal_plot(s,t)
% Find the sampling frequency
Fs = 1/(t(2)-t(1));
% Find the number of samples
N = length(s);
% Calculate the Fourier transform (FFT)
S = fftshift(fft(s))/N;
% Define the corresponding frequency vector
f   = (-N/2:N/2-1)*(Fs/N);
% Get a new figure window
figure
% Plot the real part of the signal in the time domain
subplot(2,1,1)
plot(t*1e6,real(s))
xlabel('Time [\mus]')
ylabel('Real part of the signal')
grid
% Plot the magnitude of the Fourier transform
subplot(2,1,2)
plot(f*1e-6,abs(S)/max(abs(S)))
xlabel('Frequency [MHz]')
ylabel('Magnitude of the Fourier transform.')
grid
end
        