% function [si,ti]=interpolate_signal(s,t,L)
% This function interpolates the signal s(t) by the interpolation factor L,
% using zero-padding in the frequency domain. The time vector t must be
% evenly spaced. 
% The output si is the interpolated signal at the imterpolation times ti. 
function [si,ti]=interpolate_signal(s,t,L)
% Find the length of the signal
N=length(s);
% Find the time spacing
dt = t(2) - t(1);
% Calculate the length of the interpolated signal
M=N*L;
% Calculate the time vector for the interpolated signal
ti = min(t)+(0:M-1)*dt/L;
% Calculate the interpolated signal
si = interpft(s,M);