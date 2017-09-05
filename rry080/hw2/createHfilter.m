function out=createHfilter(nl, PRF_v)
% createHfilter

Ronear = 10446;       %Distance to near range [m]
delta_R = 4;          %Pixel size in range direction [m]
NOfRows = 350;        %Number of rows
lambda = 0.0566;      %Wavelength [m]
# n = -721:720;         %Update n after you determine the filter length!
n = -nl/2:nl/2-1;         %Update n after you determine the filter length!
% PRF_v = 2.33;       %PRF divided by aircraft speed [m^-1]


r = Ronear + delta_R*(0:NOfRows-1);
[N, R] = meshgrid(n, r);
out = exp(2j*pi ./ (lambda*R*PRF_v^2) .* N.^2);