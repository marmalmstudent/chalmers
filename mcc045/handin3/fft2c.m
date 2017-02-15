function out=fft2c(x)
    % shifted fourier transform
    out = fftshift(fft2(fftshift(x)));
end