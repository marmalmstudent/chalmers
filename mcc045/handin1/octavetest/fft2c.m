function out =  fft2c(x)
    out = fftshift(fft2(fftshift(x)));
end