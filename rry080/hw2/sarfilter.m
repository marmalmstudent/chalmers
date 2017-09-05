function output=sarfilter(H, signal)

% image = sarfilter(H, signal);
%
% Apply the SAR filter H to the radar data signal
% using Fourier transforms.
% H should be a matched filter.

% Patrik Dammert 1998-02-23
%
% Smal modifications made to avoid alliasing effects.
% Bj�rn Hallberg 2002-12-27
% get right linear convolution from ifft.
% Wiebke Aldenhoff 2016-05-24

[rows,cols]=size(H);


if (rows ~= 350) || (cols > 2000)
	error('Wrong size of filter input!');
end

image = zeros(350,4096);

% Choose between
% 1 - Fourier transform of whole matrix (uses MUCH memory)
% 2 - Fourier transform of row-by-row (uses little memory)

choice = 2;

if choice == 1
	% This implementation will be affected by alliasing 
    image = ifftrows( fftrows(H) .* fftrows(signal) );

elseif choice == 2

	h = waitbar(0,'Please wait...');
	for s=1:350,

        temp = ifft( fft(H(s,:), 8192) .* fft(signal(s,:),8192) );
        a = floor(size(H,2)/2)+1; % center portion of the convolution
        b = floor(size(H,2)/2)+4096;
        image(s,:)= temp(a:b);

        waitbar(s/350);
	end
	close(h);
end


output = image;
