% createHfilter
%
% This is an incomplete script that creates the Hfilter.
%

Ronear=10446;       %Distance to near range [m]
delta_R=4;          %Pixel size in range direction [m]
NOfRows=350;        %Number of rows


lambda=0.0566;      %Wavelength [m]

n=-700:700;         %Update n after you determine the filter length!


PRF_v=   (Put your PRF_v here.)  ;         %PRF divided by aircraft speed [m^-1]


Hfilter=zeros(NOfRows,length(n));

for row=1:NOfRows
    R= (R should be dependent on row number.);
    Hfilter(row,:)= (Put your filter expression here.);
end;

