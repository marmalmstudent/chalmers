% u_filt=matched_filter(u,h)
% Calculates the filtered signal u_filt. u is the signal to be filtered,
% and h is the (matched) filter.
% Hint: 
% - The time reverse xrev=x(-t) of a vector x(t) can be obtained by the
% command xrev = x(end:-1:1); 
function u_filt=matched_filter(u,h)
% Convolve with h(t)
u_filt=conv(u,h); 
% Remove the edge-parts of the filtered signal
u_filt = u_filt(floor(length(h)/2)+(1:length(u))) ; 
