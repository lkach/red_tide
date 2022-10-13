% [R, freq, S_R] = R_make(IN, N, Format, Cov_cutoff, Window)
% 
% Accessory function to red_tide. This builds the spectrum corresponding to
% the R matrix (which is in the time domain). This functionality exists in
% R_make.m, but for cases where that functionality was not called but the
% spectrum corresponding to R is needed after the fact, this function
% does that.
% 
% INPUT:    R = The R-matrix that is given by R_make or red_tide (whoch
%                calls R_make)
% 
% OUTPUT:   freq = frequency vector corresponding to S_R below.
%                Note that here, df = 1/(2N) [because double the
%                autocovariance is used to build the spectrum] and
%                freq(end) = 0.5. That is, the units of "freq" are not
%                1/hours but rather 1/(dt/hour).
%                E.g. if dt = 1/3 hours, then multiply "freq" by 3/hour to
%                put it in units of 1/hour.
% 
%           S_R = spectrum of a process with the covariance of the first
%                column of R. This is what red_tide actually "sees" for
%                model covariance in the spectral domain when it's
%                estimating model coefficients.

function [freq, S_R] = Spectrum_from_R(R)

if size(R,1) == size(R,2)
else
    error('R must be square')
end

N = size(R,1);

freq = [(1/((N-1))):(1/((N-1))):0.5]';
R_col = R(1:round(N),1);
R_col = full(R_col);
S_R = ifft([R_col; flip(R_col(2:(end-1)))]);
if rms(real(S_R))/rms(imag(S_R)) < 10^6
    disp(['rms real = ',num2str(rms(real(S_R))),',   rms imag = ',num2str(rms(imag(S_R)))]);
else
end
S_R = real(S_R);
S_R = 2*S_R(2:(length(freq)+1));

end