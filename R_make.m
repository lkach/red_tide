% [R, freq, S_R] = R_make(IN, N, Format, Cov_cutoff, Window)
% 
% Accessory function to red_tide. This builds the R matrix.
% 
% INPUT:    IN = vector, either a spectrum if "Format" is one of 's',
%                'spec', or 'spectrum', or an (auto)covariance if "Format"
%                is one of 'c', 'cov', 'covariance'.
%                - If "IN" is a spectrum, it is assumed that IN(end)
%                corresponds to the Nyquist frequency and IN(1) corresponds
%                to the fundamental frequency, with spacing = df = 1/T,
%                both in units of 1/hour, and where T > Cov_cutoff.
%                - If "IN" is an autocovariance, it is assumed that IN(1)
%                corresponds to zero lag (i.e. variance), and that each
%                successive element corresponds to one additional time lag
%                equal to the time step of the time series being analyzed.
% 
%           N = the dimension of R, whch will be NxN
% 
%           Format = string, see above for accepted values in "IN".
% 
%           Cov_cutoff = (Optional) the index of the element of the
%                covariance (either given or internally calculated from
%                spectrum) after which zero is assumed. Default is
%                length(IN).
% 
%           WindowF = (Optional) Windowing function by which the covariance
%                is multiplied (this can significantly reduce spectral
%                ringing while maintaining its spectral shape). Default is
%                'rectwin' (i.e. no windowing). The optional variable
%                "Cov_cutoff" must also be included if "WindowF" is used.
% 
% 
% OUTPUT:   R = The R-matrix (error covariance) to be used in red_tide
% 
%           freq = (OPTIONAL) frequency vector corresponding to S_R below.
%                Note that here, df = 1/(2N) [because double the
%                autocovariance is used to build the spectrum] and
%                freq(end) = 0.5. That is, the units of "freq" are not
%                1/hours but rather 1/(dt/hour).
%                E.g. if dt = 1/3 hours, then multiply "freq" by 3/hour to
%                put it in units of 1/hour.
% 
%           S_R = (OPTIONAL) spectrum of a process with the covariance of
%                a column of R. This is what red_tide actually "sees" for
%                model covariance in the spectral domain when it's
%                estimating model coefficients.

function [R, varargout] = R_make(IN,N,Format,varargin)

if nargin == 3
    Cov_cutoff = length(IN);
    WindowF = 'rectwin';
elseif nargin == 4
    Cov_cutoff = varargin{1};
    WindowF = 'rectwin';
elseif nargin == 5
    Cov_cutoff = varargin{1};
    WindowF = varargin{2};
else
    error('Wrong number of inputs.')
end

Window = eval([WindowF,'(Cov_cutoff*2 - 1)']);
Window = Window((Cov_cutoff):end);
Window = Window';

if iscolumn(IN); IN = IN'; else; end

if strcmp(Format,'s') || strcmp(Format,'spec') || strcmp(Format,'spectrum')
    S = [0,IN]; % 0 for zero mean, because fft treats the first element as zero frequency
    SS = [S flip(S(2:(end-1)))];
    fSS = ifft(SS)/length(SS);
    if rms(real(fSS)) < 10000*rms(imag(fSS))
        error('Non-trivial imaginary component to the Fourier transform of the given spectrum.')
    else
        fSS = real(fSS);
    end
    C = fSS(1:Cov_cutoff);
    C = sum(S)*C/C(1);
    C = C.*Window;
    R = spdiags(repmat([flip(C(2:end)),C],N,1), (1-length(C)):(length(C)-1), speye(N));
    % figure;plot(0:(length(C)-1),C,'.-')
elseif strcmp(Format,'c') || strcmp(Format,'cov') || strcmp(Format,'covariance')
    C = IN(1:Cov_cutoff);
    C = C.*Window;
    R = spdiags(repmat([flip(C(2:end)),C],N,1), (1-length(C)):(length(C)-1), speye(N));
    % figure;plot(0:(length(C)-1),C,'.-')
else
    error('Invalid second input "Format".')
end

if nargout == 1
elseif nargout == 3
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
    
    varargout{1} = freq;
    varargout{2} = S_R;
else
    warning('Either 1 or 3 outputs expected (see documentation).')
end

end
