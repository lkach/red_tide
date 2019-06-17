% [X_Coef,F,X_Coef_Unc,X_modeled,H,m,R,P] =
% red_tide(X,T,Spec,n_lowNO,df_NO,n_lowO,tide_f,sideband_centers,n_sidebands,df_sidebands,inertial)
% OR
% red_tide(X,T,{Spec,F,H,R})
% 
% NOTE: H,m,R,P are optional outputs
%
% Use "USER-FRIENDLY VERSION" if you want a single or a few time series
% analyzed, without needing to know a ton about what goes on under the
% hood. Many inputs, lots of automation.
% Use the "PRO VERSION" if you are analyzing many time series in the same
% way, and/or if you want to directly control what goes on under the hood.
% Few inputs, lots of direct control.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
%                 NECESSARY INPUTS
% 
% IN:   X = Time series (column vector, length N)
% IN:   T = Corresponding time in hours (column vector, length N)
%           This should be evenly spaced (if inclusing NaN) if "Spec" is a
%           scalar. If "Spec" is the 2-column matrix as given, then you
%           have great freedom in the kind of time series/time bases you
%           can analyze.
% 
%                 USER-FRIENDLY VERSION INPUTS (AFTER X & T)
% 
% IN:   Spec = Either a scalar or two-column matrix (size = N_spec by 2):
%           if scalar: This is the maximum number of time steps long we
%                      want our WKT power spectrum to look. If this option
%                      is used, then T must be evenly spaced.
%           if matrix: The first column is the f_spec corresponding to the
%                      power spectrum given in the second column
%           if cell:   First element must be the [f_spec, spec] like above,
%                      while the second cell element is R_frac, the
%                      residual variance fraction (see the "PRO" version
%                      inputs below).
% IN:   df_NO = The frequency step of the low non-orthogonal frequencies,
%               e.g. 1/(8760) 1/hr for df = 1/(1 non-leap yr).
% IN:   n_lowNO = The number of low non-orthogonal frequiencies modeled
% IN:   n_lowO = The number of low orthogonal frequiencies modeled, all of
%                which are higher than any low non-orthogonal frequencies.
%                Note that for this range, df = 1/(T_last - T_first)
% IN:   tide_f = Cell array of tide constituent characters formatted as
%                follows: {"S1","S2","M2",...}. Available pre-definied
%                frequencies are:
%                If you want to specify one not listed, simply write it as
%                a number instead of a string, with units 1/hr,
%                e.g. {"M2", 1/(8)} for only the M2 and third duirnal
%                harmonic (not a realistic choice for tide_f, but a good
%                example).
% IN:   sideband_centers = The tidal frequencies about which you wish to
%                          make sidebands, e.g. {"M2", 1/8}. Note that the
%                          contents of this cell must be a subset of tide_f
% IN:   n_sidebands = cell array (size = size(sideband_centers)) of
%                     whole numbers, denoting how many frequencies to look
%                     at on either side of a peak. Coordinate in array
%                     corresponds to the coordinate in "sideband_centers".
% IN:   df_sidebands = the frequency step of the sidebands frequencies
%                      about tidal peaks
% IN:   inertial = Either 0 if we don't want any inertial frequencies
%                  modeled, or a three-element vector if we do:
%                  [lat, df, N] where lat = latitude in degrees,
%                  df = frequency step for the near-inertial band (lowest
%                  frequency being f_cor), and N = the number of
%                  frequencies desired to be modeled in the NI band
% 
%                 PRO VERSION INPUT (AFTER X & T)
% 
% IN:   PRO = {Spec,F,H,R_frac}
%            Spec   = Same format as in the user-friendly version.
%            F      = Model frequencies prescribed by user (may be
%                     arbitrary, even greater than the Nyquist frequency if
%                     desired.
%            H      = OPTIONAL Model matrix (regression matrix), whose
%                     columns are sines/cosines at frequencies in F. These
%                     must correspond with those in F. This optional input
%                     is intended for repetitive cases where generating H
%                     every time from the same F would be wasteful for
%                     memory and/or CPU.
%            R_frac = OPTIONAL Residual variance fraction, i.e. prescribe
%                     this parameter instead of automatically calculate
%                     from the spectrum of X given as the prior (in
%                     "Spec"). NOTE: This overrules the optional "R" given
%                     explicitly in "Spec" if that is provided.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% OUT:  X_Coef = Matrix (M by 2) of the odd (sine, 1st column) and even
%                (cosine, 2nd column) coefficients
% OUT:  F = Vector (M by 1) of the frequencies corresponding to the
%           coefficients in X_Coef
% OUT:  X_Unc =  Matrix (M by 2) of the corresponding uncertainty (std) in
%                X_Coef
% OUT:  X_modeled = modeled X (i.e. only at frequencies in X)
% 
% OPTIONAL:
% OUT:  H = Large matrix (N by 2M) by which X_Coef (rearranged slightly)
%           can be multiplied to get X_modeled. ONLY DECLARE THIS OPTIONAL
%           OUTPUT IF YOU REALLY NEED THIS MATRIX, IT CAN BE BIG.
% OUT:  m = Same information as in X_Coef, but in the format where you
%           could calculate X_modeled = H*m as necessary.
% OUT:  R = Residual energy
% OUT:  P = Prior as used when solving the inverse problem
%
%

function [X_Coef,F,X_Coef_Unc,X_modeled,varargout] = ...
    red_tide(X,T,varargin)
% red_tide(X,T,Spec,n_lowNO,df_NO,n_lowO,tide_f,sideband_centers,n_sidebands,df_sidebands,inertial)
% OR
% red_tide(X,T,{Spec,F,H,R})

% Demean the data (and detrend if necessary):
X_nonan = X(isfinite(X));
T_nonan = T(isfinite(X));
H_detrend = [ones(length(T),1),T];
H_detrend_nonan = [ones(length(T_nonan),1),T_nonan];
m_detrend = H_detrend_nonan\X_nonan;
% X_detrend = X - H_detrend*m_detrend; % uncomment if this is necessary
% X = X_detrend; % uncomment if this is necessary
nanmeanX = nanmean(X);
X = X - nanmeanX;
% X_nonan = X(isfinite(X)); % uncomment if this is necessary
% % The line near the end where "X_modeled" is defined should also be
% % uncommented if this is necessary.

df = 1/(T(end) - T(1));
f_Ny = 1/(2*mean(diff(T)));
% ^ "T" will be evenly-spaced if the user follows the directions in the documentation

% Predefined tidal frequencies, for convenience:
f_m2 = 1/12.4206012;
f_s2 = 1/12;
f_n2 = 1/12.65834751;
f_k1 = 1/23.93447213;
f_o1 = 1/25.81933871;
f_s1 = 1/24;
f_oo1 = 1/22.30608083;
f_m4 = 1/6.210300601;
f_m6 = 1/4.140200401;
f_mk3 = 1/8.177140247;
f_s4 = 1/6;
f_mn4 = 1/6.269173724;
f_s6 = 1/4;
f_Mm = 1/661.3111655;

if nargin == 11 % no H_in or F_in or F_full
    % Build F (and then H)
    Spec = varargin{1};
    n_lowNO = varargin{2};
    df_NO = varargin{3};
    n_lowO = varargin{4};
    tide_f = varargin{5};
    sideband_centers = varargin{6};
    n_sidebands = varargin{7};
    df_sidebands = varargin{8};
    inertial = varargin{9};
    
    if size(n_sidebands) ~= size(sideband_centers)
        error('size(n_sidebands) must == size(sideband_centers)')
    else
    end
    
    Tide_Cell = {'M2','S2','N2','K1','O1','S1','OO1','M4','M6','MK3','S4',...
        'MN4','S6','Mm'};
    Full_Tide_Vec = [f_m2,f_s2,f_n2,f_k1,f_o1,f_s1,f_oo1,f_m4,f_m6,f_mk3,f_s4,f_mn4,f_s6,f_Mm];
    
    F_user_defined_tide = [];
    F_user_defined_center = [];
    if iscell(tide_f)
        F_tidal_boolean = zeros(size(Tide_Cell));
        F_cusp_boolean = zeros(size(Tide_Cell));
        for i = 1:length(tide_f)
            for j = 1:length(F_tidal_boolean)
                if ischar(tide_f{i})
                    F_tidal_boolean(j) = strcmp(Tide_Cell{j},tide_f{i}) + F_tidal_boolean(j);
                else
                    F_user_defined_tide = [F_user_defined_tide, tide_f{i}];
                end
            end
            for k = 1:length(F_cusp_boolean)
                if isempty(sideband_centers)
                    F_user_defined_center = [];
                elseif ischar(sideband_centers{i})
                    F_cusp_boolean(k) = strcmp(Tide_Cell{k},sideband_centers{i}) + F_cusp_boolean(k);
                else
                    F_user_defined_center = [F_user_defined_center, sideband_centers{i}];
                end
            end
        end
        F_tidal = Full_Tide_Vec(logical(F_tidal_boolean));
            F_tidal = [F_tidal, F_user_defined_tide];
            F_tidal = unique(F_tidal);
        F_cuspcenters = Full_Tide_Vec(logical(F_cusp_boolean));
            F_cuspcenters = [F_cuspcenters, F_user_defined_center];
            F_cuspcenters = unique(F_cuspcenters);
    else
        F_tidal = [];
        F_cuspcenters = [];
    end
    
    % Build the lowfrequency non-orthogonal and orthogonal bands:
    F_lowNO = df_NO:df_NO:(n_lowNO*df_NO);

    if n_lowNO == 0
        F_lowO = df:df:(n_lowO*df);
    else
        F_lowO = (F_lowNO(end) + df):df:(F_lowNO(end) + n_lowO*df);
    end
    
    % Build the sideband cusps
    F_cusps = [];
    for i = 1:length(F_cuspcenters)
        F_cusps = [F_cusps,...
            (F_cuspcenters(i) - n_sidebands{i}*df_sidebands):df_sidebands:(F_cuspcenters(i) - df_sidebands),...
            (F_cuspcenters(i) + df_sidebands):df_sidebands:(F_cuspcenters(i) + n_sidebands{i}*df_sidebands)];
    end
    
    if length(inertial) == 1
        F_inertial = []; % i.e. no inertial modeling
    elseif length(inertial) == 3
        % [lat, df, N]
        Om = 1/(23.9344699);
        f_cor = 2*Om*sind(inertial(1));
        df_inertial = inertial(2);
        F_inertial = f_cor:df_inertial:(inertial(3)*df_inertial);
    else
        error('The input variable "inertial" is formatted incorrectly.')
    end
    
    F = [F_lowNO, F_lowO, F_inertial, F_cusps, F_tidal];

    F = unique(F); % Removes duplicate frequencies which appear when
    % building the cusps about tidal frequencies, and also
    % sorts the vector of frequencies, just in case anything
    % is out of order (it doesn't actually matter, but it
    % does make life easier).
    
    diff_F = diff(F);
    F_unmodeled = []; % will necessarily change size in the loop
    for i=1:length(diff_F)
        if diff_F(i) > 2*df
            F_unmodeled = [F_unmodeled, (F(i) + df):df:(F(i+1) - df)];
        end
        % ^ A bit crude, but it works in one line.
    end
    
    F_full = [F, F_unmodeled, (F(end) + df):df:f_Ny]; % all frequencies
    F_full = unique(F_full); % sort and eliminate duplicates
    F_full = F_full'; % Make it into a column vector.
    F = F'; % Make it into a column vector.
    
    if sum(diff_F < 0.99*df_NO) % make sure that there are no carelessly close frequencies
        warning(['Somewhere in your F vector, you have a frequency',...
            ' difference less than the lowest you prescribed in',...
            ' df_NO. This could lead to very high uncertainty',...
            ' at this and nearby frequencies.'])
    else
    end
    if sum(diff(F_full) < 0.99*df_NO) % make sure that there are no carelessly close frequencies
        warning(['Somewhere in your F_full vector, you have a frequency',...
            ' difference less than the lowest you prescribed in',...
            ' df_NO. This could lead to very high uncertainty',...
            ' at this and nearby frequencies.'])
    else
    end
    
    % Make H from F:
    H = zeros(length(X),2*length(F));
    for i=1:length(F)
%         H(:,2*i - 1) = sin(2*pi*24*F(i)*T_nonan); % I like sine (an odd function ) having an odd  index
%         H(:,2*i) =     cos(2*pi*24*F(i)*T_nonan); % and  cosine (an even function) having an even index
        % Consider implementing H_ for getting the full X_modeled, even at
        % times when there's no data, or simply define H as the whole thing
        % and then the appropriate subset thereof for fitting
        H(:,2*i - 1) = sin(2*pi*F(i)*T); % I like sine (an odd function ) having an odd  index
        H(:,2*i) =     cos(2*pi*F(i)*T); % and  cosine (an even function) having an even index
    end
    
elseif nargin == 3
    if length(varargin{1}) == 3 || length(varargin{1}) == 4
        Spec = varargin{1}{1}; % = [f_spec, spec] or scalar (see documentation)
        F = unique(varargin{1}{2});
        H = varargin{1}{3};
        % R determined later.
    elseif length(varargin{1}) == 2
        Spec = varargin{1}{1}; % = [f_spec, spec] or scalar (see documentation)
        F = unique(varargin{1}{2});
        H = [];
        % R determined later.
    else
        error('Not enough input information given.')
    end
    F_ = F;
    if F_(1) > df; F_ = [df;F_]; else; end
    if F_(end) < f_Ny; F_ = [F_;0.5]; else; end
    diff_F = diff(F_);
    F_unmodeled = []; % will necessarily change size in the loop
    for i=1:length(diff_F)
        if diff_F(i) > 2*df
            F_unmodeled = [F_unmodeled, (F_(i) + df):df:(F_(i+1) - df)];
        end
    end
    F_full = [F; F_unmodeled'];
    F_full = unique(F_full);
else
    error('Wrong number of input arguments. See documentation.')
end

if size(H) == [length(X),2*length(F)]
else
    H = zeros(length(X),2*length(F));
    for i=1:length(F)
        H(:,2*i - 1) = sin(2*pi*F(i)*T); % I like sine (an odd function ) having an odd  index
        H(:,2*i) =     cos(2*pi*F(i)*T); % and  cosine (an even function) having an even index
    end
end

% Make spectrum, or use the one given if one is given:
N_x = length(X);
if length(Spec) == 1 && ~iscell(Spec) % scalar (max time lag, so that the spectrum need to be calculated here)
    autocov_x = zeros(1 + Spec,1); % the autocovariance
    autocov_std_x = zeros(1 + Spec,1); % the std of the points that went into each element in pre-shaped "autocov_x" (starting with zero-lag)
    n_autocov_x = zeros(1 + Spec,2); % how many pairs were finite and how many pairs could have been finite if perfect
    for i=1:(1 + Spec) % lag = (i-1)*dt
        autocov_x_i = zeros(N_x - i + 1,1);
        for j=1:(N_x - i + 1)
            autocov_x_i(j) = X(j)*X(j+i-1);
        end
        autocov_x(i) = nanmean(autocov_x_i); % autocov
        autocov_std_x(i) = nanstd(autocov_x_i); % std of points that made each element of autocov
        n_autocov_x(i,1) = sum(isfinite(autocov_x_i)); % how many pairs were finite
        n_autocov_x(i,2) = N_x - i + 1; % how many pairs could have been finite if perfect
    end
    % NOTE: "autocov_std_x" and "n_autocov_x" may be unused, and are
    % retained here for troubleshooting and possible future inclusion.
    autocov_x = [autocov_x;flip(autocov_x(2:(end-1)))]; % autocov_x = [autocov_x;0;flip(autocov_x(2:end))];
    Window = fftshift(hanning(length(autocov_x),'periodic'));
    spec = fft(autocov_x.*Window); % Window applied for smoothing
    spec = spec/length(spec); % Satisfy Parseval's Theorem
    if var(real(spec)) < 10000*var(imag(spec)) % make sure that the autocovariance is periodic and thus that its fft is real
        warning(['Imaginary part [removed] of fft(autocov_x) was ',num2str(var(imag(spec))/var(real(spec)))],' times that of the real part.')
    end
    spec = real(spec);
    if sum(spec < 0) % if the WKT spectrum is negative at any frequency
        warning(['The internally calculated spectrum (using the WKT) results',...
            ' in a partly negative spectrum. Please estimate the',...
            ' spectrum separately and use the vector option for "Spec".']);
        error('See warning.') % Comment this line if you really want to
                              % use this (nothing good will likely come of
                              % it though).
    else
    end
    % Fold the spectrum in on itself due to symmetry:
    if mod(length(spec),2) % odd length
        spec = spec(1:ceil(length(spec)/2));%%%
    elseif ~mod(length(spec),2) % even length
        spec = spec(1:ceil(length(spec)/2 + 0.5));%%%
    end
    f_Ny = 1/(2*mean(diff(T))); % "T" will be evenly-spaced if the user follows the directions in the documentation
    df_cov = f_Ny/length(spec);
    f_spec = [df_cov:df_cov:f_Ny]';
    spec(2:end) = 2*spec(2:end);
elseif size(Spec,2) == 2 && size(Spec,1) > 1 % (Spectrum and frequencies are given)
    f_spec = Spec(:,1);
    spec = Spec(:,2);
elseif length(Spec) == 1 && iscell(Spec) % same as above, this just adds more formatting options
    f_spec = Spec{1}(:,1);
    spec = Spec{1}(:,2);
elseif length(Spec) == 2 && iscell(Spec) && size(Spec{1},2) == 2 % (Spectrum and frequencies are given and the noise fraction is user-defined)
    f_spec = Spec{1}(:,1);
    spec = Spec{1}(:,2);
    R_user_defined = Spec{2};
else
    error('The 3rd input variable "Spec" is not formatted correctly. Please see the documentation.')
end

f_Ny = 1/(2*mean(diff(T))); % redone in case it's not defined (previously it only appears in an if-statements)
df = 1/(T(end) - T(1)); % Same as above
Sxx_Prior = interp1([f_spec;(f_Ny+df)],[spec;spec(end)],F_full,'linear');

% ^ Still not done: there will be NaN's in the frequencies of F outside of
% the min and max of f_vec. Make it the same as the lowest-frequency value
% of "Sxx_avg".

if isfinite(spec(1)) % make sure that you aren't replacing NaN with NaN
else
    disp(spec(1))
    error('Your averaged spectral estimates begin with a non-finite value, which absolutely should not be the case.')
end


Sxx_Prior(isnan(Sxx_Prior)) = spec(1);
% STILL not done: this Priors will need to be scaled so that they sum to
% the average variance of the power series (unlike the Fourier spectra,
% these autocovariance spectra are not automatically scaled by the
% frequency resolution, i.e. they will have a different scale depending on
% the size of max_lags, unlike Fourier spectra which are always at the same
% scale, regardless of segment count).

diffF_full = diff(F_full); % Probably better
    diffF_full = [diffF_full(1); diffF_full];
Sxx_Prior = Sxx_Prior.*(diffF_full./mean(diff(f_spec)));

% % Define the Prior covariance matrix
% We are assuming an extreme form of localization (otherwise a non-diagonal
% P would imply spurious off-diagonal correlation). The spectrum is already
% squared, so we don't have to worry about it; we just need to insert it in
% the diagonal. We only use the part of Sxx_Prior at the frequencies we are
% going to model (so really, THIS is now our true prior):

Sxx_ModelPrior = Sxx_Prior(ismember(F_full,F));

P_diag = zeros(2*length(Sxx_ModelPrior),1);
P_diag(1:2:(2*length(Sxx_ModelPrior))) = Sxx_ModelPrior;
P_diag(2:2:(2*length(Sxx_ModelPrior))) = Sxx_ModelPrior;

P = spdiags(P_diag,0,length(P_diag),length(P_diag));

% ^ P is the prior model covariance matrix. I think that the scaling is
% right, given this example:
% 
% x(t) = 5*sin(om*t) - 3*cos(om*t)
% <x^2> = 17
% autocovariance spectrum should sum to 17 then (0 everywhere except at
%       f = om/(2*pi), i.e. Sxx = [0 ... 17 ... 0]')
% true model is m = [0 0 ... 5 -3 ... 0 0]'
% diagonal of m*m' is [0 0 ... 25 9 ... 0 0]'
% Mod square of the sin and cos components (Xs and Xc) are:
%       Xs.^2 + Xc.^2 = [0...5...0]'.^2 + [0...3...0]'.^2 = [0...34...0]'
% This is the same as if Xs = Xc = \pm[0...sqrt(17)...0]', which is we need to
%       assume for a prior (that is, Xs.^2 = Xc.^2 because we have no
%       additiobnal information from Sxx).
% We only need Xs.^2 & Xc.^2, both of which we will have to assume for
%       simplicity are [0 ... 17 ... 0]'
% Thus P = a diagonal matrix with a diagonal of [0 0 ... 17 17 ... 0 0]'
% I think this is correct, I'll have to check this more rigorously.

% From Zaron (2017) JGR "MAPPING THE NONSTATIONARY INTERNAL TIDE":
% 
% "Note that the size metrics used here are computed from SSH variance, and
% the variance is equal to one-half the squared amplitude of a harmonic
% constant."

% % Obtain the data error covariance matrix
% Partition {the variance in Sxx* at F} versus {the variance in
% Sxx* NOT at F}. The second of these will be the data error
% covariance matrix's diagonal (R). Just calculate the second (We don't
% need the first, because it's in P).
% NOTE: R (data error covariance matrix) would be a constant diagonal whose
% value is estimated by integrating the spectral estimate except
% for the frequencies we intend to fit. Basically, get all the variance we
% don't expect to fit, based on our prior, and use that for R; the
% remaining variance should go into our fitted amplitudes.

R = sum(Sxx_Prior(~ismember(F_full,F))); % residual variance of X

% User-defined R:
if iscell(Spec) && length(Spec) == 2
    R = R_user_defined*sum(Sxx_Prior);
    P = P.*(sum(Sxx_Prior)./sum(Sxx_ModelPrior)).*(1 - R_user_defined);
    skip_adjustment = 1;
else
    skip_adjustment = 0;
end
% Overrule above if "R" is included in the PRO input format
if nargin == 3 && length(varargin{1}) == 4
    %   This wastes a little computation, but only the trivial amount for
    %   defining R before:
    R = varargin{1}{4}*sum(Sxx_Prior);
    P = P.*(sum(Sxx_Prior)./sum(Sxx_ModelPrior)).*(1 - varargin{1}{4});
    skip_adjustment = 1; % Unecessary, but included for clarity.
else
    % Adjust R to a value that is reasonable, if it was not user-defined:
    if skip_adjustment
    else
        if R/sum(Sxx_Prior) > 0.9
            warning(['The residual was estiated to be ',num2str(R*100/sum(Sxx_Prior)),...
                '% of the total variance expected.',...
                ' Adjusted down to ',num2str(0.9*100),...
                ' * total variance (even now the model still has limited predictive power).'])
            R = 0.9*sum(Sxx_Prior);
        elseif R <= 0
            R = 0.05*sum(Sxx_Prior);
            warning(['The residual was estiated to be 0, which leads to a singularity. Adjusted up to ',...
                num2str(R),' * total variance.'])
        else
        end
    end
end
% % % Uncomment to see what the noise-to-signal ratio is under the hood
% % % (should be approximately R/nanvar(X)).
% disp(['The residual prior after adjustment is ',num2str(R*100/sum(Sxx_Prior)),'% of the total variance']);
% disp(['The residual prior after adjustment is ',num2str( 100 * R/(R + sum(sum(P))/2) ),'% of the total variance']);

%% The fitting procedure

% % Perform the fit: (Wunsch 3.6.20)
    HRHPinv = (H(isfinite(X),:)'/R)*H(isfinite(X),:) + inv(P); % Puritan: H'*inv(R)*H + inv(P)
    m = (HRHPinv\H(isfinite(X),:)')*(R\X(isfinite(X))); % 
% % If you didn't define "H" for all times:
%     HRHPinv = (H'/R)*H + inv(P); % Puritan: H'*inv(R)*H + inv(P)
%     m = (HRHPinv\H')*(R\X_nonan);
HRHPinvinv = inv(HRHPinv);

% I originally used "HRHPinv\H'*inv(R)*TS", where "HRHPinv = H'*R*H + inv(P);",
% but I need to explicitly invert H'*R*H + inv(P) for the error estimate on
% the coeficients later, so I might as well do it now. Note: I tested "\"
% versus "inv" for simple inversion and inv was faster (note that below I'm
% using a randn matrix, whereas we will end up inverting a symmetrical one,
% so I don't know how that will matter if at all):
% 
%     A = randn(2000);
%     tic
%     Ai1 = A\speye(2000);
%     toc
%     tic
%     Ai2 = inv(A);
%     toc
%         Elapsed time is 0.435341 seconds.
%         Elapsed time is 0.340039 seconds.
% 
% 
% The puritan way would be: m = inv(H'*inv(R)*H + inv(P))*H'*inv(R)*TS;
% I tested the "Puritan" way, and it differs from the fast way by only a
% factor of 10^-13 or less (i.e. difference over either), so you can
% definitely get away with the fast way.
% 
% Maybe a more efficient way would be: (HRHPinv\H')*(R\TS)

m = m*(length(F)/(length(F) - 1))^0;
% ^ Observed to be a necessary correction

X_Coef = nan(length(F),2);
X_Coef(:,1) = m(1:2:end); % sine coefficients
X_Coef(:,2) = m(2:2:end); % cosine coefficients

m_uncertainty = sqrt(diag(HRHPinvinv));
X_Coef_Unc = nan(length(F),2);
X_Coef_Unc(:,1) = m_uncertainty(1:2:end); % sine coefficients uncertainty
X_Coef_Unc(:,2) = m_uncertainty(2:2:end); % cosine coefficients uncertainty


X_modeled = H*m;% + H_detrend*m_detrend;
                % ^ Reintroduce trend-line

if nargout == 4
elseif nargout == 5
    varargout{1} = H;
elseif nargout == 6
    varargout{1} = H;
    varargout{2} = m;
elseif nargout == 7
    varargout{1} = H;
    varargout{2} = m;
    varargout{3} = R;
elseif nargout == 8
    varargout{1} = H;
    varargout{2} = m;
    varargout{3} = R;
    varargout{4} = P;
else
end

end
