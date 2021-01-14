% P = P_make(f,S,F,SamplePeriod,L)
% 
% Accessory function to red_tide. This builds the P matrix, which is
% especially useful for when the internal function that does this is not
% wanted. This assumes no off-diagonal model covariance (i.e. independent
% behavior at different frequencies). If you want to include off-diagonal
% variance in your model, you can add that to this matrix.
% 
% INPUT:    f = (vector) the frequency basis corresponding with S
% 
%           S = (vector) the spectrum to be used as a prior for model
%               coefficient covariance in red_tide. This should NOT be in
%               units of spectral power density, but rather variance, i.e.
%               sum(S) = var(time series), or possibly ~= if you want to
%               use a prior spectral distribution different from that of
%               the time series.
% 
%           F = (vector) the frequencies at which red_tide is to model
%               (i.e. the ones from which the sinusoids in "H" are made).
%               Same units as that of f.
% 
%           SamplePeriod = the period of sampling the time series in the
%               inverse units of f and F. E.g. sampling once per 2 hours
%               corresponds to SamplePeriod = 2
% 
%           L = time series length (also in inverse units of f and F). This
%               is especially important to give explicitly when "F" does
%               not have a frequency step the same size as 1/record length.
% 
%           INTERP_METHOD = string, may be any of the interpolation methods
%               available to "interp1", e.g. 'linear' or 'spline'. The
%               default option is 'loglinear', which is not an option for
%               "interp1" but rather an option that tells interp1 to
%               linearly interpolate the log of "S", and then run it
%               through exp(). This reduces the effect of a single sharp
%               peak in "S" giving too much energy a priori to neighboring
%               frequencies. It also prevents the already smeared spectrum
%               from being even more biased in favor of overestimating
%               amplitudes near what are actually sharp peaks.

function P = P_make(f,S,F,varargin)

if nargin == 3
    SamplePeriod = 1;
    L = 1/min(diff(F));
    INTERP_METHOD = 'loglinear';
elseif nargin == 5
    SamplePeriod = varargin{1};
    L = varargin{2};
    INTERP_METHOD = 'loglinear';
elseif nargin == 6
    SamplePeriod = varargin{1};
    L = varargin{2};
    INTERP_METHOD = varargin{3};
else
    error(['3 or 5 or 6 inputs expected, given ',num2str(nargin),'. See documentation.'])
end

f_Ny = 1/(2*SamplePeriod);
df = 1/L;

F_ = F;
if F_(1) > df; F_ = [df;F_]; else; end
if F_(end) < f_Ny; F_ = [F_;f_Ny]; else; end
diff_F = diff(F_);
F_unmodeled = []; % will necessarily change size in the loop
for i=1:length(diff_F)
    if diff_F(i) > 2*df
        % i.e. no frequencies in "F" should be spaced greater than or equal
        % to 2*df, otherwise this normalization won't work. Luckily, this
        % will not be the case in any properly-posed red_tide-applicable
        % problem. If, however, such an arrangement is truly resired,
        % simply make "P" with this function and manually adjust the values
        % up.
        F_unmodeled = [F_unmodeled, (F_(i) + df):df:(F_(i+1) - df)];
    end
end
F_full = [F; F_unmodeled'];
F_full = unique(F_full);


if strcmp(INTERP_METHOD,'loglinear')
    Sxx_Prior = exp(interp1([f;(f_Ny+df)],log([S;S(end)]),F_full,'linear'));
else
    Sxx_Prior = interp1([f;(f_Ny+df)],[S;S(end)],F_full,INTERP_METHOD);
end
Sxx_Prior(Sxx_Prior<=0) = min(Sxx_Prior(Sxx_Prior>0));

% ^ Still not done: there will be NaN's in the frequencies of F outside of
% the min and max of f_vec. Make it the same as the lowest-frequency value
% of "Sxx_avg".

if isfinite(S(1)) % make sure that you aren't replacing NaN with NaN
else
    disp(S(1))
    error('Your averaged spectral estimates begin with a non-finite value, which absolutely should not be the case.')
end


Sxx_Prior(isnan(Sxx_Prior)) = S(1);
% STILL not done: this Priors will need to be scaled so that they sum to
% the average variance of the power series (unlike the Fourier spectra,
% these autocovariance spectra are not automatically scaled by the
% frequency resolution, i.e. they will have a different scale depending on
% the size of max_lags, unlike Fourier spectra which are always at the same
% scale, regardless of segment count).

diffF_full = diff(F_full); % Probably better
    diffF_full = [diffF_full(1); diffF_full];
Sxx_Prior = Sxx_Prior.*(diffF_full./mean(diff(f)));

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

end