% [F, Coef, Coef_Unc, Amp, y_modeled, H, x, R, P] = red_tide(t, y, FSR_cell, 'var1', val1, etc.)
% 
% Input for red_tide can get complicated and requires some level of prior
% knowledge or guesswork concerning the statistics and physics of the time
% series to be analyzed. Tidal coefficients must be chosen manually, even
% ubiquitous tides like M2, unless a preset is used. The signal-to-noise
% ratio assumed by red_tide may vary with frequency and is determined by
% the prior (assumed) autocovariances of both the model coefficients and
% the residual time series.
% 
% The equation to be solved is:         y = H*x + r
% Prior model autocovariance is:        P = <xx'>
% Prior residual autocovariance is:     R = <rr'>
% "x" is estimated by "~x":            ~x = (H' R^-1 H + P^-1)^-1 H' R^-1 y
% 
% N = number of data, M = number of frequencies
% 
% % % % % % % % % % % % % % % % % INPUTS  % % % % % % % % % % % % % % % % %
% 
% t             Nx1 vector, times corresponding to data in "y". This must
%               be in units of hours.
% 
% y             Nx1 vector, time series from which tidal information is to
%               be calculated
% 
% FSR_cell      Cell, may be either empty {} or contain three sub-cells
%               {F_cell,S_cell,R_cell}:
% 
%                   F_cell contains frequency information to build "F",
%                   "H", and "P"
%                   S_cell contains spectral information to build "P"
%                   R_cell contains information to build "R"
% 
%               Any of these except F_cell may be left empty, in which case
%               "R" and/or "P" are white noise (or are given directly, see
%               below). These are best defined in advance, as they can be
%               long and cumbersome.
% 
%               F_cell = {Preset}   (see end of documentation)
%                     OR {n_lowNO, df_NO, n_lowO, tide_f, fband_centers, n_sidebands, df_sidebands, inertial}
%               where     n_lowNO = number of non-orthogonally spaced
%                                   (low) frequencies, starting at df_NO
%                         df_NO   = arbitrary spacing of low non-orthogonal
%                                   frequencies
%                         n_lowO  = number of orthogonally-spaced (df) low
%                                   frequencies, starting just above the
%                                   highest of the non-orthogonal low
%                                   frequencies (or at df if n_lowNO = 0)
%                         tide_f  = cell with all the tidal frequencies to
%                                   be fit, formatted as characters (e.g.
%                                   'S2') or numerically (e.g. 1/12) or a
%                                   mixture, order is unimportant.
%                         fband_centers = cell, the centers of frequency
%                                   bands, may be the same as "tide_f"
%                         n_sidebands = cell (same size as "fband_centers")
%                                   with the number of frequencies away
%                                   from the corresponding center to be fit
%                         df_sidebands = scalar, the spacing of frequencies
%                                   in the bands centered around the
%                                   frequencies in "fband_centers"
%                         inertial = (Optional) three-element vector:
%                                   [lat, df, N] where lat = latitude in
%                                   degrees, df = frequency step for the
%                                   near-inertial band (lowest frequency
%                                   being f_coriolis), and N = the number
%                                   of frequencies to be modeled in the
%                                   near-inerial band.
% 
%               S_cell = {f_spec, spec}
%               where     f_spec  = frequency vector corresponding to "spec"
%                         spec    = spectrum to be used as the model prior,
%                                   where sum(spec) = variance, not
%                                   spectral power
% 
%               R_cell = {R_input, R_format, Cov_cutoff, Window}
%               where     R_input = scalar, pair, or vector:
%                                   if scalar: 0 < R_input < 1 is fraction
%                                   of variance in residual prior, with
%                                   white noise assumed
%                                   if pair: R_input = [slope,frac], where
%                                   slope <= 0 is the spectral slope of the
%                                   spectral distribution assumed for R
%                                   if vector: covariance or spectrum
%                                   corresponding to R
%                         R_format = character, one of 'c' (for
%                                   [auto]covariance) or 's' (for
%                                   'spectrum'), which tells which form
%                                   "R_input" takes. Choice only matters
%                                   for vector "R_input".
%                         Cov_cutoff = cutoff time lag beyond which zero
%                                   covariance is assumed for R.
%                         Window  = character (e.g. 'triang') for window
%                                   function used to spectrally smooth the
%                                   autocovariance to enable inversion of
%                                   R. Default is 'hanning'.
% 
% Variable/value pairs may be used after "FSR_cell" to insert 'F', 'H',
% 'P', and/or 'R' if they are already constructed. This option is available
% for extra control over what red_tide uses, as well as for computational
% efficiency when running multiple time series through red_tide using the
% same prior assumptions and modeling at the same frequencies, otherwise
% these matrices would need to be recalculated each time.
% Additionally, the option 'Fig' may be set to 'on' in order to produce
% a cursory figure that includes time- and frequency-domain information
% about the input and output, including: data (time series), model time
% series, data periodogram, model prior spectrum, residual prior spectrum,
% and model amplitude.
% 
% 
% % % % % % % % % % % % % % % % % OUTPUTS % % % % % % % % % % % % % % % % %
% 
% F             Mx1 vector of frequencies to which the time series is fit.
%               This will has units of cycles per hour.
% 
% Coef          Mx2 matrix of model coefficients, first column are sine
%               coefficients, second column are cosine coefficients. This
%               is the sine/cosine sorted version of "x", an optional
%               output (see below).
% 
% Coef_Unc      Mx2 matrix of the uncertainty of the model coefficients.
%               Formally, this is the sine/cosine sorted version of
%               "model spread about mean" <(x - ~x)(x - ~x)'>, but only the
%               diagonal elements (off-diagonal covariance is assumed to be
%               negligible).
% 
% Amp           (Optional) Mx1 vector of tidal amplitudes =
%               0.5*(Coef(:,1).^2 + Coef(:,2).^2)
% 
% y_modeled     (Optional) Nx1 modeled time series, = H*x
% 
% H             (Optional) Nx2M regressor matrix, columns are alternating
%               sines and cosine of frequencies contained in "F". This and
%               the optional outputs "R" and "P" may be reused in
%               additional applications of red_tide to prevent needless
%               recalculation.
% 
% x             (Optional) 2Mx1 model coefficients in a form that can
%               multiply by "H" to give "y_modeled"
% 
% R             (Optional) NxN prior error (residual) covariance matrix
%               (time domain)
% 
% P             (Optional) 2Mx2M prior model covariance matrix (frequency
%               domain)
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% Presets for "F_cell" may be any number from 1 to 37. There are 37 tidal
% constutents
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% [F, Coef, Coef_Unc, Amp, y_modeled, H, x, R, P] = red_tide(t, y, FSR_Cell, 'var1', val1, etc.)
%                                                                     |
%    _________________________________________________________________|
%   |
%   v
% FSR_Cell = {{n_lowNO, df_NO, n_lowO, tide_f, fband_centers, n_sidebands, df_sidebands, inertial},...
%             {f_spec, spec},...
%             {R_input, R_format, Cov_cutoff, Window}}
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [F, Coef, Coef_Unc, varargout] = red_tide(t, y, FSR_cell, varargin)

if nargin == 0
    disp('[F, Coef, Coef_Unc, Amp, y_modeled, H, x, R, P] = ')
    disp('red_tide(t, y, FSR_Cell, ''var1'', val1, etc.)')
else


% Option to detrend (1 = yes, 0 = no, only de-mean). Because detrending
% will usually be advantageous (removing broadband contamination by a
% trend), this value is adjustable only here and not as an option. Note
% that the output "y_modeled" has the trend added back in so that
% y_modeled ~= H*x in general.
detrendBoolean = 1;

if (iscell(FSR_cell) && isempty(FSR_cell)) || (iscell(FSR_cell) && (length(FSR_cell) == 3))
else
    error(['The third argument must be a {cell}, which may be empty or',...
           ' have three elements, each of which are a cell (see',...
           ' documentation for accepted formats).']);
end

% Assign variables from the contents of these cells:
if isempty(FSR_cell)
    % Rely on varargin to build F, H, P, and R
else % i.e. FSR_cell = {{...},{...},{...}}
    F_cell = FSR_cell{1};
    S_cell = FSR_cell{2};
    R_cell = FSR_cell{3};
    % Any of the above may be empty, in which case the relevant variable
    % will have to be given after FSR_Cell to avoid an error.
    
    F_vars = {'n_lowNO', 'df_NO', 'n_lowO', 'tide_f', 'fband_centers', 'n_sidebands', 'df_sidebands', 'inertial'};
    S_vars = {'f_spec', 'spec'};
    R_vars = {'R_input', 'R_format', 'Cov_cutoff', 'Window'};
    % ^ "Window" is not a required input like all the others
    
    if length(F_cell) > 1 % if length(F_cell) == 1, then preset frequencies will be used
        for i=1:length(F_cell)
            eval([F_vars{i},' = F_cell{i};'])
        end
    else
    end
    for i=1:length(S_cell)
        eval([S_vars{i},' = S_cell{i};'])
    end
    for i=1:length(R_cell)
        eval([R_vars{i},' = R_cell{i};'])
    end
end
% Note: dynamically assigning variables like above is generally bad
% practice; however, in this case the variables are constrained to the few
% named in "F_vars", "S_vars", and "R_vars", so there's no chance of this
% spiraling out of control.

if ~isempty(varargin)
    % Turn option/value pairs into a structure:
    Str = struct(varargin{:});
    Names = fieldnames(Str);
    % Verify that variables are only the ones that are allowed:
    AllowedVars = {'F','H','P','R','Fig'}';
    for i=1:length(Names)
        if ismember(Names{i},AllowedVars)
        else
            error(['''',Names{i},''' is not a possible option for red_tide. Please see documentation.'])
        end
    end
    % If that test is passed, reassign to variables without the structure
    % prefix (this should not affect memory in MATLAB because assigning two
    % variables to the same array does not copy the array until one of them is
    % changed). Nevertheless, clear "Str" after.
    for i=1:length(Names)
        eval([Names{i},' = Str.',Names{i},';'])
    end
    % Note: dynamically assigning variables like above is generally bad
    % practice; however, in this case the variables are constrained to the few
    % named in "AllowedVars", so there's no chance of this spiraling out of
    % control.
    
else
end

clear Str

%% Tidal Presets

% The ordering of the tides is based on NOAA's harmonic constituent order,
% e.g. at:
% <https://tidesandcurrents.noaa.gov/stations.html?type=Harmonic+Constituents>
% or any other site (though the amplitudes vary geographically, NOAA
% maintains a 37-element list of constituents in descending order of
% expected prominence. The order (first column of "NOAA_HARM_CONST") may be
% modified below, but this is only recommended for special cases, otherwise
% individual tidal constituents should be specified directly in "tide_f",
% an argument in "F_cell", which is in turn an argument in "FSR_cell".

NOAA_HARM_CONST = [...
    1       1/(360/28.9841042)  ;...M2
    2       1/(360/30)          ;...S2
    3       1/(360/28.4397295)  ;...N2
    4       1/(360/15.0410686)  ;...K1
    5       1/(360/57.9682084)  ;...M4
    6       1/(360/13.9430356)  ;...O1
    7       1/(360/86.9523127)  ;...M6
    8       1/(360/44.0251729)  ;...MK3
    9       1/(360/60        )  ;...S4
    10      1/(360/57.4238337)  ;...MN4
    11      1/(360/28.5125831)  ;...NU2
    12      1/(360/90        )  ;...S6
    13      1/(360/27.9682084)  ;...MU2
    14      1/(360/27.8953548)  ;...2N2
    15      1/(360/16.1391017)  ;...OO1
    16      1/(360/29.4556253)  ;...LAM2
    17      1/(360/15)          ;...S1
    18      1/(360/14.4920521)  ;...M1
    19      1/(360/15.5854433)  ;...J1
    20      1/(360/0.5443747)   ;...MM
    21      1/(360/0.0821373)   ;...SSA
    22      1/(360/0.0410686)   ;...SA
    23      1/(360/1.0158958)   ;...MSF
    24      1/(360/1.0980331)   ;...MF
    25      1/(360/13.4715145)  ;...RHO
    26      1/(360/13.3986609)  ;...Q1
    27      1/(360/29.9589333)  ;...T2
    28      1/(360/30.0410667)  ;...R2
    29      1/(360/12.8542862)  ;...2Q1
    30      1/(360/14.9589314)  ;...P1
    31      1/(360/31.0158958)  ;...2SM2
    32      1/(360/43.4761563)  ;...M3
    33      1/(360/29.5284789)  ;...L2
    34      1/(360/42.9271398)  ;...2MK3
    35      1/(360/30.0821373)  ;...K2
    36      1/(360/115.9364169) ;...M8
    37      1/(360/58.9841042) ];%..MS4
    %   order   1/(period in hr)

if ~exist('F','var') && exist('F_cell','var')
    if length(F_cell) == 1
        for i=1:F_cell{1}
            F(i) = NOAA_HARM_CONST(NOAA_HARM_CONST(:,1)==i, 2);
        end
        F = sort(F)';
    else
    end
else
end

%% Demean and detrend the data:

y_in = y;
y_nonan = y(isfinite(y));
t_nonan = t(isfinite(y));
if detrendBoolean
    H_detrend = [ones(length(t),1),t];
    H_detrend_nonan = [ones(length(t_nonan),1),t_nonan];
    x_detrend = H_detrend_nonan\y_nonan;
    y_trend = H_detrend*x_detrend;
    y = y - y_trend;
else
    y_trend = nanmean(y)*ones(size(y));
    y = y - y_trend;
end

SamplePeriod = median(diff(t_nonan)); % "mean(diff(T));" may be better.
% This choice was made to account for, say, instruments sampling only every
% other observation time while also preventing a few large gaps from making
% "SamplePeriod" too large.
Leng = (t_nonan(end) - t_nonan(1)); % In case T(1:...) and/or T(...:end) are nan
df = 1/Leng;
f_Ny = 1/SamplePeriod;

%% Sort out the flexible input for building R

if exist('R_input','var') && length(R_input) == 1  % i.e. R_input = [residual variance fraction], white noise assumption
    if R_input >= 0 && R_input <=1
        R_input = [R_input*nanvar(y); zeros(Cov_cutoff,1)];
        R_format = 'c'; % i.e. choice is an illusion in this case
    else
        error('"R_input" is formatted incorrectly.')
    end
elseif exist('R_input','var') && length(R_input) == 2
%     'R_input', 'R_format', 'Cov_cutoff'
    if R_input(1) < 0 % i.e. R_input = [negative spectral slope, residual variance fraction]
        f_R = df:df:f_Ny;
        S_R = f_R.^R_input(1);
        S_R = S_R * R_input(2)*nanvar(y)/sum(S_R);
        R_input = S_R;
        R_format = 's'; % i.e. choice is an illusion in this case
    else
        error('"R_input" is formatted incorrectly.')
    end
elseif exist('R_input','var') && isvector(R_input)
elseif ~exist('R_input','var') && exist('R_cell','var')
else
    error('"R_input" is formatted incorrectly.')
end

%% Check for prior existence of F, H, R, and P

if exist('F','var') % "F" is already defined
else
    if exist('inertial','var')
    else
        inertial = [];
    end
    F = F_make(df, n_lowNO, df_NO, n_lowO, tide_f, fband_centers, n_sidebands, df_sidebands, inertial);
end

if exist('P','var') % "P" is already defined
elseif isempty(S_cell)
    % P = P_make([df,f_Ny]',nanvar(y)*[0.5 0.5]',F,SamplePeriod,Leng);
    f_spec = [df:df:f_Ny]'; spec = nanvar(y)*ones(size(f_spec))/length(f_spec);
    P = P_make(f_spec,spec,F,SamplePeriod,Leng);
    % i.e. same model covariance at all frequencies
else
    P = P_make(f_spec,spec,F,SamplePeriod,Leng);
    % i.e. model covariance is interpolated from the given spectrum "spec"
end

if exist('H','var') % "H" is already defined
else
    H = H_make(t_nonan,F);
end

if exist('R','var') % "R" is already defined
elseif isempty(R_cell)
    R_white_frac = 0.1; % i.e. fraction of data variance expected to not be modeled at fitted frequencies
    [R,f_R,spec_R] = R_make([R_white_frac*nanvar(y); zeros(10,1)], length(y_nonan), 'c', 5, 'rectwin');
else
    if exist('Window','var')
    else
        Window = 'hanning';
    end
    % ^ "Window" is not a required input like all the others
    [R,f_R,spec_R] = R_make(R_input, length(t), R_format, Cov_cutoff, Window); % R_format = 'c' or 's'
    % Omit columns corresponding to gaps:
    R = R(isfinite(y),isfinite(y));
end

if exist('Fig','var') % "Fig" is already defined
else
    Fig = 'off';
end

%% Least Squares Problem

% This chain of try-catch's serves to whiten R by adding to the diagonal
% (zero-lag) of R when doing the Cholesky factorization until no more error
% is given by ichol (ichol only works for positive definite matrices, and
% the particular structure of R can make this condition difficult to
% achieve, especially for a noise spectrum with spectral slope <= -2).
% 
% To avoid unnecessary computation when running many time series through
% red_tide but reusing the same "R", make sure that "R" is positive
% definite.
OPTS.diagcomp = 0;
try L = ichol(R,OPTS); catch;
    warning(['"R" is not positive definite. If this calculation is to',...
    ' be repeated many times for the same "R", consider first',...
    ' verifying that "R" is positive definite to avoid unnecessary',...
    ' computation. Whiten "R" by increasing the diagonal with the',...
    ' option "diagcomp" (see MATLAB''s ichol for more info).'])
    OPTS.diagcomp = 0.01;
    disp(['"diagcomp" = ',num2str(OPTS.diagcomp)])
    try L = ichol(R,OPTS); catch;
        OPTS.diagcomp = 0.1;
        disp(['"diagcomp" = ',num2str(OPTS.diagcomp)])
        try L = ichol(R,OPTS); catch;
            OPTS.diagcomp = 1;
            disp(['"diagcomp" = ',num2str(OPTS.diagcomp)])
            try L = ichol(R,OPTS); catch;
                OPTS.diagcomp = 10;
                disp(['"diagcomp" = ',num2str(OPTS.diagcomp)])
                try L = ichol(R,OPTS); catch;
                    OPTS.diagcomp = 20;
                    disp(['"diagcomp" = ',num2str(OPTS.diagcomp)])
                    try L = ichol(R,OPTS); catch;
                        error('ichol option "diagcomp" was not enough to get L from R.')
                    end
                end
            end
        end
    end
end

yw = L\y(isfinite(y));
Hw = L\H(isfinite(y),:);
HHPinv = Hw'*Hw + inv(P);
x = (HHPinv\Hw')*yw;
HRHPinvinv = inv(HHPinv);

x = x*(length(F)/(length(F) - 1));
% ^ Observed to be a necessary correction
Coef = nan(length(F),2);
Coef(:,1) = x(1:2:end); % sine coefficients
Coef(:,2) = x(2:2:end); % cosine coefficients

x_uncertainty = sqrt(diag(HRHPinvinv));
Coef_Unc = nan(length(F),2);
Coef_Unc(:,1) = x_uncertainty(1:2:end); % sine coefficients uncertainty
Coef_Unc(:,2) = x_uncertainty(2:2:end); % cosine coefficients uncertainty


y_modeled = H*x + y_trend;
                % ^ Reintroduce trend-line
Amp = 0.5*(Coef(:,1).^2 + Coef(:,2).^2);

if strcmp(Fig,'on')
    y_0padded = y;
    y_0padded(~isfinite(y)) = nanmean(y);
    Periodogram_y = fft(y_0padded)/length(y);
    f_periodogram = [(0):(1/length(y)):(1 - 1/length(y))]';
    
    figure('Color',[1 1 1])
    subplot(2,1,1)
    plot(t,y_in-nanmean(y_in),'.-k'); hold on
    plot(t_nonan,y_modeled-nanmean(y_in),'.-r')
    plot(t,y_in-y_modeled,'b.-')
    legend(['Data, var = ',num2str(nanvar(y_in))],...
           ['Fit, var = ',num2str(nanvar(y_modeled))],...
           ['Residual, var = ',num2str(nanvar(y_in-y_modeled))])
    xlabel('Time'); ylabel('Data'); title('Data minus mean')
    subplot(2,1,2)
    loglog(24*f_periodogram(2:end),2*abs(Periodogram_y(2:end)).^2,'color',[0.5 0.5 0.5]); hold on
    loglog(24*f_spec,spec,'g.-')
    loglog(24*f_R,spec_R,'b.-')
    loglog(24*F,Amp,'r.-')
    legend(['|FFT(Data)|^2, \Sigma = ',num2str(sum(abs(Periodogram_y(2:end)).^2))],...
           ['Given spectral prior, \Sigma = ',num2str(sum(spec))],...
           ['Residual prior, \Sigma = ',num2str(sum(spec_R))],...
           ['Model coefficients squared, \Sigma = ',num2str(sum(0.5*(Coef(:,1).^2 + Coef(:,2).^2)))])
    xlabel('Frequency (cpd)'); ylabel('Variance')
    set(gcf,'Position',[440   127   560   671])
else
end

% varargout = {Amp,y_modeled,H,x,R,P}
if nargout == 3
elseif nargout == 4
    varargout{1} = Amp;
elseif nargout == 5
    varargout{1} = Amp;
    varargout{2} = y_modeled;
elseif nargout == 6
    varargout{1} = Amp;
    varargout{2} = y_modeled;
    varargout{3} = H;
elseif nargout == 7
    varargout{1} = Amp;
    varargout{2} = y_modeled;
    varargout{3} = H;
    varargout{4} = x;
elseif nargout == 8
    varargout{1} = Amp;
    varargout{2} = y_modeled;
    varargout{3} = H;
    varargout{4} = x;
    varargout{5} = R;
elseif nargout == 9
    varargout{1} = Amp;
    varargout{2} = y_modeled;
    varargout{3} = H;
    varargout{4} = x;
    varargout{5} = R;
    varargout{6} = P;
else
end

end

end
