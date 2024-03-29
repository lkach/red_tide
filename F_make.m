% F = F_make(N) where N is an integer from 1 to 37, indicating the first N
% tidal constituents to fit (ordere by NOAA, see F_make.m).
% 
% OR
% 
% F = F_make(df, n_lowNO, df_NO, n_lowO, tide_f, sideband_centers, n_sidebands, df_sidebands, inertial)
%
% "inertial" may be [] if it's not desired


%            F_make(df,n_lowNO,df_NO,n_lowO,tide_f,sideband_centers,...
%                     n_sidebands,df_sidebands,inertial)
function F = F_make(varargin)

% Predefined tidal frequencies, expand as necessary:
f_m2 =        1/(360/28.9841042 ); % M2
f_s2 =        1/(360/30         ); % S2
f_n2 =        1/(360/28.4397295 ); % N2
f_k1 =        1/(360/15.0410686 ); % K1
f_m4 =        1/(360/57.9682084 ); % M4
f_o1 =        1/(360/13.9430356 ); % O1
f_m6 =        1/(360/86.9523127 ); % M6
f_mk3 =       1/(360/44.0251729 ); % MK3
f_s4 =        1/(360/60         ); % S4
f_mn4 =       1/(360/57.4238337 ); % MN4
f_nu2 =       1/(360/28.5125831 ); % NU2
f_s6 =        1/(360/90         ); % S6
f_mu2 =       1/(360/27.9682084 ); % MU2
f_2n2 =       1/(360/27.8953548 ); % 2N2
f_oo1 =       1/(360/16.1391017 ); % OO1
f_lam2 =      1/(360/29.4556253 ); % LAM2
f_s1 =        1/(360/15         ); % S1
f_m1 =        1/(360/14.4920521 ); % M1
f_j1 =        1/(360/15.5854433 ); % J1
f_mm =        1/(360/0.5443747  ); % MM
f_ssa =       1/(360/0.0821373  ); % SSA
f_sa =        1/(360/0.0410686  ); % SA
f_msf =       1/(360/1.0158958  ); % MSF
f_mf =        1/(360/1.0980331  ); % MF
f_rho =       1/(360/13.4715145 ); % RHO
f_q1 =        1/(360/13.3986609 ); % Q1
f_t2 =        1/(360/29.9589333 ); % T2
f_r2 =        1/(360/30.0410667 ); % R2
f_2q1 =       1/(360/12.8542862 ); % 2Q1
f_p1 =        1/(360/14.9589314 ); % P1
f_2sm2 =      1/(360/31.0158958 ); % 2SM2
f_m3 =        1/(360/43.4761563 ); % M3
f_l2 =        1/(360/29.5284789 ); % L2
f_2mk3 =      1/(360/42.9271398 ); % 2MK3
f_k2 =        1/(360/30.0821373 ); % K2
f_m8 =        1/(360/115.9364169); % M8
f_ms4 =       1/(360/58.9841042 ); % MS4

if nargin == 9
    % The more elaborate, customizable input format.
    df = varargin{1};
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
    
    Tide_Cell =     {'M2','S2','N2','K1','M4','O1','M6','MK3','S4','MN4',...
        'Nu2','S6','MU2','2N2','OO1','Lam2','S1','M1','J1','Mm',...
        'Ssa','Sa','Msf','Mf','Rho','Q1','T2','R2','2Q1','P1',...
        '2SM2','M3','L2','2MK3','K2','M8','MS4'};
    Full_Tide_Vec = [f_m2, f_s2, f_n2, f_k1, f_m4, f_o1, f_m6, f_mk3, f_s4, f_mn4,...
        f_nu2, f_s6, f_mu2, f_2n2, f_oo1, f_lam2, f_s1, f_m1, f_j1, f_mm,...
        f_ssa, f_sa, f_msf, f_mf, f_rho, f_q1, f_t2, f_r2, f_2q1, f_p1,...
        f_2sm2, f_m3, f_l2, f_2mk3, f_k2, f_m8, f_ms4];
    
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
        end
        for i = 1:length(sideband_centers)
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
        F_cuspcenters = [F_cuspcenters, unique(F_user_defined_center)];
    else
        F_tidal = [];
        F_cuspcenters = [];
    end
    
    % Build the low-frequency non-orthogonal and orthogonal bands:
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
    
    if isempty(inertial)
        F_inertial = []; % i.e. no inertial modeling
    elseif length(inertial) == 3
        % [lat, df, N]
        Om = 1/(23.9344699);
        f_cor = 2*Om*sind(inertial(1));
        df_inertial = inertial(2);
        F_inertial = f_cor:df_inertial:(f_cor + inertial(3)*df_inertial);
    else
        error('The input variable "inertial" is formatted incorrectly.')
    end
    
    F = [F_lowNO, F_lowO, F_inertial, F_cusps, F_cuspcenters, F_tidal]';
    
    F = unique(F); % Removes duplicate frequencies which appear when
    % building the cusps about tidal frequencies, and also
    % sorts the vector of frequencies, just in case anything
    % is out of order (it doesn't actually matter, but it
    % does make life easier).
    
    % Make sure that there are no carelessly close frequencies
    if df_sidebands==0 || ~isfinite(df_sidebands) % in case df_sidebands was not meant to be used when building F
        df_sidebands = df;
    else
    end
    diff_F = diff(F);
    if isempty(F_lowNO)
        if sum(diff_F < 0.99*df_sidebands) % if the gap is less than df_sidebands (because df_lowNO is not given)
            warning(['In your F vector, you have a frequency',...
                ' difference less than 1/(',num2str(1/df_sidebands),...
                ' hr) = df_sidebands (or the fundamental frequency',...
                ' df = 1/(record length) if df_sidebands was given',...
                ' as 0, i.e. skipped). This could lead to very high',...
                ' uncertainty at this and nearby frequencies.',...
                ' The corresponding F indices are:'])
            F_Indices = 1:length(F);
            disp(F_Indices(diff_F < 0.99*df_sidebands))
        else
        end
    else
        if sum(diff_F < 0.99*df_NO) % if the gap is less than df_lowNO
            warning(['In your F vector, you have a frequency',...
                ' difference less than the lowest you prescribed in',...
                ' df_NO = 1/(',num2str(1/df_NO),' hr). This could lead',...
                ' to very high uncertainty at this and nearby',...
                ' frequencies. The corresponding F indices are:'])
            F_Indices = 1:length(F);
            disp(F_Indices(diff_F < 0.99*df_NO))
        else
        end
    end
    
    
elseif nargin == 1
    
    % The ordering of the tides is based on NOAA's harmonic constituent
    % order, e.g. at:
    % <https://tidesandcurrents.noaa.gov/stations.html?type=Harmonic+Constituents>
    % or any other site (though the amplitudes vary geographically, NOAA
    % maintains a 37-element list of constituents in descending order of
    % expected prominence. The order (first column of "NOAA_HARM_CONST")
    % may be modified below, but this is only recommended for special
    % cases, otherwise individual tidal constituents should be specified
    % directly in "tide_f", an argument in "F_cell", which is in turn an
    % argument in "FSR_cell".
    
    NOAA_HARM_CONST = [...
    1       f_m2  ;...M2
    2       f_s2  ;...S2
    3       f_n2  ;...N2
    4       f_k1  ;...K1
    5       f_m4  ;...M4
    6       f_o1  ;...O1
    7       f_m6  ;...M6
    8       f_mk3 ;...MK3
    9       f_s4  ;...S4
    10      f_mn4 ;...MN4
    11      f_nu2 ;...NU2
    12      f_s6  ;...S6
    13      f_mu2 ;...MU2
    14      f_2n2 ;...2N2
    15      f_oo1 ;...OO1
    16      f_lam2;...LAM2
    17      f_s1  ;...S1
    18      f_m1  ;...M1
    19      f_j1  ;...J1
    20      f_mm  ;...MM
    21      f_ssa ;...SSA
    22      f_sa  ;...SA
    23      f_msf ;...MSF
    24      f_mf  ;...MF
    25      f_rho ;...RHO
    26      f_q1  ;...Q1
    27      f_t2  ;...T2
    28      f_r2  ;...R2
    29      f_2q1 ;...2Q1
    30      f_p1  ;...P1
    31      f_2sm2;...2SM2
    32      f_m3  ;...M3
    33      f_l2  ;...L2
    34      f_2mk3;...2MK3
    35      f_k2  ;...K2
    36      f_m8  ;...M8
    37      f_ms4];%..MS4
    %   order   1/(period in hr)
    
    F = nan(varargin{1},1);
    for i=1:varargin{1}
        F(i) = NOAA_HARM_CONST(NOAA_HARM_CONST(:,1)==i, 2);
    end
    
    F = unique(F);
    
else
    error('Wrong number of inputs.')
end



end
