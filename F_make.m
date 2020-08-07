% F = F_make(df, n_lowNO, df_NO, n_lowO, tide_f, sideband_centers, n_sidebands, df_sidebands, inertial)
% 
% "inertial" may be [] if it's not desired


function F = F_make(df,n_lowNO,df_NO,n_lowO,tide_f,sideband_centers,...
                    n_sidebands,df_sidebands,inertial)

% Predefined tidal frequencies, expand as necessary:
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

if size(n_sidebands) ~= size(sideband_centers)
    error('size(n_sidebands) must == size(sideband_centers)')
else
end

Tide_Cell =     { 'M2','S2','N2','K1','O1','S1','OO1','M4','M6','MK3','S4','MN4','S6','Mm'};
Full_Tide_Vec = [f_m2,f_s2,f_n2,f_k1,f_o1,f_s1,f_oo1,f_m4,f_m6,f_mk3,f_s4,f_mn4,f_s6,f_Mm ];
% ^ Make sure that these two variables match.

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

F = [F_lowNO, F_lowO, F_inertial, F_cusps, F_tidal]';

F = unique(F); % Removes duplicate frequencies which appear when
% building the cusps about tidal frequencies, and also
% sorts the vector of frequencies, just in case anything
% is out of order (it doesn't actually matter, but it
% does make life easier).

end