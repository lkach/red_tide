% QUANTITIES = red_tide_derived_quantities(Coef)
% 
% Given the "Coef" output(s) from red_tide, get all derived values mentioned
% in the accompanying PDF (e.g. phase, rotary quantities, etc.). Output
% depends on in the input is real (one quantity) or complex (u and v,
% rotary)
% 
% IN:
%   Coef = Mx2 coefficients (1st column sine coef, second column cosine
%          coef). If you want derived quantities for a scalar observed
%          quantity (e.g. pressure) then "Coef" is real. If you want rotary
%          quantities for a vector quantity (e.g. u and v) then "Coef" is
%          complex and defined as follows: Coef = u_Coef + 1i*v_Coef
% 
% OUT:
%   QUANTITIES = Output structure containing the following quantities if
%          "Coef" is real:
%            - phaseSIN (for adjustments to starting time, use "red_tide_phase")
%            - phaseCOS
%            - var
%          If "Coef" is complex:
%            - u_phaseSIN, v_phaseSIN
%            - u_phaseCOS, v_phaseCOS
%            - u_var, v_var
%            - Rcw, Rccw
%            - theta
%            - rotary_phase
%            - SM_axis, sm_axis (Semi-Major and semi-minor axes)
%            - eccentricity
% 

function QUANTITIES = red_tide_derived_quantities(Coef)
if isreal(Coef)
    aa = Coef(:,1);
    bb = Coef(:,2);
    
    QUANTITIES.phaseSIN = atan2( bb,aa);
    QUANTITIES.phaseCOS = atan2(-aa,bb);
    QUANTITIES.var = 0.5*[aa.^2 + bb.^2];
else
    aa = real(Coef(:,1));
    bb = real(Coef(:,2));
    cc = imag(Coef(:,1));
    dd = imag(Coef(:,2));

    QUANTITIES.u_phaseSIN = atan2( bb,aa);
    QUANTITIES.u_phaseCOS = atan2(-aa,bb);
    QUANTITIES.u_var = 0.5*[aa.^2 + bb.^2];

    QUANTITIES.v_phaseSIN = atan2( dd,cc);
    QUANTITIES.v_phaseCOS = atan2(-cc,dd);
    QUANTITIES.v_var = 0.5*[cc.^2 + dd.^2];

    QUANTITIES.Rcw  = 0.5*[bb - cc + 1i*[dd + aa]];
    QUANTITIES.Rccw = 0.5*[bb + cc + 1i*[dd - aa]];

    QUANTITIES.theta = 0.5*[atan2( dd-aa , bb+cc ) + atan2( dd+aa , bb-cc )];
    QUANTITIES.rotary_phase = ...
                       0.5*[atan2( dd-aa , bb+cc ) - atan2( dd+aa , bb-cc )];
    QUANTITIES.SM_axis =     abs(QUANTITIES.Rccw) + abs(QUANTITIES.Rcw);
    QUANTITIES.sm_axis = abs(abs(QUANTITIES.Rccw) - abs(QUANTITIES.Rcw));
    QUANTITIES.eccentricity = sqrt(1 - [QUANTITIES.sm_axis ./ QUANTITIES.SM_axis]);
end
end
