% [phaseSIN, phaseCOS] = red_tide_phase(X_Coef)
% OR
% % [phaseSIN, phaseCOS] = red_tide_phase(X_Coef,{F,Dateformat,T0,T0_new})
%
% Given the "X_Coef" output from red_tide, get the phase shift at each
% frequency. There is also the option to set the zero-point for time (e.g.
% if the phase relative to a certain time is needed).
%
% The model used in red_tide can be written as:
%
%       X = a*sin(om*t) + b*cos(om*t)
%
% This may be useful when rewritten as a phase-shifted sin or cosine:
%
%       X = A*sin(om*t + phaseSIN)
%       X = A*cos(om*t + phaseCOS)
%
% where:
%
%       A = sqrt(a^2 + b^2)
%       phaseSIN = atan2( b,a)
%       phaseCOS = atan2(-a,b)
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% This assumes that "X_Coef" is formatted exactly as it comes out of
% red_tide, i.e. Mx2 in size, where M = the number of modeled frequencies,
% with the first column being sine coefficients and the second column being
% cosine coefficients.
%
% IN:   X_Coef   = Mx2 matrix, output from red_tide
% IN:   OptionIn = Optional cell input for the case where the zero-time for
%                  phase calculation needs to be specified (as opposed to
%                  simply using the time corresponding with T==0 for the
%                  input "T" in red_tide). The cell's elements are:
%
%                  F = Mx1 column of frequencies (not angular frequencies)
%                      used in red_tide (as output and sometimes input).
%                      Therefore, this must have units of hr^-1.
%                  Dateformat = format that the internal call of MATLAB's
%                      datenum will use, e.g. 'yyyy-mm-dd HH:MM:SS' (default)
%                      Enter the empty string '' for the default; if this
%                      is used, then the following arguments must be
%                      formatted appropriately.
%                  T0 = The time corresponding to zero for the red_tide
%                       input "T". Format according to Dateformat.
%                  T0_new = The shifted time (e.g. 3 years earlier than T0)
%                           at which phases are to be calculated. This too
%                           must be formatted according to Dateformat.
%
% OUT:  phaseSIN = Mx1 matrix, phase shift of a sine wave at corresponding
%                  frequencies to those of X_Coef (i.e. "F" from red_tide)
% OUT:  phaseCOS = Mx1 matrix, phase shift of a cosine wave at corresponding
%                  frequencies to those of X_Coef (i.e. "F" from red_tide)
%

% Note: if this is to be translated into a different scientific programming
% language, be aware that the implementation of the arctangent may be
% different. For example, Mathematica's ArcTan[a,b] gives the same answer
% as MATLAB's atan2(b,a)

function [phaseSIN, phaseCOS] = red_tide_phase(X_Coef,varargin)
if nargin == 1
    phaseSIN = atan2( X_Coef(:,2),X_Coef(:,1));
    phaseCOS = atan2(-X_Coef(:,1),X_Coef(:,2));
elseif nargin == 2
    F = varargin{1}{1};
    Dateformat = varargin{1}{2};
    T0 = varargin{1}{3};
    T0_new = varargin{1}{4};
    
    % Numerical form (in "days from January 0, 0000")
    t0 = datenum(T0,Dateformat);
    t0_new = datenum(T0_new,Dateformat);
    dt_hours = (t0 - t0_new)*24;
    
    phaseSIN = mod(atan2( X_Coef(:,2),X_Coef(:,1)) + dt_hours*2*pi*F, 2*pi);
    phaseCOS = mod(atan2(-X_Coef(:,1),X_Coef(:,2)) + dt_hours*2*pi*F, 2*pi);
else
    error('Incorrect number of inputs, please see documentation.')
end
end

%% arctangent demonstration, to aid in translation
% close all
% a=-10:0.1:10;
% b=-10:0.1:10;
% [A,B]=meshgrid(a,b);
% C=atan2(A,B);
% figure;surf(A,B,C);colorbar('southoutside');colormap hsv;
% xlabel('a');ylabel('b');shading interp
