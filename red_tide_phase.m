% [phaseSIN, phaseCOS] = red_tide_phase(X_Coef)
% 
% Given the "X_Coef" output from red_tide, get the phase shift at each
% frequency.
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

function [phaseSIN, phaseCOS] = red_tide_phase(X_Coef)
phaseSIN = atan2( X_Coef(:,2),X_Coef(:,1));
phaseCOS = atan2(-X_Coef(:,1),X_Coef(:,2));
end

%% arctangent demonstration, to aid in translation
% close all
% a=-10:0.1:10;
% b=-10:0.1:10;
% [A,B]=meshgrid(a,b);
% C=atan2(A,B);
% figure;surf(A,B,C);colorbar('southoutside');colormap hsv;
% xlabel('a');ylabel('b');shading interp
