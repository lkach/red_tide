% red_tide_plot(F,X_Coef)
% OR
% red_tide_plot(F,X_Coef,xScale)
% 
% Takes the red_tide outputs "F" and "X_Coef" and makes a new figure that
% displays the energy and phase of the harmonic decomposition calculated by
% red_tide.
% 
% F      is Mx1, units of hr^-1
% X_Coef is Mx2, units of "X" (input for red_tide)
% xScale is optional, either 'log' (default) or 'linear', and sets the
%        x-axis scaling (y-scale is always log for the top plot, linear for
%        the bottom one).
% 

function red_tide_plot(F,X_Coef,varargin)

[phaseS,phaseC] = red_tide_phase(X_Coef);

figure
subplot(211)
loglog(F,0.5*(X_Coef(:,1).^2 + X_Coef(:,2).^2),'.:')
xlabel('Frequency (cph)','interpreter','latex')
ylabel({'Variance contribution';'(units of time series squared)'},'interpreter','latex')
if nargin == 3
    set(gca,'xscale',varargin{1})
end
subplot(212)
semilogx(F,phaseS,'.:');hold on
semilogx(F,phaseC,'.:')
legend('\phi_s_i_n','\phi_c_o_s')
xlabel('Frequency (cph)','interpreter','latex')
ylabel('Phase w.r.t. $t = 0$','interpreter','latex')
set(gca,'ylim',[-pi,pi])
if nargin == 3
    set(gca,'xscale',varargin{1})
end

end
