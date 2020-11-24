%boundary conditions for heat transfer
%in form of p(x,t,u)+q(x,t)f(x,t,u,DuDx) = 0
%example: -D DMDx = h(Pwv_inf - Pwv_surface) at x=0
% x=0, p = q', q=1;
function [pl,ql,pr,qr] = heat1_bc(xl,ul,xr,ur,t)
global T_inf h rho cp
pl = h*(T_inf - ul);
ql = rho*cp; %K DT/Dx, K = alpha*rho*cp
pr = -h*(T_inf - ur);
qr = rho*cp;
end