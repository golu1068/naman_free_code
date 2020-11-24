%boundary conditions
%in form of p(x,t,u)+q(x,t)f(x,t,u,DuDx) = 0
%example: -D DMDx = h(Pwv_inf - Pwv_surface) at x=0
% x=0, p = q', q=1;
function [pl,ql,pr,qr] = mass1_bc(xl,ul,xr,ur,t)
global Pwv_inf Pwv_sat RH
global hm rho
%RH = hheqn(ul)
%RH = hheqn(ur)
pl = -hm*(Pwv_inf - RH*Pwv_sat);
ql = rho; %ql*D*DMDx
pr = hm*(Pwv_inf - RH*Pwv_sat);
qr = rho;
end