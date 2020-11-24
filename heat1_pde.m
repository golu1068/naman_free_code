%store Parabolic (one-D) partial differential equation
% c*DuDt = f*DuDx + s
function [c,f,s]=heat1_pde(x,t,u,DuDx)
global rho cp K
c = 1;
f = K/(rho*cp)*DuDx;
s = 0;
end 