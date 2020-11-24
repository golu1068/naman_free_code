%store Parabolic (one-D) partial differential equation
% c*DuDt = f*DuDx) + s
function [c,f,s]=mass1_pde(x,t,u,DuDx)
global Deff
c = 1;
f = Deff*DuDx;
s = 0;
end 