%Function returns saturation pressure of water (Pa)
%given temperature
function pwv_sat = antoine_eqn(T)
A = 8.07131;
B = 1730.63;
C = 233.426;
psat_mmHg = 10^(A-B/(T+C)); %mmHg
pwv_sat = psat_mmHg/760*101325; %Pa
end
