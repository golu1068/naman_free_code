%Function calculates relative humidity, phi from 
%moisture content, u
function phi = hheqn(u)
A = 1.2733890;
B = 0.121262;
C = 0.001009;
phi = ((B*u-1)+sqrt((1-B*u)^2+4*A*C*u^2))/(2*C*u);
     
end
