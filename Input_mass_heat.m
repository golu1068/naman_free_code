clc
close all
clear
format long
%Input file for heat and mass transfer
global rho K cp h Deff hm
global T_inf Pwv_inf Pwv_sat RH


L = 2e-2 %m thickness of wood
tend = 3600*360; %s length of exposure time
node_x = 20;
node_t = 50;

m=0; %0 is cartesian, 1 is cylindrical, 2 is spherical
%Heat Transfer - Inputs 
K = 0.2 %W/m.K
rho = 600 %kg/m^3
cp = 2500 %J/kg-K
h = 5 %convection heat transfer coefficient of air W/m^2.K
T_inf = 12+273%K

%Mass transfer- Inputs
Pa = 101325; %air Pressure Pa
Cpa = 1010; %Heat capacity of air J/kg.K
Ka = 0.0285; %Thermal conductivity of air W/m.K
Ra = 287.058; %Gas Constant J/kg.K
T = T_inf; %K
Deff = 1e-10; % moisture diffusivity, [m^2/s]
hm = 1E-10;
0.622*h*Deff^(2/3)/(Pa*Cpa^(1/3)*(Ka*Ra*T)^(2/3)) %mass transfer coefficient [s/m]

Pwv_inf = 1170; %pa
Pwv_sat = antoine_eqn(T-273);
RH = 0.6;

%Solving mass transfer PDE
x = linspace(0,L, node_x);
t = linspace(0,tend,node_t);
sol = pdepe(m, @mass1_pde, @mass1_ic,@mass1_bc,x,t);
Moisture = sol(:,:,1);

%Plot results of mass PDE
figure, plot(x, Moisture(1,:),x,Moisture(11,:),x,Moisture(21,:),x,Moisture(31,:),x,Moisture(41,:),x,Moisture(end,:)) %plot temp profile at all coordinate at end time
xlabel('x (m)')
  ylabel('Moisture (kg moisture/kg solid)')
  legend('t= 0 day','t=3 day','t=6 day','t=9 day','t=12 day','t=15 day')
  %legend('t= 0 h','t=2 h','t=4 h','t=6 h','t=8 h','t=10 h')
figure, plot(t,Moisture(:,1),t,Moisture(:,5),t,Moisture(:,9),t,Moisture(:,13),t,Moisture(:,17),t,Moisture(:,end)) %plot tem profile at all time at x = 1;
xlabel('t (s)')
  ylabel('Moisture (kg moisture/kg solid)')
  legend('x = 0.000 m','x=0.004 m','x=0.008 m','x=0.012 m','x=0.016 m','x=0.020 m')
%figure, surfl(Moisture)
%xlabel('x - node nos.')
  %ylabel('Time - node nos.'), zlabel('Moisture (kg moisture/kg solid)')
  

 %Solving heat transfer PDE
x = linspace(0,L,node_x);
t = linspace(0,tend,node_t);
sol = pdepe(m, @heat1_pde, @heat1_ic,@heat1_bc,x,t);
Temperature = sol(:,:,1);

%Plot results of heat PDE
figure, plot(x,Temperature(1,:),x,Temperature(5,:),x,Temperature(11,:),x,Temperature(21,:),x,Temperature(31,:),x,Temperature(41,:),x,Temperature(end,:)) %plot temp profile at all coordinate at end time
xlabel('x (m)')
  ylabel('Temperature (K)')
  legend('t= 0 h','t=1 h','t=2 h','t=4 h','t=6 h','t=8 h','t=10 h')
figure, plot(t,Temperature(:,1),t,Temperature(:,5),t,Temperature(:,9),t,Temperature(:,13),t,Temperature(:,17),t,Temperature(:,end)) %plot tem profile at all time at x = 1;
xlabel('t (s)')
  ylabel('Temperature (K)')
  legend('x = 0.000 m','x=0.004 m','x=0.008 m','x=0.012 m','x=0.016 m','x=0.020 m')
%figure, surfl(Temperature)
%xlabel('x - node nos.')
  %ylabel('Time - node nos.'), zlabel('Temperature (K)')
  
 