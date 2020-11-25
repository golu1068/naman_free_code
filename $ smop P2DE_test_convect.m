 clc
 clear
global Pwv_sat RH rho
  
thick = 0.018;%Thickness in metre 
tfinal = 3600; %total time of 22,000 seconds
nx = 10; hx = thick/nx; %number of spatial step and spatial step size
nt = 10; ht = tfinal/nt; %number of temporal step and temporal step size
init_T = 100*ones(1,nx+1); 
init_M = 0.25*ones(1,nx+1);
%Heat Transfer - Inputs 
k = 0.2; %W/m.K
rho = 600; %kg/m^3
cp = 2500; %J/kg-K
K = k/(rho*cp);
h = 5; %convection heat transfer coefficient of air W/m^2.K
T_inf = 12+273;%K

%Mass transfer- Inputs
Pa = 101325; %air Pressure Pa
Cpa = 1010; %Heat capacity of air J/kg.K
Ka = 0.0285; %Thermal conductivity of air W/m.K
Ra = 287.058; %Gas Constant J/kg.K
T = T_inf; %K
Deff = 1e-10; % moisture diffusivity, [m^2/s]
hm = 0.622*h*Deff^(2/3)/(Pa*Cpa^(1/3)*(Ka*Ra*T)^(2/3)); %mass transfer coefficient [s/m]
Pwv_inf = 1170; %pa
Pwv_sat = antoine_eqn(T-273);

  [Temperature] = pde1D_convect(nx,hx,nt,ht,init_T,T_inf,K,h);
  
  [Moisture] = pde1D_convect_M(nx,hx,nt,ht,init_M,Pwv_inf,Deff,hm);
  
  