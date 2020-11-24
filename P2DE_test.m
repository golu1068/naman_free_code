%P2DE_test
%Test m file for parabolic second order partial equation
 % e3s604.m
 % K*D2u/dx2 = du/dt
 
 %Problem
 %We now use the function heat to study how the temperature distribution 
 %in a brick wall varies with time. The 
 %wall is 0.3 m thick and is initially at a uniform temperature of 100�C. 
 %For the brickwork, K = 5 � 10?7 m/s2. If the temperature of 
 %both surfaces is suddenly lowered to 20� C and kept at this temperature, we wish 
 %to plot the subsequent variation of temperature through the wall 
 %at 440 s (7.33 min) intervals for 22,000 s (366.67 min).
 clc
 clear
 
  K = 5e-7;%thermal diffussivity of brick [5 x 10^-7 m/s2]
  thick = 0.3;%Thickness in metre 
  tfinal = 3600; %total time of 22,000 seconds
  
  nx = 10; hx = thick/nx; %number of spatial step and spatial step size
  nt = 10; ht = tfinal/nt; %number of temporal step and temporal step size
  init = 100*ones(1,nx+1); lowb = 20; hib = 20; %initial value (t = 0), BC at x=0, x = 0.3 m
  
  %Evaluate temperature, u at 15 spatial nodes and 50 temporal nodes
  %And evaluate thermal diffissivity, al
  [u al] = pde1D(nx,hx,nt,ht,init,lowb,hib,K);
  alpha = al, 
  % Generate three-dimensional shaded surfaces 
  surfl(u)
  axis([0 nx+1 0 nt+1 0 120])
  view([-217 30]), xlabel('x - node nos.')
  ylabel('Time - node nos.'), zlabel('Temperature')