function [u] = pde1D_convect_M(nx,hx,nt,ht,init,u_inf, K,h)
  %Based on codes in Lindfield & Penny (2012) 
  %Numerical Methods Using Matlab, Third Edition
  global Pwv_sat rho
  alpha = K*ht/hx^2;
  beta = h*hx/(K*rho); 
    A = zeros(nx-1,nx-1); u = zeros(nt+1,nx+1);
    u(:,1) = (beta*(u_inf + Pwv_sat*hheqn(u(:,2))))+u(:,2);
    u(:,nx+1) = (beta*(u_inf + Pwv_sat*hheqn(u(:,nx))))+u(:,nx);
  u(1,:) = init;
  A(1,1) = 1+2*alpha+2*alpha*beta;
  A(1,2) = -2*alpha;
  
  for i = 2:nx-2
      A(i,i) = 1+2*alpha;
      A(i,i-1) = -alpha; A(i,i+1) = -alpha;
  end
  
   A(nx-1,nx-2) = -2*alpha; A(nx-1,nx-1) = 1+2*alpha;
  b(1,1) = init(2)+init(1)*alpha;
  
  for i = 2:nx-2, b(i,1) = init(i+1); end
   b(nx-1,1) = init(nx)+init(nx+1)*2*alpha;
  [L,U] = lu(A);
  
  for j = 2:nt+1
      y = L\b; x = U\y;
      u(j,2:nx) = x'; b = x;
      b(1,1) = b(1,1)+2*alpha*beta*u_inf;
      b(nx-1,1) = b(nx-1,1)+2*alpha*beta*u_inf;
  end
