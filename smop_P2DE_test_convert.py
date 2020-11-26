import numpy as np
import scipy.linalg
###################################################################
def antoine_eqn(T):
    global A, B, C
    A = 8.07131;
    B = 1730.63;
    C = 233.426;
    psat_mmHg = 10**(A-B/(T+C)); #mmHg
    pwv_sat = psat_mmHg/760*101325; #Pa
    return pwv_sat

def hheqn(u):
    np.seterr(divide='ignore', invalid='ignore')
    A = 1.2733890;
    B = 0.121262;
    C = 0.001009;
    n = len(u)
    phi = np.zeros(np.shape(u))
    for i in range(0,n):
        phi[i] = ((B*u[i]-1)+np.sqrt(((1-B*u[i]))**2+4*A*C*((u[i]))**2))/(2*C*u[i])
        
    return phi


def pde1D_convect(nx,hx,nt,ht,init,u_inf, K,h):
    alpha = K*ht/(hx**2)
    beta = h*hx/K
    A = np.zeros((nx-1, nx-1))
    u = np.zeros((nt+1,nx+1));
    u[0,:] = init
    A[0][0] = 1+2*alpha+2*alpha*beta;
    A[0][1] = -2*alpha;
    
    for i in range(1, 8):
        A[i][i] = 1+2*alpha;
        A[i][i-1] = -alpha;
        A[i][i+1] = -alpha;
    A[nx-2][nx-3] = -2*alpha; 
    A[nx-2][nx-2] = 1+2*alpha;
    
    b = np.zeros((nx-1, nx-1))
    b[0][0] = init[1]+init[0]*alpha;
    
    for i in range(1, 8):
        b[i][0] = init[i+1]
    b[nx-2][0] = init[nx-1] + init[nx-1]*2*alpha
    
    P, L, U = scipy.linalg.lu(A)
    L[[7,8]] = L[[8,7]]


    for j in range(1, 11):
        y,resid,rank,s = np.linalg.lstsq(L,b,rcond=None)
        x,resid,rank,s = np.linalg.lstsq(U,y, rcond=None)
        u[j][1:nx] = np.transpose(x)[0]
        b = np.array(x)
        b[0] = b[0] + 2*alpha*beta*u_inf;
        b[nx-2] = b[nx-2] + 2*alpha*beta*u_inf;
    
    u[:,0] = (beta*u_inf + u[:,1])/(1+beta);
    u[:,nx] = (beta*u_inf + u[:,nx-1])/(1+beta);

    return u

def pde1D_convect_M(nx,hx,nt,ht,init,u_inf, K,h):
    alpha = K*ht/(hx**2)
    beta = h*hx/(K*rho)
    A = np.zeros((nx-1, nx-1))
    u = np.zeros((nt+1,nx+1));
    
    u[0,:] = init
    A[0][0] = 1+2*alpha+2*alpha*beta;
    A[0][1] = -2*alpha;
    
    for i in range(1, 8):
        A[i][i] = 1+2*alpha;
        A[i][i-1] = -alpha; 
        A[i][i+1] = -alpha;
        
    A[nx-2][nx-3] = -2*alpha; 
    A[nx-2][nx-2] = 1+2*alpha;
    b = np.zeros((nx-1, nx-1))
    b[0][0] = init[1]+init[0]*alpha;
    
    for i in range(1, 8):
        b[i][0] = init[i+1]
    b[nx-2][0] = init[nx-1] + init[nx-1]*2*alpha
    
    P, L, U = scipy.linalg.lu(A)
    
    for j in range(1, 11):
        y,resid,rank,s = np.linalg.lstsq(L,b,rcond=None)
        x,resid,rank,s = np.linalg.lstsq(U,y, rcond=None)
        u[j][1:nx] = np.transpose(x)[0]
        b = np.array(x)
        b[0] = b[0] + 2*alpha*beta*u_inf;
        b[nx-2] = b[nx-2] + 2*alpha*beta*u_inf;
    u[:,0] = (beta*(u_inf + Pwv_sat*hheqn(u[:,1])))+u[:,1];
    u[:,nx] = (beta*(u_inf + Pwv_sat*hheqn(u[:,nx-1])))+u[:,nx-1];

    return u
#################################################################
class main():
    def __init__(self, **kwargs):
        super(main, self).__init__(**kwargs)
        global rho, Pwv_sat
        thick = 0.018;##Thickness in metre 
        tfinal = 3600; ##total time of 22,000 seconds
        nx = 10; 
        hx = thick/nx; ##number of spatial step and spatial step size
        nt = 10; 
        ht = tfinal/nt; ##number of temporal step and temporal step size
        init_T = 100*np.ones(nx+1); 
        init_M = 0.25*np.ones(nx+1);
        ##Heat Transfer - Inputs 
        k = 0.2; ##W/m.K
        rho = 600; ##kg/m^3
        cp = 2500; ##J/kg-K
        K = k/(rho*cp);
        h = 5; ##convection heat transfer coefficient of air W/m^2.K
        T_inf = 12+273;##K
        
        ##Mass transfer- Inputs
        Pa = 101325; ##air Pressure Pa
        Cpa = 1010; ##Heat capacity of air J/kg.K
        Ka = 0.0285; ##Thermal conductivity of air W/m.K
        Ra = 287.058; ##Gas Constant J/kg.K
        T = T_inf; ##K
        Deff = 1e-10; ## moisture diffusivity, [m^2/s]
        hm = 0.622*h*Deff**(2/3)/(Pa*Cpa**(1/3)*(Ka*Ra*T)**(2/3)); ##mass transfer coefficient [s/m]
        Pwv_inf = 1170; ##pa
        
        Pwv_sat = antoine_eqn(T-273);
        
        Temperature = pde1D_convect(nx,hx,nt,ht,init_T,T_inf, K,h)
        
        Moisture = pde1D_convect_M(nx,hx,nt,ht,init_M,Pwv_inf,Deff,hm);
        
        ##########################################################################
        x = np.linspace(0,L,nx+1);
        
        Moisture = np.array(Moisture)

        fig1 = plt.figure(1)
        plt.plot(x, Moisture[0])
        plt.plot(x, Moisture[2])
        plt.plot(x, Moisture[4])
        plt.plot(x, Moisture[6])
        plt.plot(x, Moisture[8])
        plt.plot(x, Moisture[10])
        plt.legend(['t= 0 h', 't= 2 h' , 't=4 h', 't=6 h', 't=8 h', 't=10 h'])
        plt.xlabel('x (m)')
        plt.ylabel('Moisture (kg moisture/kg solid)')
        plt.xlim(0,0.02)

if __name__ == "__main__":
    main()
