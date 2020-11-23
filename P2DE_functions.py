import numpy as np
import scipy.linalg

def pde1D(nx, hx, nt, ht, init, lowb, hib, K):
    ## Solves parabolic equ'n.
      ## e.g. heat flow equation.
      ## Example call: [u alpha] = heat(nx,hx,nt,ht,init,lowb,hib,K)
      ## nx, hx are number and size of x panels
      ## nt, ht are number and size of t panels
      ## init is a row vector of nx+1 initial values of the function.
      ## lowb & hib are boundaries at low and hi values of x.
      ## Note that lowb and hib are scalar values.
      ## K is a constant in the parabolic equation.
    alpha = K*ht/(hx**2);
    A = np.zeros((nx-1, nx-1))
    u = np.zeros((nx+1, nx+1))
    u[:,0] = lowb##*np.ones((nt+1,1))
    u[:,nx] = hib##ones((nt+1, 1))
    u[0,:]= init[0]
    A[0][0] = 1+2*alpha
    A[0][1] = -alpha

    for i in range(1, nx-2):
        A[i][i] = 1+2*alpha
        A[i][i-1] = -alpha
        A[i][i+1] = -alpha
    

    A[nx-2][nx-3] = -alpha
    A[nx-2][nx-2] = 1+2*alpha
    
##    print(A)
    b=[];
    b.append(init[1]+init[0]*alpha)
    for i in range(1, 8):
        b.append(init[i])

    b.append(init[nx-1]+init[nx]*alpha)
##    print(b)

    P, L, U = scipy.linalg.lu(A)
##    print(L)
##    print(b)
    
    for j in range(1, nt+1):
        y,resid,rank,s = np.linalg.lstsq(L,b)
        x,resid,rank,s = np.linalg.lstsq(U,y)
        u[j][1:nx] = np.transpose(x)
        b = np.array(x)
        b[0] = b[0] + lowb*alpha
        b[nx-2] = b[nx-2] + hib*alpha
    print('\n')
    print('u', u)

    return u, alpha

    
######################################################################        
K = 5e-7;##thermal diffussivity of brick [5 x 10^-7 m/s2]
thick = 0.3;##Thickness in metre 
tfinal = 3600; ##total time of 22,000 seconds
nx = 10;
hx = thick/nx; ##number of spatial step and spatial step size
nt = 10;
ht = tfinal/nt; ##number of temporal step and temporal step size
init = 100*np.ones(nx+1);
lowb = 20;
hib = 20; ##initial value (t = 0), BC at x=0, x = 0.3 m
