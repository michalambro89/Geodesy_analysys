import numpy as np

a = 6378137
f = 1/298.257223563
b = a*(1-f)
e2 = (a**2-b**2)/a**2
e2p = (a**2-b**2)/b**2

def Np(B: 'float', a : 'int' = 6378137, e2 : 'float' = 0.00669438002290) -> 'float':
    # TODO: description
    N = a/(1- e2*(np.sin(B)**2))**0.5
    return N

def hirv(x,y,z):
    p = np.sqrt(x**2+y**2)
    
    phi = np.arctan(z/(p*(1-e2)))
    
    while True:
        N = Np(phi)
        h = p/np.cos(phi) - N
        phi_poprzednie = phi
        phi = np.arctan((z/p) * (1-(N*e2)/(N+h))**(-1))
        if abs(phi_poprzednie-phi)<(0.000001/60/60/360):
            break
    lam = np.arctan2(y,x)
    return phi, lam, h