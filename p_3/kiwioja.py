import numpy as np


def Np(phi):
    a = 6378137 # wielka polos elipsoidy
    b = 6356752.3142 # mala polos elipsoidy
    e2 =  0.00669438002290 # kwadrat pierwszego mimośrodu elipsoidy
    N = a / np.sqrt(1-(e2*np.sin(phi)**2)) # promien krzywizny N
    return N

def Mp(phi):
    a = 6378137 # wielka polos elipsoidy
    b = 6356752.3142 # mala polos elipsoidy
    e2 =  0.00669438002290 # kwadrat pierwszego mimośrodu elipsoidy
    M = a * (1 - e2) / (1 - e2 * np.sin(phi)**2)**(3/2) # promien krzywizny M
    return M

phi = np.deg2rad(52)
lam = np.deg2rad(21)

s = 41234
az = np.deg2rad(45)


def kivioj(phi, lam, s, az):
    # podzielenie linii geodezyjnej na n elementow ds
    n = round(s/1000)
    ds = s/n
    
    for i in range(n):
        # 1. Obliczamy główne promienie krzywizny N i M w punkcie wyjściowym P1 oraz stałą c linii geodezyjnej.
        N_i = Np(phi) # promien krzywizny N
        M_i = Mp(phi) # promien krzywizny M

        # 3. Pierwsze przybliżenie przyrostu szerokości i azymutu:
        dphi_i = ds*np.cos(az)/M_i
        dA_i = ds*np.sin(az)/np.tan(phi)/N_i

        # 4. Obliczenie szerokości i azymutu w punkcie środkowym (m) odcinka, na podstawie przyrostów:
        phi_im = phi + dphi_i/2
        az_im = az + dA_i/2

        # 5. Obliczenie promieni krzywizn w kierunkach głównych w punkcie m:
        N_im = Np(phi_im)
        M_im = Mp(phi_im)

        # 6. Ostateczne przyrosty szerokości, długości i azymutu:
        dphi_i = ds*np.cos(az_im)/M_im
        dA_i = ds*np.sin(az_im)/(np.tan(phi_im)/N_im)
        dlam_i = ds*np.sin(az_im)/(N_im/np.cos(phi_im))

        phi = phi + dphi_i
        az = az + dA_i
        lam = lam + dlam_i
        az_odw = az + np.pi

    return phi, lam, az_odw

from pyproj import geod

g = geod.Geod(ellps='WGS84')

phi_k, lam_k, az_odw_k = kivioj(phi, lam, s, az)

phi_p, lam_p, az_odw_p = g.fwd(lam, phi, np.rad2deg(az), s)

import sys
sys.path.insert(1, 'C:/Users/wojci/Documents/Python Scripts/Geodezja')
from funkcje import deg2dms, rad2dms, dms2rad

print(f'phi Kivioj: {rad2dms(phi_k)}')