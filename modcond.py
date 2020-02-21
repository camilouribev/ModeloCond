## Importar librerías 

import matplotlib.pyplot as plt
import math
import numpy as np
import scipy
import scipy.integrate as integrate
from functools import partial



def R_dif_Wrikamanayake (T_amb, T_int, t_c, h_fg):
#----------------Resistencia a la difusión del vapor en la capa de no condensables [K/W]----------------
    R_dif=(T_amb-T_int)/(t_c*h_fg)
    return R_dif

def R_cond_aire_Wrikamanayake (r,k_aire, delta):
#--------------Resistencia a la conducción del aire [K/W]---------------------
    R_cond_aire=delta/(k_aire*math.pi*(r**2))
    return R_cond_aire
#--Definición de h_i (coeficiente interfacial)
#h_i=Coeficiente interfacial de transferencia de Calor

h_i=((2*alfa)/(2-alfa))*(1/math.sqrt(2*math.pi*R_g*T_sat))*((h_fg**2)/(v_g*T_s))# Calcular de otra manera

def R_int_liq_Wrikamanayake (r, theta, alfa, T_s):


#--------Resistencia de Interfase líquido-gas [K/W]----------------------------

    R_int=1/(h_i*2*math.pi*(r**2)*(1-math.cos(theta)))

    return R_int

def R_cond_gotas_Wrikamanayake(r,theta, k1):
    """"Resistencia a la conducción de las gotas [K/W]"""

    R_drop=theta/(4*math.pi*k1*r*math.sin(theta))

    return R_drop

def Q_total_n(r,R_dif,T_amb,T_s):
 "Resistencia total por gotas pequeñas"""
 Q_total_n=((T_amb-T_s)/(((1/R_dif+(1/R_cond_dif(r)))**(-1))+R_int_liq(r)+R_got(r)))*n(r)
 return Q_total_n
    
def Q_total_N(r,R_dif,T_amb,T_s,R_cond_dif,R_int_liq,R_got):
#Resistencia totala por gotas grandes'''
 Q_total_N = ((T_amb-T_s)/(((1/R_dif+(1/R_cond_dif(r)))**(-1))+R_int_liq(r)+R_got(r)))*N(r)
 return Q_total_N


def N_LeFevre(r, r_max):
 """ Distribución de tamaño de gotas para gotas grandes """
 N = 1/(3*math.pi*r**2*r_max) * (r/r_max)**(-2/3)
 return N

    
def r_e_KimKim(N_s):
 """ Radio efectivo de la gota """
 r_e = (4*N_s)**(-0.5)
 return r_e

r_e=r_e_KimKim(N_s)
    
def n_KimKim(r, deltaT_sub, r_min, delta_coat, k_coat, k_c, theta, h_i, rho, h_fg, r_e, r_max):
 """ Distribución de tamaño de gotas para gotas pequeñas """
 A_1 = deltaT_sub / (2*rho*h_fg)
 A_2 = theta * (1-math.cos(theta)) / (4*k_c*math.sin(theta))
 A_3 = 1/(2*h_i) + delta_coat*(1-math.cos(theta)) / (k_coat*(math.sin(theta))**2)
 tau = 3*r_e**2 * (A_2*r_e + A_3)**2 / (A_1*(11*A_2*r_e**2 - 14*A_2*r_e*r_min + 8*A_3*r_e - 11*A_3*r_min))
 B_1 = A_2/(tau*A_1) * ((r_e**2-r**2)/2 + r_min*(r_e-r) - r_min**2*math.log((r-r_min)/(r_e-r_min))) 
 B_2 = A_3/(tau*A_1) * (r_e-r - r_min*math.log((r-r_min)/(r_e-r_min)))
 n = 1/(3*math.pi*r_e**3*r_max) * (r_e/r_max)**(-2/3) * r*(r_e-r_min)/(r-r_min) * (A_2*r+A_3)/(A_2*r_e+A_3) * math.exp(B_1+B_2)
 return n 


#-----------RESISTENCIA TÉRMICA TOTAL [K/W]---------------

R_dif=R_dif_Wrikamanayake(T_amb, T_int ,t_c, h_fg)
R_cond_aire=R_cond_aire_Wrikamanayake(r,k_aire, delta)
R_int= R_int_liq_Wrikamanayake (r, theta, alfa, T_s)
R_drop=R_cond_gotas_Wrikamanayake(r, theta, k1)



deltaT_sub=T_amb-T_s
# Definir funciones para la tasa de transferencia de calor



R_cond_dif = partial(R_cond_aire_Wrikamanayake, k_aire=k_aire, delta=delta)
R_cond_dif.__doc__ = "Resistencia para conducción del aire  en cada gota dependiente del radio r en m"  
   
R_int_liq = partial(R_int_liq_Wrikamanayake, theta=theta, alfa=alfa, T_s=T_s)
R_int_liq.__doc__ = "Resistencia para interfase líquida dependiente del radio en m" 
    
    
R_got= partial(R_cond_gotas_Wrikamanayake,theta=theta, k1=k1)
R_got.__doc__ = "Resistencia conducción en las gotas, dependiente del r en m" 
    
    

n = partial(n_KimKim, deltaT_sub=deltaT_sub, r_min=r_min, delta_coat=delta_coat, k_coat=k_coat, k_c=k_c,
                theta=theta, h_i=h_i, rho=rho, h_fg=h_fg, r_e=r_e, r_max=r_max)
n.__doc__ = "drop size distribution for small drops depending on drop radius r in m"
N = partial(N_LeFevre, r_max=r_max)
N.__doc__ = "drop size distribution for large drops depending on drop radius r in m"

    # integrate and calculate heat flux density
    

"""Integración de los calores de acuerdo al número de gotas"""

    
q_n, q_n_interr = integrate.quad(Q_total_n, r_min, r_max,(R_dif, T_amb, T_s))
q_N, q_N_interr = integrate.quad(Q_total_N, r_min, r_max,( R_dif, T_amb, T_s, R_cond_dif, R_int_liq, R_got))
    
    
q = q_n + q_N

    
print("R_dif:", R_dif,"\n"
"R_cond_aire:",R_cond_aire,"\n"
      "R_int:", R_int,"\n"
      "R_drop:",R_drop,"\n"
      "R_dif:",R_dif,"\n"
      "T_amb:",T_amb,"\n"
      "q_N:",q_N,"W/m2\n"
      "q_n",q_n,"W/m2\n"
      "q:", round(q,2),"W/m2\n"
        "Relacion q_N/q:",q_N/q,"\n")

print ("h_i:", round(h_i,1), "\n")


