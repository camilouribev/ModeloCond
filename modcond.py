
import math
import numpy as np
import scipy
import scipy.integrate as integrate


#--------------------Nomenclatura---------------------------------
#R_dif=resistencia térmica por la capa de difusión

#T_amb=Temperatura ambiente [K]

#T_int=Temperatura de la gota [K]

#t_c=Tasa de condensación. Hallada experimentalmente

#h_fg= Calor latente de condensación [J/kg]

# K_aire=Conductividad térmica del K_aire [W/m K]

#delta= grosor de la capa límite [m]

#S_r=Área proyectada de la gota [m²]

#R_g: Constante específica de los gases J/[Kg K]

#v_g: Volumen especifico del vapor de agua [m²/Kg]
#h_fg: Calor latente de vaporización [J/Kg]

#alfa=Coeficiente de acomodación o condensación;
# fracción de moléculas de vapor que se mueven en la gota líquida durante el cambio de fase
#debe ser de entre 0.02 a 0.04 para agua. En Wrikamanayake et al. se usó 0.02

#T_s=Temperatura de superficie de pared [K]

#h_i=Coeficiente interfacial de transferencia de Calor

#r=radio de la gota [m]

#theta=Ángulo de contacto [°]

T_amb=
T_int=
t_c=
h_fg=
K_aire=
delta=
S_r=
R_g=
v_g=
h_fg=
alfa=
T_s=


#Método de resistencias
#Se irán haciendo las resistencias desde el exterior de la gota hasta el sustrato
def R_dif_Wrikamaniyake (T_amb, T_int, t_c, h_fg):
#----------------Resistencia a la difusión del vapor en la capa de no condensables [K/W]----------------
    R_dif=(T_amb-T_int)/(t_c*h_fg)
    return R_dif

def R_cond_aire_Wrikamanayake (K_aire, delta, S_r):
#--------------Resistencia a la conducción del aire [K/W]---------------------
    R_cond_aire=delta/(k_aire*S_r)
    return R_cond_aire

#------------------SECCIÓN PENDIENTE------------------------------

# def Capa_limite (h_fg, R, T_int, T_boil, RH,P_sat, P_amb, D, M, c, x_i, x_a, Delta):
# #------Inicialmente, Ecuación Clausius-Clapeyron para determinar x_i, termino necesario------------------
#     #T_boil=Punto de fusión a presión atmosférica (P_atm)
#     #R=
#
#     x_i=math.exp(-(h_fg/R)*((1/T_int)-(1/T_boil)))
#     #x_a=Fracción molar en el agua, determinada por las condiciones ambientales
#
#     x_a=RH**(P_sat/P_amb)
# #-------Igualdad de la tasa de condensación a la de difusión de masa para hallar la aproximación de -----
#     #D= coeficiente de difusión del vapor de agua en el aire
#     #M= peso molecular del agua
#     #x_i= Fraccion molar en la interfase líquido-vapor de la gota determinada por la ecuación de Clausius-Clapeyron
#
#     #c= concentracion
#
#     delta=-D*M*c*((x_i-x_a)/t_c
#     return delta

#----------------FIN DE SECCIÓN PENDIENTE-----------------------


def R_int_liq_Wrikamanayake (r, h_i, theta, alfa, T_s):

#--Definición de h_i (coeficiente interfacial)

    h_i=((2*alfa)/(2-alfa))*(1/math.sqrt(2*math.pi*R_g*T_s))*((h_fg**2)/v_g*T_s)


#--------Resistencia de Interfase líquido-gas [K/W]----------------------------

    R_int=1/(h_i*2*math.pi*(r**2)*(1-math.cos(theta)))



    return R_int

def R_cond_gotas_Wrikamaniyake():
#-----Resistencia a la conducción de las gotas-------

R_drop=theta/(4*math.pi*k1*r*math.sin(theta))

#k1=Conductividad térmica del agua

#-----------RESISTENCIA TÉRMICA TOTAL---------------

R_total=(((1/R_dif)+(1/R_cond_aire))**(-1))+R_int+R_drop

#-----TRANSFERENCIA DE CALOR DEPENDIENTE DEL RADIO DE LA GOTA----------

q=(T_amb-T_s)/R_total

#-----FLUJO DE CALOR EN LA SUPERFICIE CONOCIENDO LA DISTRIBUCIÓN DE GOTAS

from scipy.integrate import quad

def f(r):
    return q(r)*N(r) # está incompleta

i, err= quad(f, r_min, r_max)
print (i)
print (err)
