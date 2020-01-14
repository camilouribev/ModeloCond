import numpy as np
import scipy
import scipy.integrate as integrate
import math

#Método de resistencias
#Se irán haciendo las resistencias desde el exterior de la gota hasta el sustrato
def R_dif_Wrikamaniyake (T_amb, T_int, t_c, h_fg):
#----------------Resistencia a la difusión del vapor en la capa de no condensables----------------

#R_dif=resistencia térmica por la capa de difusión
#T_amb=Temperatura ambiente
#T_int=Temperatura de la gota
#t_c=Tasa de condensación. Hallada experimentalmente
#h_fg= Calor latente de condensación

    R_dif=(T_amb-T_int)/(t_c*h_fg)
    return R_dif

def R_cond_aire_Wrikamanayake (K_aire, delta, S_r):
#--------------Resistencia a la conducción del aire---------------------
# K_aire=Conductividad térmica del K_aire
#delta= grosor de la capa límite
#S_r=Área proyectada de la gota

    R_cond_aire=delta/(k_aire*S_r)
    return R_cond_aire

def Capa_limite (h_fg, R, T_int, T_boil,D, M, c, x_i, x_a, Delta):
#------Inicialmente, Ecuación Clausius-Clapeyron para determinar x_i, termino necesario------------------
    #T_boil=Punto de fusión a presión atmosférica (P_atm)
    #R= 

    x_i=math.exp(-(h_fg/R)*((1/T_int)-(1/T_boil)))

#-------Igualdad de la tasa de condensación a la de difusión de masa para hallar la aproximación de -----
    #D= coeficiente de difusión del vapor de agua en el aire
    #M= peso molecular del agua
    #x_i= Fraccion molar en la interfase líquido-vapor de la gota determinada por la ecuación de Clausius-Clapeyron
    #x_a=Fracción molar en el agua, determinada por las condiciones ambiente
    #c= concentracion

    delta=-D*M*c*((x_i-x_a)/t_c
    return delta

#--------Resistencia de Interfase líquido-gas----------------------------

R_int=1/(h_i*2*math.pi*(r**2)*(1-math.cos(theta)))

#h_i=Coeficiente interfacial de transferencia de Calor
#r=radio de la gota
#theta=Ángulo de contacto

#--Definición de h_i (coeficiente interfacial)

h_i=((2*alfa)/(2-alfa))*(1/math.sqrt(2*math.pi*R_g*T_s))*((h_fg**2)/v_g*T_s)

#alfa=Coeficiente de acomodación o condensación;
# fracción de moléculas de vapor que se mueven en la gota líquida durante el cambio de fase
#debe ser de entre 0.02 a 0.04 para agua. En Wrikamanayake et al. se usó 0.02
#T_s=Temperatura de superficie(no está claro)

#-----Resistencia a la conducción de las gotas-------

R_drop=theta/(4*math.pi*k1*r*math.sin(theta))

#k=Conductividad térmica del agua

#-----------RESISTENCIA TÉRMICA TOTAL---------------

R_total=(((1/R_dif)+(1/R_cond_aire))**(-1))+R_int+R_drop

#-----TRANSFERENCIA DE CALOR DEPENDIENTE DEL RADIO DE LA GOTA----------

q=(T_amb-T_s)/R_total

#-----FLUJO DE CALOR EN LA SUPERFICIE CONOCIENDO LA DISTRIBUCIÓN DE GOTAS

from scipy.intgrate import quad

def f(r):
    return q(r)*N(r) # está incompleta

i, err= quad(f, r_min, r_max)
print (i)
print (err)
