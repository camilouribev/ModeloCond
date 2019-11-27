import numpy as np
import pandas as pd
import math
#Método de resistencias
#Se irán haciendo las resistencias desde el exterior de la gota hasta el sustrato

#----------------Resistencia a la difusión del vapor en la capa de no condensables----------------

R_dif=(T_amb-T_int)/(m_p*h_fg)

#R_dif=resistencia térmica por la capa de difusión
#T_amb=Temperatura ambiente
#T_int=Temperatura de la gota
#t_c=Tasa de condensación. Hallada experimentalmente
#h_fg= Calor latente de condensación

#--------------Resistencia a la conducción del aire---------------------

R_cond_aire=delta/(k_aire*S_r)

# K_aire=Conductividad térmica del K_aire
#delta= grosor de la capa límite
#S_r=Área proyectada de la gota

#-------Igualdad de la tasa de condensación a la de difusión de masa para hallar la aproximación de delta-----

t_c=-D*M*c*((x_i-x_a)/delta)

#D= coeficiente de difusión del vapor de agua en el aire
#M=peso molecular del agua
#x_i= Fraccion molar en la interfase líquido-vapor de la gota determinada por la ecuación de Clausius-Clapeyron
#x_a=Fracción molar en el agua, determinada por las condiciones ambiente

#------Ecuación Clausius-Clapeyron para determinar x_i------------------

x_i=math.exp(-(h_fg/R)*((1/T_int)-(1/T_boil)))

#T_boil=Punto de fusión a presión atmosférica (P_atm)
#

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

R_cond
