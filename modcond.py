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
