#import Programación.funciones as funciones

import RungeKutta as RungeKutta
import funciones as funciones
import numpy as np
#Importamos lo necesario para dibujar las gráficas:
import DibujafasesYPoincare as Dibuja
import matplotlib.pyplot as plt
import subprocess
subprocess.run('clear', shell = True) #Para que se limpie la consola cada vez que se ejecute el programa.

#Ahora mismo estoy haciendo pruebas y lo que voy a hacer es comprobar que funciona el método Runge Kutta:

#Declaro las variables que necesito para llamar al programa Runge-Kutta:
PInicial= 0
PFinal =  100 #25 para el péndulo girando como loco
NumeroPasos = 10000

Tolerancia = 0.001

CondicionInicial = np.array([1.036,-2.150])


#Ahora llamo al método Runge-Kutta:

Tiempo, Variables = RungeKutta.RungeKutta (PInicial, PFinal, CondicionInicial, NumeroPasos)

"""
#Ahora voy a imprimir los valores del Runge-Kutta:

if len(CondicionInicial) == 2:
    #Para el péndulo simple
    print("| t | Theta | MomentoTheta |")
    print("-"*20)
    for i in range(NumeroPasos):
        print(f"{i:6d} | {Tiempo[i]:6} | {Variables[i,0]:6} | {Variables[i,1]:6}") #Theta2[i]:6 es para que al menos tenga 6 números y si no los rellena.

    print("-"*20)
else:
    #Para el péndulo Doble
    print("| t | Theta1 | Theta2 | MomentoTheta1 | MomentoTheta2 ")
    print("-"*40)
    for i in range(NumeroPasos):
        print(f"{i:6d} | {Tiempo[i]:6} | {Variables[i,0]:6} | {Variables[i,1]:6} | {Variables[i,2]:6} | {Variables[i,3]:6}") #Theta2[i]:6 es para que al menos tenga 6 números y si no los rellena.

    print("-"*40)
"""


#Guardo los datos obtenidos del Runge-Kutta en un archivo csv:
np.savetxt('PosicionVelocidad.csv', np.column_stack((Tiempo, Variables)), delimiter='  ' )

#Guardo los datos para el mapa de Poincaré
funciones.PoincarePenduloForzado(Tiempo, Variables,Tolerancia)

"""
#Ahota vamos a dibujar el diagrama de fases del sistema:
Dibuja.DiagramaFases(Variables)
plt.tight_layout()
plt.show()

#Aqui dibujamos Theta frente a T:
Dibuja.ThetaFrenteaT(Tiempo,Variables)
plt.show()

#Mapa de Poincaré:
Dibuja.PoincarePenduloForzado(Tiempo,Variables,Tolerancia,PFinal,PInicial)
plt.show()
"""

#plot 'PosicionVelocidad.csv' u 2:3 w l lc rgb "blue", 'PosicionVelocidadPoincare.csv' u 1:2 w p lc rgb "red"