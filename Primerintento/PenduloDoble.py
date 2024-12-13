#Este es el programa principal del péndulo doble.
#Importo las librerias y los programas necesarios:
import numpy as np
import funciones as funciones
import RungeKutta as RungeKutta
import DibujafasesYPoincare as Dibuja
import matplotlib.pylab as plt
import subprocess
subprocess.run('clear', shell = True) #Para que se limpie la consola cada vez que se ejecute el programa.
#Declaro las variables que necesito para el RungeKutta:
PInicial= 0
PFinal =  100 
NumeroPasos = 5000

Tolerancia = 0.9 #Determina que error admito para seleccionar los puentos del mapa de Poincaré.


#Las referidas al péndulo doble:
EnergiaSistema = -20 #Energía del sistema del péndulo doble
CondicionInicial = np.array([1,0,0,0]) 

#Con esto selecciono para el pendulo doble el valor inicial del Momento generalizado 2 para el valor de energía dado
CondicionInicial[3] = funciones.EcuacionSegundoGrado(EnergiaSistema,CondicionInicial)

#Ahora llamo al método Runge-Kutta:

Tiempo, Variables = RungeKutta.RungeKutta (PInicial, PFinal, CondicionInicial, NumeroPasos)

#Guardo los datos obtenidos del Runge-Kutta en un archivo csv:
np.savetxt('PDPosicionVelocidad.csv', np.column_stack((Tiempo,Variables)),delimiter='  ' )

#Guardo los datos para el mapa de Poincaré
funciones.PoincarePenduloDoble(Variables, EnergiaSistema,Tolerancia)

"""
#Ahota vamos a dibujar el diagrama de fases del sistema:
Dibuja.DiagramaFases(Variables)
plt.tight_layout()
plt.show()

#Aqui dibujamos Theta frente a T:
Dibuja.ThetaFrenteaT(Tiempo,Variables)
plt.show()

#Dibuja Mapa de Poincaré:
Dibuja.PoincarePenduloDoble(Variables,EnergiaSistema,Tolerancia)
plt.tight_layout()
plt.show()
"""

#Comando del gnuplot
#plot 'PDPosicionVelocidad.csv' u 2:4 w l lc rgb "blue", 'PDPosicionVelocidad.csv' u 3:5 w l lc rgb "pink", 'PDPosicionVelocidadPoincare.csv' u 1:2 w p lc rgb "red"