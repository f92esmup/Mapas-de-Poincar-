import numpy as np
import Funciones as Funciones
import subprocess
subprocess.run('clear',shell=True)


"""
    TENGO QUE EXPLICAR ESTE CÓDIGO:
"""


#Declaro las variables que necesito para llamar al programa Runge-Kutta:
#Intervalo de Tiempo:
TInicial= 0
TFinal =  60
"""
En lugar de dar el número de pasos al nuestro método le vamos a dar el incremento de tiempo.
"""
TamañoPaso = 0.001 #Según el Capitulo 6 del libro.
#Establecemos la Tolerancia para la que consideramos que un punto pertenece al plano de Poincaré.
Tolerancia = 0.1

#Ene----
EnergiaSistema= 40  -29.4

#Establezco las condiciones inciales de la posición y la velocidad en un vector fila.
CondicionesIniciales = np.array([0.0, 0.0, 0.0, 0.08])

#Con esto selecciono para el pendulo doble el valor inicial del Momento generalizado 2 para el valor de energía dado
CondicionesIniciales[3] = Funciones.EcuacionSegundoGrado(EnergiaSistema,CondicionesIniciales)


#Aplico el método de Runge-Kutta:
Tiempo, Variables, Theta = Funciones.RungeKutta(lambda Tiempo, Variables: Funciones.PenduloDoble(Tiempo,Variables), CondicionesIniciales, TInicial, TFinal, TamañoPaso)

#Guardo la lista de resultados en un fichero.csv para imprimirlos por pantalla posteriormente con gnuplot:
Funciones.PoincarePenduloDoble(Tiempo, Variables,EnergiaSistema, Tolerancia, Theta)

# TiempoPosicionVelocidad = [Tiempo, Theta1, Theta2, Momento 1, Momento 2, Thetaprima1, Thetaprima2]
"""
Para graficar con gnuplot el diagrama de fases y mapa de poincare:
Momento1 vs Theta1:
plot 'TiempoPosicionVelocidad.csv' u 2:4 w p lc rgb "blue",'PuntosPoincare.csv' u 1:3 w p lc rgb "red"

Momento2 vs Theta2:
plot 'TiempoPosicionVelocidad.csv' u 3:5 w l lc rgb "blue",'PuntosPoincare.csv' u 2:4 w p lc rgb "red"

Theta1 vs tiempo:
plot 'TiempoPosicionVelocidad.csv' u 1:6 w l lc rgb "blue" 

Theta2 vs tiempo:
plot 'TiempoPosicionVelocidad.csv' u 1:7 w l lc rgb "blue" 

Momento1 vs tiempo: 
plot 'TiempoPosicionVelocidad.csv' u 1:4 w l lc rgb "blue" 

Momento2 vs tiempo: 
plot 'TiempoPosicionVelocidad.csv' u 1:5 w l lc rgb "blue" 
"""

