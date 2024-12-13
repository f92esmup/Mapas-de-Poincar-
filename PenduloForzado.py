import numpy as np
import Funciones as Funciones
import subprocess
subprocess.run('clear', shell=True) #Para limpiar la terminal cada vez que se ejecute el programa.

#Declaro las variables que necesito para llamar al programa Runge-Kutta:
#Intervalo de Tiempo:
TInicial = 0
TFinal =  600

"""
En lugar de dar el número de pasos al nuestro método le vamos a dar el incremento de tiempo como función del periodo
de la fuerza externa. De esta manera nos aseguramos que los tiempo elegidos recorreran en algun momento por lo menos un periodo
de la fuerza externa.
"""
TamañoPaso = Funciones.PeriodoFuerzaExterna/150 #150 es un número aleatoriamente alto para que el incremeto de tiempo sea considerablemente pequeño.

#Establezco las condiciones inciales de la posición y la velocidad en un vector fila.
CondicionesIniciales = np.array([0.8,0.8])

#Establecemos la Tolerancia para la que consideramos que un punto pertenece al plano de Poincaré.
Tolerancia = 0.01

#Aplico el método de Runge-Kutta:
Tiempo, Variables, Theta = Funciones.RungeKutta(lambda Tiempo, Variables: Funciones.PenduloAmortiguadoForzado(Tiempo,Variables), CondicionesIniciales, TInicial, TFinal, TamañoPaso)

#Guardo la lista de resultados en un fichero.csv para imprimirlos por pantalla posteriormente con gnuplot:
Funciones.PoincarePenduloForzado(Tiempo, Variables, Tolerancia, Theta)
"""
Para graficar con gnuplot el diagrama de fases y mapa de poincare: NO USAR LINEAS, PORQUE UNEN PUNTOS EN EXTREMOS OPUESTOS Y ENSUCIA LA GRÁFICA
plot 'TiempoPosicionVelocidad.csv' u 2:3 w p lc rgb "blue",'PuntosPoincare.csv' u 1:2 w p lc rgb "red"

Theta vs tiempo:
plot 'TiempoPosicionVelocidad.csv' u 1:4 w l lc rgb "blue" 

Velocidad vs tiempo: 
plot 'TiempoPosicionVelocidad.csv' u 2:4 w l lc rgb "blue" 
"""