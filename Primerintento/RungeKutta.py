
#En este programa voy a aplicar el método Runge-Kutta de 4º orden.
#La idea es definir un método general al cual se le llamara desde otro programa
#para no tener que repetirlo constantemente.

#Voy a importar el archivo que tiene la expresión de la función a la que le voy a aplicar el método:
import funciones as funciones

#Importo la lista numpy para poder crear una lista de zeros:
import numpy as np


#Definimos la función:
    #Las entradas de este método son: intervalo de la función [PInicial, Pfinal], condición inicial (CondicionInicial) y
    # el número de pasos (NumeroPasos)
    
def RungeKutta(PInicial, PFinal, CondicionInicial, NumeroPasos):

    #Calculamos el tamaño del paso "h":

    h= (PFinal - PInicial) / NumeroPasos

    #Para un función genérica f= f(x,y) a la que le aplicamos el método RungeKutta, inicializamos las variables:
    #Creo una lista vacia porque quiero guardar los valores de todos los pasos:
    
    x = np.zeros(NumeroPasos + 1) #El +1 es por que si no en el bucle for tenemos un error de dimensiones 
    y = np.zeros((NumeroPasos + 1, len(CondicionInicial)))
    
    #Como la posición de los vectores comienza desde 0
    x[0] = PInicial
    y[0,:] = CondicionInicial
    

    for i in range(NumeroPasos):
        if len(y[0]) == 2: #Pendulo forzado
            K1 = np.array (funciones.PenduloAmortiguamientoForzado(x[i],y[i,:]))
            K2 = np.array (funciones.PenduloAmortiguamientoForzado(x[i] + h/2, y[i,:] + (h/2)*K1))
            K3 = np.array (funciones.PenduloAmortiguamientoForzado(x[i] + h/2, y[i,:] + (h/2)*K2))
            K4 = np.array (funciones.PenduloAmortiguamientoForzado(x[i] + h, y[i,:]+ h*K3))
        else:
            K1 = np.array (funciones.PenduloDoble(x[i],y[i,:]))
            K2 = np.array (funciones.PenduloDoble(x[i] + h/2, y[i,:] + (h/2)*K1))
            K3 = np.array (funciones.PenduloDoble(x[i] + h/2, y[i,:] + (h/2)*K2))
            K4 = np.array (funciones.PenduloDoble(x[i] + h, y[i,:]+ h*K3))
        #Actualizamos los valores de x e y:
        
        y[i+1,:] = y[i,:] + (h/6)*(K1 + 2*K2 + 2*K3 + K4)
        x[i+1] = x[i] + h
        
    return [x,y]
# x representa el paso del tiempo, es el "salto" en los que se descretiza el intervalor [a,b]. 
#La x es la variable independiente que en nuestro caso es el tiempo, por eso le voy a cambiar el nombre
#a x= t 

# y es el valor de la función para cada tiempo x.
# La y son las variables dependientes, en nuestro caso tenemos varias
#Definimos k1,k2,k3,k4, la función, las condiciones iniciales y Y como un array para poder emplear el método de Runge-Kutta para un sistema de ecuaciones.

