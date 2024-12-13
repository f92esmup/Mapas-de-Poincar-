import numpy as np
"""
Establecemos los valores de las constantes que se usan en las expresiones:
"""
#Valores para el péndulo amortiguado forzado:

CoeficienteAmortiguamiento = 0.2
FrecuenciaNatural = 1.0 #np.sqrt(Gravedad/Longitud[0]) #Frecuencia de oscilación del péndulo cuando no está sometido a ninguna fuerza.
Amplitud = 1.3 #Amplitud de la fuerza externa
FrecuenciaFuerzaExterna = 2.0  #Frecuencia de oscilación de la fuerza externa.

#Subamortiguado: CoeficientAmortiguamiento <  FrecuenciaNatural 
#Sobreamortiguadoy», CoeficientAmortiguamiento >  FrecuenciaNatural 
#Amortiguamiento critico CoeficientAmortiguamiento =  FrecuenciaNatural 

PeriodoFuerzaExterna = 2.0 * np.pi/ FrecuenciaFuerzaExterna

#Valores para el péndulo doble:
#Péndulo doble:
Masa = np.array([1,1]) #[Masa1,Masa2]
Longitud = np.array([1,1]) #Longitud de [varilla1, varilla2]
Gravedad = 9.8


"""
Vamos a definir todas las expresiones matemáticas que describen nuestros casos de estudio. Además también
incluimos el método Runge-Kutta de 4º orden y la función que selecciona los puntos pertenecientes al Mapa de Poincare
"""

#Función que implementa el método Runge-Kutta:
def RungeKutta(Funcion, CondicionesIniciales, TInicial, TFinal, TamañoPaso):
    """
    Entradas:

    Función -> Función de primer orden a la que se le aplicará el método.
    CondicionesIniciales -> Incluye la posición y velocidad inicial de nuestro sistema.
    TInicial -> Punto inicial, donde empezamos a evaluar la función
    TFinal -> Punto final, hasta donde evaluamos la función
    TamañoPaso -> Es el incremento de tiempo, la 'cantidad' de tiempo que avanzamos en cada paso. A partir
    de este valor y de los tiempos inicial y final es que vamos a obterner el número de pasos

    Salidas: 
    Tiempo -> Vector fila que agrupa todos los tiempos en los que se ha evaluado la función, 
    va desde TInicial hata TFinal (equiespaciados) con una dimensión igual al numero de pasos.
    """

    #Definimos el número de pasos:

    NumeroPasos = int((TFinal - TInicial) / TamañoPaso)

    #Incializamos las variables: 
    Tiempo = np.array(np.linspace(TInicial,TFinal,NumeroPasos)) #De esta forma nos evitamos tener que actualizar el valor en cada iteración del bucle for.
    Variables = np.array(np.zeros((NumeroPasos, len(CondicionesIniciales))))
    Theta = np.array(np.zeros((NumeroPasos,int(len(CondicionesIniciales)/2)))) 

    #Aplico las condiciones iniciales:
    Variables[0] = CondicionesIniciales
    Theta[0] = CondicionesIniciales[0]
    #Aplico el método Runge-Kutta:
    for i in range(NumeroPasos-1): #Ponemos NumeroPasos-1 por que sino al actualizar los datos del valor i+1 nos da un error de dimensiones.
        K1=np.array(Funcion(Tiempo[i], Variables[i]))
        K2=np.array(Funcion(Tiempo[i] + TamañoPaso/2, Variables[i] + (TamañoPaso/2) * K1))
        K3=np.array(Funcion(Tiempo[i] + TamañoPaso/2, Variables[i] + (TamañoPaso/2) * K2))
        K4=np.array(Funcion(Tiempo[i] + TamañoPaso, Variables[i] + TamañoPaso * K3))
        
        Variables[i+1,:] = Variables[i,:] + (TamañoPaso/6) * (K1+2*K2+2*K3+K4)     
        #Valores de Theta sin aplicarle el 'cambio de base'?
        Theta[i+1] =Theta[i] + (TamañoPaso/6) * (K1[0:1]+2*K2[0:1]+2*K3[0:1]+K4[0:1]) 
        #Ahora actualizamos los valores para la próima iteración. Tenemos que tener en cuenta la restricción de que theta pertenece a [-pi, pi]:
        #   (Variables[i+1,0] + np.pi) --> Este término nos asegura que, independientemente del ángulo, se encuentra en [0,2pi]
        #   con % calculamo el resto, es decir 'cuanto' se pasa del rango[0,2pi]
        #   Por último  se resta pi para estar en el rango [-pi,pi].
        
        Variables[i+1,0] = (Variables[i+1,0] + np.pi) % (2 * np.pi) - np.pi

        #print(Theta[i] - Variables[i,0])
        
    return Tiempo,Variables,Theta

#Función con la expresión que describe al péndulo forzado amortiguado:
def PenduloAmortiguadoForzado(Tiempo, Variables):
    """
    Entradas:
    Tiempo --> es el instante en que calculamos la función
    
    Variables --> Es un vector fila con [Theta, dTheta/dt]
    
    Salidas:
    
    Du1 -> Resultado de la primera ecuacion de primer grado
    Du2 -> Resultado de la segunda ecuación de sefundo grado
    
    -------------------------------------------------------------
    Como la ecuación del péndulo forzado amortiguado es de segundo orden,
    para aplicarle el método de Runge-Kutta tenemos que reducirla a un sistema 
    de dos ecuaciones de primer orden. Siendo
    Du1 = Theta
    Du2 = dTheta/dt
    """
    
    Du1= Variables[1]

    Du2 = -CoeficienteAmortiguamiento * Variables[1] - (FrecuenciaNatural**2  + 2* Amplitud * np.cos((FrecuenciaFuerzaExterna * Tiempo))) * np.sin(Variables[0])
    return np.array([Du1,Du2])
     
#Función que recoge los valores del mapa de Poincare:  
def PoincarePenduloForzado(Tiempo, Variables, Tolerancia, Theta):
    """
    Para obtener los puntos de poincare primero comprobamos que el tiempo es un múltiplo del periodo de la fuerza externa. 
    Para eso dividimos el tiempo entre el periodo con la division normal y la division que da un número entero, la diferencia
    entre estas dos valores debería ser cero para el caso en el que t = mT. Le aplicamos una tolerancia para no descartar todos los 
    puntos.
    
    Ademas, en lugar de usar el formato array de la librería Numpy usamos la listas de Phyton para ahorranos tener que establecer
    con anterioridad la dimensión del array, ya que no sabemos cual será.
    """
    t=0
    PuntosPoincare = []
    for i in range(len(Tiempo)):
        if abs(Tiempo[i] / PeriodoFuerzaExterna - Tiempo[i] // PeriodoFuerzaExterna) <= Tolerancia:
            #Voy a añadir un contador para eliminar el primer punto del que cumple la condicion:
            if t>1:
                PuntosPoincare.append(Variables[i, :])
            else:
                t=t+1
            
    #Guardo los datos en un fichero csv:
    np.savetxt('TiempoPosicionVelocidad.csv', np.column_stack((Tiempo, Variables,Theta)),delimiter='  ')
    np.savetxt('PuntosPoincare.csv', PuntosPoincare,delimiter='  ' )
    return

"""
"""  
      
#Función con la expresión que describe al péndulo doble:
def PenduloDoble(Tiempo,Variables):
    """
     Entradas:
     t -> Tiempo/paso en el que Runge-Kutta evalua la función. No depende de esta variable, por eso no se usa.
     Variables -> Vector de dos dimensiones con el valor de theta de la masa 1, theta de la masa 2, Momento de la masa 1 y Momento de la masa 2, respectivamente.
     
     Salidas:
     Vector con componentes:
     DTheta1 -> valor de la posición
     DTheta2 -> 
     DMomentoThta1 ->
     DMomentoTheta2 ->
    """
    DTheta1 = (Longitud[1] * Variables[2] - Longitud[0] * Variables[3] * np.cos(Variables[0]-Variables[1])) / (Longitud[0]**2 * Longitud[1] * (Masa[0] + Masa[1] * (np.sin(Variables[0] - Variables[1]))**2 ))
    
    DTheta2 = ((Masa[0] + Masa[1]) * Longitud[0] * Variables[3] - Masa[1] * Longitud[1] * Variables[2] * (np.cos(Variables[0] - Variables[1]))**2) / (Longitud[0] * Longitud[1]**2 * Masa[1] * (Masa[0] + Masa[1] * (np.sin(Variables[0] - Variables[1]))**2 ))
    
    A = (Variables[2] * Variables[3] * np.sin(Variables[0] - Variables[1])) / (Longitud[0] * Longitud[1] * (Masa[0] + Masa[1] * (np.sin(Variables[0] - Variables[1])**2)))
    
    B = ((Masa[1] * Longitud[1]**2 * Variables[2]**2 + (Masa[0] + Masa[1]) * Longitud[0]**2 * Variables[3]**2 -  2 * Masa[1] * Longitud[0] * Longitud[1] * Variables[2] * Variables[3] * np.cos( Variables[0] - Variables[1] )) / (2 * Longitud[0]**2 * Longitud[1]**2 * (Masa[0] + Masa[1] * (np.sin(Variables[0] - Variables[1])**2))) ) * np.sin(2 * (Variables[0] - Variables[1]))
    
    
    DMomentoTheta1 = -(Masa[0] + Masa[1]) * Gravedad * Longitud[0] * np.sin(Variables[0]) - A + B
    
    DMomentoTheta2 = - Masa[1] * Gravedad * Longitud[1] * np.sin(Variables[1]) + A - B
    
    return [DTheta1, DTheta2, DMomentoTheta1, DMomentoTheta2]

#Fución que contiene el Hamiltoniano del péndulo doble:
def HamiltonPenduloDoble(Variables):
    H = ((Masa[1] * Longitud[1]**2 * Variables[2]**2 + (Masa[0] + Masa[1]) * Longitud[0]**2 * Variables[3]**2 - 2*Masa[1] * Longitud[0] * Longitud[1] * Variables[2] * Variables[3] * np.cos(Variables[0] - Variables[1]) ) / (2 * Longitud[0]**2 * Longitud[1]**2 * Masa[1] * (Masa[0] + Masa[1]*(np.sin(Variables[0]-Variables[1]))**2)) ) - (Masa[0] + Masa[1]) * Gravedad * Longitud[0] * np.cos(Variables[0]) - Masa[1] * Gravedad * Longitud[1] * np.cos(Variables[1]) 
    return H

#Función que recoge los valores del mapa de Poincaré:
def PoincarePenduloDoble(Tiempo,Variables, EnergiaSistema,Tolerancia,Theta):
    PuntosPoincare = []

    for i in range(len(Variables)):
        H = HamiltonPenduloDoble(Variables[i, :])
        #print(H,"\n")
        if abs(EnergiaSistema-H) <= Tolerancia:
            PuntosPoincare.append(Variables[i, :])
    
    #Guardo los datos en un fichero csv:
    np.savetxt('TiempoPosicionVelocidad.csv', np.column_stack((Tiempo, Variables,Theta)),delimiter='  ')
    np.savetxt('PuntosPoincare.csv', PuntosPoincare,delimiter='  ' )

#Función para resolver una ecuación de segundo grado:
def EcuacionSegundoGrado(EnergiaSistema,CondicionInicial):
    #PARA EL CASO DE OBTENER EL VALOR DEl MomentoTheta2:
    #De la expresión general de la ecución de segundo grado sacamos los coeficientes:
    A = ((Masa[0] + Masa[1]) * Longitud[0]**2 ) / (2 * Longitud[0]**2 * Longitud[1]**2 * Masa[1] * (Masa[0] + Masa[1] * (np.sin(CondicionInicial[0] - CondicionInicial[1]))**2 ))
    
    B = - (2 * Masa[1] * Longitud[0] * Longitud[1] * CondicionInicial[2] * np.cos(CondicionInicial[0] - CondicionInicial[1])) / (2 * Longitud[0]**2 * Longitud[1]**2 * Masa[1] * (Masa[0] + Masa[1] * (np.sin(CondicionInicial[0] - CondicionInicial[1]))**2 ))
    
    C = ((Masa[1] * Longitud[1]**2 * CondicionInicial[2]**2) / (2 * Longitud[0]**2 * Longitud[1]**2 * Masa[1] * (Masa[0] + Masa[1] * (np.sin(CondicionInicial[0] - CondicionInicial[1]))**2 )) ) - EnergiaSistema - (Masa[0] + Masa[1]) * Gravedad * Longitud[0] * np.cos(CondicionInicial[0]) - Masa[1] * Gravedad * Longitud[1] * np.cos(CondicionInicial[1])
    
    ValorDentroDeLaRaiz = B**2 - 4*A*C
    #print(ValorDentroDeLaRaiz)
    Solucion1 =9999999999999
    if ValorDentroDeLaRaiz >= 0:
        Solucion1 = (- B + np.sqrt(ValorDentroDeLaRaiz)) / (2 * A)
        Solucion2 = (- B - np.sqrt(ValorDentroDeLaRaiz)) / (2 * A)
    
    return Solucion1