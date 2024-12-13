#En este archivo voy a definir la función a la que llamaré en varias ocasiones en el método Runge-Kutta:
import numpy as np



#Péndulo doble:
Masa = np.array([1,1]) #[Masa1,Masa2]
Longitud = np.array([1,1]) #[varilla1, varilla2]
Gravedad = 9.8

#Valores para el péndulo amortiguado forzado:

CoeficientAmortiguamiento = 0.2
FrecuenciaNatural = 1 #np.sqrt(Gravedad/Longitud[0]) #Frecuencia de oscilación del péndulo cuando no está sometido a ninguna fuerza.
Amplitud = 1.036 #Amplitud de la fuerza externa
FrecuenciaExterna = 2  #Frecuencia de oscilación de la fuerza externa.
#Subamortiguado: CoeficientAmortiguamiento <  FrecuenciaNatural 
#Sobreamortiguadoy», CoeficientAmortiguamiento >  FrecuenciaNatural 
#Amortiguamiento critico CoeficientAmortiguamiento =  FrecuenciaNatural 
PeriodoFuerzaExterna = 2*np.pi/ FrecuenciaExterna

def PenduloSimple(t , Variables): #Variables= [Theta, MomentoTheta]
    DTheta = Variables[1]/(Masa[0] * Longitud[0]**2)
    DMomentoTheta = - Masa[0] * Longitud[0] * Gravedad * np.sin(Variables[0])
    return [DTheta, DMomentoTheta]


#Para el péndulo doble tendremos que resolver con Runge-Kutta el sistema de EDOS de primer orden.

def PenduloDoble(t,Variables):
    DTheta1 = (Longitud[1] * Variables[2] - Longitud[0] * Variables[3] * np.cos(Variables[0]-Variables[1])) / (Longitud[0]**2 * Longitud[1] * (Masa[0] + Masa[1] * (np.sin(Variables[0] - Variables[1]))**2 ))
    
    DTheta2 = ((Masa[0] + Masa[1]) * Longitud[0] * Variables[3] - Masa[1] * Longitud[1] * Variables[2] * (np.cos(Variables[0] - Variables[1]))**2) / (Longitud[0] * Longitud[1]**2 * Masa[1] * (Masa[0] + Masa[1] * (np.sin(Variables[0] - Variables[1]))**2 ))
    
    A = (Variables[2] * Variables[3] * np.sin(Variables[0] - Variables[1])) / (Longitud[0] * Longitud[1] * (Masa[0] + Masa[1] * (np.sin(Variables[0] - Variables[1])**2)))
    
    B = ((Masa[1] * Longitud[1]**2 * Variables[2]**2 + (Masa[0] + Masa[1]) * Longitud[0]**2 * Variables[3]**2 -  2 * Masa[1] * Longitud[0] * Longitud[1] * Variables[2] * Variables[3] * np.cos( Variables[0] - Variables[1] )) / (2 * Longitud[0]**2 * Longitud[1]**2 * (Masa[0] + Masa[1] * (np.sin(Variables[0] - Variables[1])**2))) ) * np.sin(2 * (Variables[0] - Variables[1]))
    
    
    DMomentoTheta1 = -(Masa[0] + Masa[1]) * Gravedad * Longitud[0] * np.sin(Variables[0]) - A + B
    
    DMomentoTheta2 = - Masa[1] * Gravedad * Longitud[1] * np.sin(Variables[1]) + A - B
    
    return [DTheta1, DTheta2, DMomentoTheta1, DMomentoTheta2]    

def PenduloAmortiguamientoForzado(Tiempo,Variables): 
    #Hago un sistema donde reduzco el orden de la ecuación de segundo orden a una de primer orden con 
    #u1 = Theta
    #u2 = DTheta/Dt
     
    #Tiempo = t * TamañoPaso #Este sería el tiempo discreto, y no debeo usarlo aqui
    
    Du1= Variables[1]

    Du2 = - CoeficientAmortiguamiento * Variables[1] - (FrecuenciaNatural**2 + 2 * Amplitud * np.cos(FrecuenciaExterna * Tiempo) ) * np.sin(Variables[0])
    return [Du1, Du2]   

def HamiltonPenduloDoble(Variables):
    H = ((Masa[1] * Longitud[1]**2 * Variables[2]**2 + (Masa[0] + Masa[1]) * Longitud[0]**2 * Variables[3]**2 - 2*Masa[1] * Longitud[0] * Longitud[1] * Variables[2] * Variables[3] * np.cos(Variables[0] - Variables[1]) ) / (2 * Longitud[0]**2 * Longitud[1]**2 * Masa[1] * (Masa[0] + Masa[1]*(np.sin(Variables[0]-Variables[1]))**2)) ) - (Masa[0] + Masa[1]) * Gravedad * Longitud[0] * np.cos(Variables[0]) - Masa[1] * Gravedad * Longitud[1] * np.cos(Variables[1]) 
    return H

def EcuacionSegundoGrado(EnergiaSistema,CondicionInicial):
    #PARA EL CASO DE OBTENER EL VALOR DEl MomentoTheta2:
    #De la expresión general de la ecución de segundo grado sacamos los coeficientes:
    A = ((Masa[0] + Masa[1]) * Longitud[0]**2 ) / (2 * Longitud[0]**2 * Longitud[1]**2 * Masa[1] * (Masa[0] + Masa[1] * (np.sin(CondicionInicial[0] - CondicionInicial[1]))**2 ))
    
    B = - (2 * Masa[1] * Longitud[0] * Longitud[1] * CondicionInicial[2] * np.cos(CondicionInicial[0] - CondicionInicial[1])) / (2 * Longitud[0]**2 * Longitud[1]**2 * Masa[1] * (Masa[0] + Masa[1] * (np.sin(CondicionInicial[0] - CondicionInicial[1]))**2 ))
    
    C = ((Masa[1] * Longitud[1]**2 * CondicionInicial[2]**2) / (2 * Longitud[0]**2 * Longitud[1]**2 * Masa[1] * (Masa[0] + Masa[1] * (np.sin(CondicionInicial[0] - CondicionInicial[1]))**2 )) ) - EnergiaSistema - (Masa[0] + Masa[1]) * Gravedad * Longitud[0] * np.cos(CondicionInicial[0]) - Masa[1] * Gravedad * Longitud[1] * np.cos(CondicionInicial[1])
    
    ValorDentroDeLaRaiz = B**2 - 4*A*C
    #print(A,B,C)
    Solucion1 =9999
    if ValorDentroDeLaRaiz >= 0:
        Solucion1 = (- B + np.sqrt(ValorDentroDeLaRaiz)) / (2 * A)
        Solucion2 = (- B - np.sqrt(ValorDentroDeLaRaiz)) / (2 * A)
    
    return Solucion1

def PoincarePenduloDoble(Variables, EnergiaSistema,Tolerancia):
    PuntosPoincare = []

    for i in range(len(Variables)):
        H = HamiltonPenduloDoble(Variables[i, :])
        #print(H,"\n")
        if abs(EnergiaSistema-H) <= Tolerancia:
            PuntosPoincare.append(Variables[i, :])
    
    #Guardo los datos en un fichero csv:
    np.savetxt('PDPosicionVelocidadPoincare.csv', PuntosPoincare,delimiter='  ' )


    """        
    fig, Mp = plt.subplots()
    Mp.plot([punto[0] for punto in PuntosPoincare], [punto[2] for punto in PuntosPoincare], label='Primero', color='red',linestyle='none', marker='o', markersize='4')
    Mp.plot([punto[1] for punto in PuntosPoincare], [punto[3] for punto in PuntosPoincare], label='Segundo', color='blue',linestyle='none', marker='o',markersize='4')
    Mp.set_title('Mapa Poincaré')
    Mp.set_xlabel("Ángulo Theta")
    Mp.set_ylabel("Momento Lineal")
    Mp.legend()
    """
    
    
def PoincarePenduloForzado(Tiempo,Variables,Tolerancia):
    PuntosPoincare = []
    """
    TamañoPaso = (PFinal-PInicial)/(len(Tiempo)-1)
    Tiempo = np.array([t * TamañoPaso for t in Tiempo])  
    """
    #Esta es la primera forma que hice para obtener los puntos:
    for i in range(len(Tiempo)):
        if abs(Tiempo[i] / PeriodoFuerzaExterna - Tiempo[i] // PeriodoFuerzaExterna) <= Tolerancia:
            PuntosPoincare.append(Variables[i, :])
            
    #Guardo los datos en un fichero csv:
    np.savetxt('PosicionVelocidadPoincare.csv', PuntosPoincare,delimiter='  ' )
    
    """
    plt.plot([punto[0] for punto in PuntosPoincare], [punto[1] for punto in PuntosPoincare], label='Pendulo Forzado', color='red',linestyle='none', marker='o', markersize='4')
    plt.title('Mapa de Poincaré')
    plt.xlabel('Theta')
    plt.ylabel('Velocidad angular')
    plt.legend()
    """
    
