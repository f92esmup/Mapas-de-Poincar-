import numpy as np
import subprocess
subprocess.run('clear', shell=True)

# FUNCIONES:

def RungeKutta(PInicial,TamañoPaso, CondicionesIniciales, NumeroPasos):
    #Voy a definir un tamaño de paso que en funcióndel periodo de la fuerza externa para obtener valores que luego podré
    #manejar mejor para sacar el mapa de Poincare.
    
    
    Tiempo = np.array(np.zeros((NumeroPasos)))
    #Tiempo = np.array(np.linspace(PInicial, PFinal, NumeroPasos))
    Variables = np.array(np.zeros((NumeroPasos,2)))
    

    
    Variables[0] = CondicionesIniciales
    Tiempo[0] = PInicial

    for i in range(NumeroPasos-1):
        K1=np.array(Funcion(Tiempo[i], Variables[i]))
        K2=np.array(Funcion(Tiempo[i] + TamañoPaso/2, Variables[i] + (TamañoPaso/2) * K1))
        K3=np.array(Funcion(Tiempo[i] + TamañoPaso/2, Variables[i] + (TamañoPaso/2) * K2))
        K4=np.array(Funcion(Tiempo[i] + TamañoPaso, Variables[i] + TamañoPaso * K3))
        
        Variables[i+1,:] = Variables[i,:] + (TamañoPaso/6) * (K1+2*K2+2*K3+K4)
        Tiempo[i+1] = Tiempo[i] + TamañoPaso
    
    return [Tiempo,Variables]

def Funcion(Tiempo,Variables):
    Du1= Variables[1]
    
    #Du2 = -CoeficienteAmortiguamiento * Variables[1] - (FrecuenciaNatural**2  + 2* Amplitud * np.cos((FrecuenciaFuerzaExterna * Tiempo))) * np.sin(Variables[0])
    
    #Du2 = -CoeficienteAmortiguamiento * Variables[1] - FrecuenciaNatural**2 * np.sin(Variables[0])  -  2*Amplitud * np.cos(FrecuenciaFuerzaExterna * Tiempo)     
    Du2 = -CoeficienteAmortiguamiento * Du1 - FrecuenciaNatural**2 * np.sin(Variables[0]) + FrecuenciaNatural**2 * Amplitud * np.cos(FrecuenciaFuerzaExterna * Tiempo)
    return  [Du1, Du2]

def PoincarePenduloForzado(Tiempo,Variables,Tolerancia):
    #Esta es la primera forma que hice para obtener los puntos:
    PuntosPoincare = []
    for i in range(2,len(Tiempo)):
        if abs(Tiempo[i] / PeriodoFuerzaExterna - Tiempo[i] // PeriodoFuerzaExterna) <= Tolerancia:
            PuntosPoincare.append(Variables[i, :])
            
    #Guardo los datos en un fichero csv:
    np.savetxt('PenduloForzadoPoincare.csv', PuntosPoincare,delimiter='  ' )
"""
def Funcion(t , Variables): #Variables= [Theta, MomentoTheta]
    DTheta = Variables[1]/(1 * 1**2)
    DMomentoTheta = - 1 * 1 * 9.8 * np.sin(Variables[0])
    return [DTheta, DMomentoTheta]
"""

FrecuenciaNatural=1
FrecuenciaFuerzaExterna= 2/3
#Subamortiguado: CoeficientAmortiguamiento <  FrecuenciaNatural 
#Sobreamortiguadoy», CoeficientAmortiguamiento >  FrecuenciaNatural 
#Amortiguamiento critico CoeficientAmortiguamiento =  FrecuenciaNatural 
PeriodoFuerzaExterna = 2*np.pi/FrecuenciaFuerzaExterna


Amplitud=1.36  #0.3 -->Periodic state  1.23 --> Quasi-periodic state 1.36 --> Chaotic state 
CoeficienteAmortiguamiento =0.4
CondicionesIniciales = np.array([0.8,0.8])
#Oscilador Amortiguado ---> CORRECTO ->Por lo tanto el sin(Theta) NO lleva el np.radians().
#Oscilador Forzado --->  Creo que el otro si debe tener el radians.


PInicial = 0
PFinal = PeriodoFuerzaExterna/150
NumeroPasos =10000
Tolerancia = 0.01


Tiempo, Variables = RungeKutta(PInicial,PFinal, CondicionesIniciales, NumeroPasos)



np.savetxt('PenduloForzado.csv',np.column_stack((Tiempo, Variables)),delimiter='  ')

PoincarePenduloForzado(Tiempo,Variables,Tolerancia)

#plot 'PenduloForzado.csv' u 2:3 w l lc rgb "blue", 'PenduloForzadoPoincare.csv' u 1:2 w p lc rgb "red"
#plot 'PenduloForzado.csv' every ::0::0 u 2:3 w p lc rgb "yellow",'PenduloForzado.csv' every ::1 u 2:3 w l lc rgb "blue", 'PenduloForzadoPoincare.csv' u 1:2 w p lc rgb "red"

"""Solo pueden estar fallando dos cosas:
-El método Runge-Kutta --> Voy a aplicarlo al péndulo simple y ver que valores obtengo
-La función del péndulo forzado amortiguado
LO ÚNICO QUE SE ME OCURRE ES QUE LA ECUACIÓN DEL PÉNDULO FORZADO SEA INCORRECTA, NO ESTÁ BIEN CALCULADA Y ES OTRA EXPRESIÓN SIN EL SENO MULTIPLICANDO."""