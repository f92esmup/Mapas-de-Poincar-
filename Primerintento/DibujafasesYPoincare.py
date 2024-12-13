#En este programa voy a implementar el código para hacer una representación del diagrama de fases.
import funciones as funciones
import numpy as np
import matplotlib.pyplot as plt#Libreria para hacer las representaciones gráficas.
def ThetaFrenteaT(t,Variables):
     if len(Variables[0]) == 2: #Péndulo simple u amortiguado, como el simple lo he hecho para ver que funciona no es necesario ponerlo.
          plt.plot(t,Variables[:,0], label = 'Péndulo forzado')
          #Definimos las etiquetas:
          plt.title('Theta en función del tiempo')
          plt.ylabel("Ángulo Theta")
          plt.xlabel("Tiempo")
          plt.legend()
          #plt.show()
     else:
          fig, fig = plt.subplots(1,2)
          fig[0].plot(t,Variables[:, 0], label='Theta1', color = 'red')
          fig[1].plot(t, Variables[:, 1], label='Theta2', color = 'blue')
          for i in range(2):
            fig[i].set_title('Theta en función del tiempo')
            fig[i].set_xlabel("Tiempo")
            fig[i].set_ylabel("Ángulo Theta")
            fig[i].legend()



def DiagramaFases(Variables):
    if len(Variables[0]) == 2:  # Pendulo amortiguamiento forzado
        plt.plot(Variables[:, 0], Variables[:, 1], label='Péndulo forzado', color = 'blue')
        plt.title('Diagrama de fases')
        plt.xlabel("Ángulo Theta")
        plt.ylabel("Momento Lineal")
        plt.legend()
    else:  # Pendulo Doble
        fig, fases = plt.subplots(1, 2)
        fases[0].plot(Variables[:, 0], Variables[:, 2], label='Primero', color = 'red')
        fases[1].plot(Variables[:, 1], Variables[:, 3], label='Segundo', color = 'blue')
        for i in range(2):
            fases[i].set_title('Diagrama de fases')
            fases[i].set_xlabel("Ángulo Theta")
            fases[i].set_ylabel("Momento Lineal")
            fases[i].legend()
    
    #Para ver los valores por pantalla usamos:
    #plt.show()



def PoincarePenduloDoble(Variables, EnergiaSistema,Tolerancia):
    PuntosPoincare = []

    for i in range(len(Variables)):
        H = funciones.HamiltonPenduloDoble(Variables[i, :])
        #print(H,"\n")
        if abs(EnergiaSistema-H) <= Tolerancia:
            PuntosPoincare.append(Variables[i, :])
            
    fig, Mp = plt.subplots()
    Mp.plot([punto[0] for punto in PuntosPoincare], [punto[2] for punto in PuntosPoincare], label='Primero', color='red',linestyle='none', marker='o', markersize='4')
    Mp.plot([punto[1] for punto in PuntosPoincare], [punto[3] for punto in PuntosPoincare], label='Segundo', color='blue',linestyle='none', marker='o',markersize='4')
    Mp.set_title('Mapa Poincaré')
    Mp.set_xlabel("Ángulo Theta")
    Mp.set_ylabel("Momento Lineal")
    Mp.legend()



def PoincarePenduloForzado(Tiempo,Variables,Tolerancia,PFinal,PInicial):
    T = funciones.PeriodoFuerzaExterna
    PuntosPoincare = []
    
    TamañoPaso = (PFinal-PInicial)/(len(Tiempo)-1)
    Tiempo = np.array([t * TamañoPaso for t in Tiempo])
        
    for i in range(len(Tiempo)):
        if abs(Tiempo[i] / T - Tiempo[i] // T) <= Tolerancia:
            PuntosPoincare.append(Variables[i, :])
    
    plt.plot([punto[0] for punto in PuntosPoincare], [punto[1] for punto in PuntosPoincare], label='Pendulo Forzado', color='red',linestyle='none', marker='o', markersize='4')
    plt.title('Mapa de Poincaré')
    plt.xlabel('Theta')
    plt.ylabel('Velocidad angular')
    plt.legend()

    
        