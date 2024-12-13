import numpy as np
import matplotlib.pyplot as plt

def runge_kutta(f, t0, tf, y0, h):
    """
    Implementación del método de Runge-Kutta de cuarto orden.
    
    Args:
        f: Función que define el sistema de ecuaciones diferenciales.
        t0: Tiempo inicial.
        tf: Tiempo final.
        y0: Condición inicial.
        h: Paso de integración.
        
    Returns:
        t: Array con los tiempos.
        y: Array con las soluciones.
    """
    n = int((tf - t0) / h)
    t = np.linspace(t0, tf, n+1)
    y = np.zeros((n+1, len(y0)))
    y[0] = y0

    for i in range(n):
        k1 = h * f(t[i], y[i])
        k2 = h * f(t[i] + 0.5*h, y[i] + 0.5*k1)
        k3 = h * f(t[i] + 0.5*h, y[i] + 0.5*k2)
        k4 = h * f(t[i] + h, y[i] + k3)
        y_next = y[i] + (k1 + 2*k2 + 2*k3 + k4) / 6
        # Ajustar el ángulo para que esté en el intervalo [-pi, pi]
        y_next[0] = (y_next[0] + np.pi) % (2 * np.pi) - np.pi
        y[i+1] = y_next
    return t, y

def pendulo_amortiguado_forzado(t, y, omega_0, gamma, F_0, omega):
    """
    Ecuación diferencial para un péndulo amortiguado forzado.
    
    Args:
        t: Tiempo.
        y: Vector de estado [theta, theta_dot].
        omega_0: Frecuencia natural del péndulo.
        gamma: Coeficiente de amortiguamiento.
        F_0: Amplitud de la fuerza externa.
        omega: Frecuencia de la fuerza externa.
        
    Returns:
        dydt: Derivadas del vector de estado.
    """
    theta, theta_dot = y
    dydt = theta_dot
    dydotdt = -omega_0**2 * np.sin(theta) - gamma*theta_dot + omega_0**2 * F_0 * np.cos(omega*t)
    return np.array([dydt, dydotdt])

# Parámetros del péndulo
omega_0 = 1.0  # Frecuencia natural
gamma = 0.4    # Coeficiente de amortiguamiento
F_0 = 1.36     # Amplitud de la fuerza externa
omega = 2.0/3.0    # Frecuencia de la fuerza externa
T = 2*np.pi/omega

# Condiciones iniciales
theta_0 = 0.8
theta_dot_0 = 0.8
y0 = np.array([theta_0, theta_dot_0])

# Tiempo de integración
t0 = 0.0
tf = 500.0
h = T/150.0

# Resolviendo la ecuación diferencial
t, y = runge_kutta(lambda t, y: pendulo_amortiguado_forzado(t, y, omega_0, gamma, F_0, omega), t0, tf, y0, h)

# Guardar los datos en un archivo CSV
data = np.column_stack((t, y[:,0], y[:,1]))
np.savetxt('PenduloForzado.csv', data, delimiter=' ')