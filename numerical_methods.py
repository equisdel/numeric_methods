# This code implements three well known numeric methods 
#   used to approximate a function from its derivative:
#   1. Euler's Method
#   2. Improved Euler's Method
#   3. Runge-Kutta Method

from matplotlib import pyplot as plt
import math

# To test and compare the algorithms, the following problem is proposed:
#   A spinach pie is taken out of the oven. Its temperature at that moment is 160°C. 
#   It's a summer day, with an ambient temperature of 35°C. 
#   Its temperature at 12:15 PM is 120°C, while at 12:20 PM the temperature is 90°C.
#   Approximate the temperature function over time, y(t), by deducing its rate of change - or derivative - y'(t).
# However, it can be adapted to any similar problem by simply changing y'(t), initial data and graphication parameters.


# FUNCTIONS: derivative and its solution

# y'(x): rate of change of pie's temperature (obtained analitically from given data at @@)
def Y(x,y):
    k = - (math.log(11/17)) / 5
    return  - k * (y - 35)

# y(x): solution - pie's temperature over time (optional, only for comparison purposes).
def solution(x):
    k = - (math.log(11/17)) / 5
    return math.exp(-k*x) * 125 + 35

# Adjust solution to the global granularity (h timeStep, n iterations)
def g_solution(x0, y0, h, n):
    xn = [x0]
    yn = [y0]
    for t in range(n):
        x0 = x0 + h
        y0 = solution(x0)
        xn.append(x0)
        yn.append(y0)
    return xn, yn


# NUMERICAL METHODS

# [1]: Euler's Method
def euler(x0, y0, h, Y, n):
    xn = [x0]
    yn = [y0]
    for t in range(n):
        x0 = x0 + h
        y0 = y0 + h*Y(x0,y0)
        xn.append(x0)
        yn.append(y0)
    return xn, yn 

# [2]: Improved Euler's Method
def eulerMejorado(x0, y0, h, Y, n):
    xn = [x0]
    yn = [y0]
    y_aux = y0
    for t in range(n):
        x0 = x0 + h
        y_aux = y_aux + h * Y(x0,y_aux)
        y0 = y0 + ((h/2)*(Y(x0,y0)+Y(x0+h,y_aux)))
        y_aux = y0
        xn.append(x0)
        yn.append(y0)
    return xn, yn;

# [3]: Runge-Kutta Method
def rungeKutta(x0, y0, h, Y, n):
    xn=[x0]
    yn=[y0]
    for t in range(n):
        x0 = x0 + h
        k1 = Y(x0,y0)
        k2 = Y(x0+(h/2), y0+((h/2)*k1))
        k3 = Y(x0+(h/2), y0+((h/2)*k2))
        k4 = Y(x0+h, y0+(h*k3))
        y0 = y0 +((h/6)*(k1+(2*k2)+(2*k3)+k4))
        xn.append(x0)
        yn.append(y0)
    return xn, yn;


# Display graphics using MathPlot Library
def graficate(x1, y1, x2, y2, x3, y3, xs, ys):
    plt.plot(x1, y1, linewidth=1, color='red', label='Euler')      # Euler
    plt.plot(x2, y2, linewidth=1, color='orange', label='Improved Euler')   # Improved Euler
    plt.plot(x3, y3, linewidth=1, color='cyan', label='Runge-Kutta')     # Runge-Kutta
    plt.plot(xs, ys, linewidth=1, color='green', label='Solution')    # Solution
    plt.legend()
    plt.ylim(0,170)         # interval on the y axis
    plt.xlim(0,47.1)     # interval on the x axis
    plt.grid(True)
    plt.show()

# Calculate number of iterations (n) based on timestep (h) and a given interval on the x axis.
def N(xi, xf, h):
    return (int)((xf-xi)/h)


# MAIN

x0,y0 = 0,160   # y(x0) = y0  -- data given in the problem statement. 
xf = 47         # y(xf) = ¿?  -- what you want to find out (particular solution).

# Try smaller values of h for better precision in all three methods.
h = 10             # timestep / granularity of algorithms.
n = N(x0,xf,h)     # number of iterations for such granularity.

# Testing all methods
xn1, yn1 = euler(x0, y0, h, Y, n)
xn2, yn2 = eulerMejorado(x0, y0, h, Y, n)
xn3, yn3 = rungeKutta(x0, y0, h, Y, n)
xs, ys = g_solution(x0, y0, h, n)

# Displaying results with respective errors
n_digitos = 6   # used to round float values
print("\nSIMULACION {n:",n,"; h:",h,";}\n")
print("Resultados de y(",xf,"):")
print("Método Euler: (", f"{xn1[n]:.{n_digitos}f}",";", f"{yn1[n]:.{n_digitos}f}","). Error: ",abs(solution(xf)-yn1[n]))
print("Método Euler mejorado: (", f"{xn2[n]:.{n_digitos}f}",";", f"{yn2[n]:.{n_digitos}f}","). Error: ",abs(solution(xf)-yn2[n]))
print("Runge-Kutta: (", f"{xn3[n]:.{n_digitos}f}",";", f"{yn3[n]:.{n_digitos}f}","). Error: ",abs(solution(xf)-yn3[n]))
print("Solucion: ", solution(xf))
graficate(xn1, yn1, xn2, yn2, xn3, yn3, xs, ys)
