import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import quad
from time import time

lims1 = (14.9, 4534.453)
def func1(x):
    f = ((x % 1234)/100)**2
    f += np.sqrt(x)*0.2
    f += np.sin(x/58 * np.exp(-x/2000)) * 30
    return f

lims2 = (-1e-6, 3e-6)
def func2(x):
    f = 1e-5 + x**2
    return f

lims3 = (474564, 474599)
def func3(x):
    f = 1/(1+(x%12))
    f += np.floor((x%23) / 5) * 0.3
    return f

lims4 = (-50, 50)
def func4(x):
    f = np.sin(x) * x
    return f


funcs = [func1, func2, func3, func4]
lims = [lims1, lims2, lims3, lims4]

start = time()

for func, lim in zip(funcs, lims):
    integral, error = quad(func, *lim, epsrel=1e-13, limit=10000, epsabs=0)
    print(integral, error)

print(f"Elapsed {(time() - start) * 1e3:0.0f} ms")

for func, lim in zip(funcs, lims):
    x = np.linspace(*lim, 500)
    plt.plot(x, func(x))
    plt.show()
