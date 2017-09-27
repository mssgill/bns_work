import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

from scipy.constants import h,k,c

def func(wa, T):
    #print wa
    lam = 1e-6 * wa
    print lam
    return  2*h*c**2 / (lam**5 * (np.exp(h*c / (lam*k*T)) - 1))

wa = np.linspace(0.1, 2, 100)
xdata = wa

#define the data to be fit with some noise


y = func(wa,  5000)
ydata = y # + y_noise

plt.plot(xdata, ydata, 'b-', label='data')

#Fit for the parameters a, b, c of the function func

print "About to fit"
popt, pcov = curve_fit(func, xdata, ydata)

print " popt = ", popt
print "After fit"
plt.plot(xdata, func(xdata, *popt), 'r-', label='fit')

#Constrain the optimization to the region of 0 < a < 3, 0 < b < 2 and 0 < c < 1:

popt, pcov = curve_fit(func, xdata, ydata)
plt.plot(xdata, func(xdata, *popt), 'g--')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()

plt.show()

print popt
print math.sqrt(pcov[0,0])
