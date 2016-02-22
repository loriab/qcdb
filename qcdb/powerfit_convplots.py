data = [(1.1,0.9),(2.1,2),(3, 2.9)]


from __future__ import print_function
import sys
sys.path.append('/Users/loriab/linux/qcdb')
sys.path.append('/Users/loriab/linux/qcdb/databases')
import qcdb

a24 = qcdb.Database('A24')
a24.load_dilabio()
a24.load_f12dilabio()

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
mcs1 = ['CCSDT-CP-aqz','CCSDT-CP-a5z','CCSDT-CP-a6z']
rxn = 1
data1 = ([(mc,a24.hrxn[rxn].data[mc].value) for mc in mcs1])

N = len(mcs1)
x = np.arange(4, N+4)
y = [t for (s, t) in data1]

plt.plot(x, y, 'ro-')
def powerfit(x, a, b):
    return a*(x**b)
popt, pcov = curve_fit(powerfit, x, y)
plt.plot(x,powerfit(x,*popt),'b--')
plt.show()
print(popt)
print(pcov)
