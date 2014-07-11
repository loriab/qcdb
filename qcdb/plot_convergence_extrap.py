import sys
sys.path.append('/Users/loriab/linux/qcdb')
sys.path.append('/Users/loriab/linux/qcdb/databases')
import qcdb

a24 = qcdb.Database('A24')
a24.load_dilabio()
a24.load_f12dilabio()

import matplotlib.pyplot as plt
import numpy as np

mcs1 = ['CCSDT-CP-adz', 'CCSDT-CP-atz', 'CCSDT-CP-aqz','CCSDT-CP-a5z','CCSDT-CP-a6z']
mcs2 = ['CCSDTAF12-CP-adz','CCSDTAF12-CP-atz','CCSDTAF12-CP-aqz','CCSDTAF12-CP-a5z']
mcs3 = ['CCSDTBF12-CP-adz','CCSDTBF12-CP-atz','CCSDTBF12-CP-aqz','CCSDTBF12-CP-a5z']
mcs4 = ['CCSDTCF12-CP-adz','CCSDTBF12-CP-atz','CCSDTCF12-CP-aqz']
mcs5 = ['CCSDT-CP-adtz','CCSDT-CP-atqz','CCSDT-CP-aq5z','CCSDT-CP-a56z']
mcs6 = ['CCSDTAF12-CP-adtz','CCSDTAF12-CP-atqz','CCSDTAF12-CP-aq5z']
mcs7 = ['CCSDTBF12-CP-adtz','CCSDTBF12-CP-atqz','CCSDTBF12-CP-aq5z']
mcs8 = ['CCSDTCF12-CP-adtz','CCSDTCF12-CP-atqz']

#rxn = 1
#data1 = ([(mc,a24.hrxn[rxn].data[mc].value) for mc in mcs1])
#data2 = ([(mc,a24.hrxn[rxn].data[mc].value) for mc in mcs2])
#data3 = ([(mc,a24.hrxn[rxn].data[mc].value) for mc in mcs3])
#data4 = ([(mc,a24.hrxn[rxn].data[mc].value) for mc in mcs4])
#data5 = ([(mc,a24.hrxn[rxn].data[mc].value) for mc in mcs5])
#data6 = ([(mc,a24.hrxn[rxn].data[mc].value) for mc in mcs6])
#data7 = ([(mc,a24.hrxn[rxn].data[mc].value) for mc in mcs7])
#data8 = ([(mc,a24.hrxn[rxn].data[mc].value) for mc in mcs8])
def plot_convergence_extrap(title, conv, f12a, f12b, f12c, xconv, xf12a, xf12b, xf12c):
    N1 = len(conv)
    N2 = len(f12a)
    N3 = len(f12b)
    N4 = len(f12c)
    N5 = len(xconv)
    N6 = len(xf12a)
    N7 = len(xf12b)
    N8 = len(xf12c)
    shift = 6
    x1 = np.arange(1, N1+1)
    x2 = np.arange(1, N2+1)
    x3 = np.arange(1, N3+1)
    x4 = np.arange(1, N4+1)
    x5 = np.arange(shift, N5+shift)
    x6 = np.arange(shift, N6+shift)
    x7 = np.arange(shift, N7+shift)
    x8 = np.arange(shift, N8+shift)
    y1 = [-num for (s, num) in data1]
    y2 = [-num for (s, num) in data2]
    y3 = [-num for (s, num) in data3]
    y4 = [-num for (s, num) in data4]
    y5 = [-num for (s, num) in data5]
    y6 = [-num for (s, num) in data6]
    y7 = [-num for (s, num) in data7]
    y8 = [-num for (s, num) in data8]
    ax, fig = plt.subplots(figsize=(12,6))
    plt.plot(x1, y1, 'ro-',label = 'CCSD(T)')
    plt.plot(x2, y2, 'bs-',label = 'CCSD(T**)-F12a')
    plt.plot(x3, y3, 'g^-',label = 'CCSD(T**)-F12b')
    plt.plot(x4, y4, 'm*-',label = 'CCSD(T**)-F12c')
    plt.plot(x5, y5, 'ro-')
    plt.plot(x6, y6, 'bs-')
    plt.plot(x7, y7, 'g^-')
    plt.plot(x8, y8, 'm*-')
    
    labels = ['aDZ','aTZ','aQZ','a5Z','a6Z','aDTZ','aTQZ','aQ5Z','a56Z']
    xlabel = np.append(x1, x5)
    plt.axhline(-bench, color='black',linestyle='--')
    plt.title('A24-'+str(rxn))
    plt.xlabel('Level of Theory')
    plt.ylabel('Interaction Energy (kcal)')
    plt.legend(loc = 'bottom right')
    plt.xticks(xlabel, labels)
    plt.margins(0.2)
    plt.show()
    
for rxn in [1,2,3]:
    title = 'A24-'+str(rxn)
    bench = a24.hrxn[rxn].data[a24.hrxn[rxn].benchmark].value
    data1 = []
    for mc in mcs1:
        try: data1.append((mc,a24.hrxn[rxn].data[mc].value))
        except: KeyError
        pass
    data2 = []
    for mc in mcs2:
        try: data2.append((mc,a24.hrxn[rxn].data[mc].value))
        except: KeyError
        pass
    data3 = []
    for mc in mcs3:
        try: data3.append((mc,a24.hrxn[rxn].data[mc].value))
        except: KeyError
        pass
    data4 = []
    for mc in mcs4:
        try: data4.append((mc,a24.hrxn[rxn].data[mc].value))
        except: KeyError
        pass
    data5 = []
    for mc in mcs5:
        try: data5.append((mc,a24.hrxn[rxn].data[mc].value))
        except: KeyError
        pass
    data6 = []
    for mc in mcs6:
        try: data6.append((mc,a24.hrxn[rxn].data[mc].value))
        except: KeyError
        pass
    data7 = []
    for mc in mcs7:
        try: data7.append((mc,a24.hrxn[rxn].data[mc].value))
        except: KeyError
        pass
    data8 = []
    for mc in mcs8:
        try: data8.append((mc,a24.hrxn[rxn].data[mc].value))
        except: KeyError
        pass
    plot_convergence_extrap(title, data1, data2, data3, data4, data5, data6, data7, data8)

print data1
print data2
print data3
print data4
print data5
print data6
print data7
print data8
