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

mcs1 = ['CCSDT-CP-adz', 'CCSDT-CP-atz', 'CCSDT-CP-aqz','CCSDT-CP-a5z','CCSDT-CP-a6z']
mcs2 = ['CCSDTAF12-CP-adz','CCSDTAF12-CP-atz','CCSDTAF12-CP-aqz','CCSDTAF12-CP-a5z']
mcs3 = ['CCSDTBF12-CP-adz','CCSDTBF12-CP-atz','CCSDTBF12-CP-aqz','CCSDTBF12-CP-a5z']
mcs4 = ['CCSDTCF12-CP-adz','CCSDTCF12-CP-atz','CCSDTCF12-CP-aqz']
mcs2f = ['CCSDTAF12-CP-dzf12','CCSDTAF12-CP-tzf12','CCSDTAF12-CP-qzf12']
mcs3f = ['CCSDTBF12-CP-dzf12','CCSDTBF12-CP-tzf12','CCSDTBF12-CP-qzf12']
mcs4f = ['CCSDTCF12-CP-dzf12','CCSDTCF12-CP-tzf12','CCSDTCF12-CP-qzf12']
mcs5 = ['CCSDT-CP-adtz','CCSDT-CP-atqz','CCSDT-CP-aq5z','CCSDT-CP-a56z']
mcs6 = ['CCSDTAF12-CP-adtz','CCSDTAF12-CP-atqz','CCSDTAF12-CP-aq5z']
mcs7 = ['CCSDTBF12-CP-adtz','CCSDTBF12-CP-atqz','CCSDTBF12-CP-aq5z']
mcs8 = ['CCSDTCF12-CP-adtz','CCSDTCF12-CP-atqz']
mcs6f = ['CCSDTAF12-CP-dtzf12','CCSDTAF12-CP-tqzf12']
mcs7f = ['CCSDTBF12-CP-dtzf12','CCSDTBF12-CP-tqzf12']
mcs8f = ['CCSDTCF12-CP-dtzf12','CCSDTCF12-CP-tqzf12']

#rxn = 1
#data1 = ([(mc,a24.hrxn[rxn].data[mc].value) for mc in mcs1])
#data2 = ([(mc,a24.hrxn[rxn].data[mc].value) for mc in mcs2])
#data3 = ([(mc,a24.hrxn[rxn].data[mc].value) for mc in mcs3])
#data4 = ([(mc,a24.hrxn[rxn].data[mc].value) for mc in mcs4])
#data5 = ([(mc,a24.hrxn[rxn].data[mc].value) for mc in mcs5])
#data6 = ([(mc,a24.hrxn[rxn].data[mc].value) for mc in mcs6])
#data7 = ([(mc,a24.hrxn[rxn].data[mc].value) for mc in mcs7])
#data8 = ([(mc,a24.hrxn[rxn].data[mc].value) for mc in mcs8])

def plot_convergence_aMNZ(title, conv, f12a, f12b, f12c, xconv, xf12a, xf12b, xf12c):
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
    plot_convergence_aMNZ(title, data1, data2, data3, data4, data5, data6, data7, data8)
    print(data1)
    print(data2)
    print(data3)
    print(data4)
    print(data5)
    print(data6)
    print(data7)
    print(data8)

def plot_convergence_MNZF12(title, conv, f12a, f12b, f12c, xconv, xf12a, xf12b, xf12c):
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
    
    labels = ['DZ-F12','TZ-F12','QZ-F12','a5Z','a6Z','DTZ-F12','TQZ-F12','aQ5Z','a56Z']
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
    title = 'A24-'+str(rxn)+'/XZ-F12'
    bench = a24.hrxn[rxn].data[a24.hrxn[rxn].benchmark].value
    data1 = []
    for mc in mcs1:
        try: data1.append((mc,a24.hrxn[rxn].data[mc].value))
        except: KeyError
        pass
    data2 = []
    for mc in mcs2f:
        try: data2.append((mc,a24.hrxn[rxn].data[mc].value))
        except: KeyError
        pass
    data3 = []
    for mc in mcs3f:
        try: data3.append((mc,a24.hrxn[rxn].data[mc].value))
        except: KeyError
        pass
    data4 = []
    for mc in mcs4f:
        try: data4.append((mc,a24.hrxn[rxn].data[mc].value))
        except: KeyError
        pass
    data5 = []
    for mc in mcs5:
        try: data5.append((mc,a24.hrxn[rxn].data[mc].value))
        except: KeyError
        pass
    data6 = []
    for mc in mcs6f:
        try: data6.append((mc,a24.hrxn[rxn].data[mc].value))
        except: KeyError
        pass
    data7 = []
    for mc in mcs7f:
        try: data7.append((mc,a24.hrxn[rxn].data[mc].value))
        except: KeyError
        pass
    data8 = []
    for mc in mcs8f:
        try: data8.append((mc,a24.hrxn[rxn].data[mc].value))
        except: KeyError
        pass
    plot_convergence_aMNZ(title, data1, data2, data3, data4, data5, data6, data7, data8)
    print(data1)
    print(data2)
    print(data3)
    print(data4)
    print(data5)
    print(data6)
    print(data7)
    print(data8)


#def plot_convergence(title, conv, f12a, f12b, f12c, f12aF, f12bF, f12cF, xconv, xf12a, xf12b, xf12c, xf12aF, xf12bF, xf12cF):
#    N1 = len(conv)
#    N2 = len(f12a)
#    N3 = len(f12b)
#    N4 = len(f12c)
#    N2F = len(f12aF)
#    N3F = len(f12bF)
#    N4F = len(f12cF)
#    N5 = len(xconv)
#    N6 = len(xf12a)
#    N7 = len(xf12b)
#    N8 = len(xf12c)
#    N6F = len(xf12aF)
#    N7F = len(xf12bF)
#    N8F = len(xf12cF)
#    shift = 6
#    x1 = np.arange(1, N1+1)
#    x2 = np.arange(1, N2+1)
#    x3 = np.arange(1, N3+1)
#    x4 = np.arange(1, N4+1)
#    x2F = np.arange(1, N2F+1)
#    x3F = np.arange(1, N3F+1)
#    x4F = np.arange(1, N4F+1)
#    x5 = np.arange(shift, N5+shift)
#    x6 = np.arange(shift, N6+shift)
#    x7 = np.arange(shift, N7+shift)
#    x8 = np.arange(shift, N8+shift)
#    x6F = np.arange(shift, N6F+shift)
#    x7F = np.arange(shift, N7F+shift)
#    x8F = np.arange(shift, N8F+shift)
#    y1 = [-num for (s, num) in data1]
#    y2 = [-num for (s, num) in data2]
#    y3 = [-num for (s, num) in data3]
#    y4 = [-num for (s, num) in data4]
#    y2F = [-num for (s, num) in data2F]
#    y3F = [-num for (s, num) in data3F]
#    y4F = [-num for (s, num) in data4F]
#    y5 = [-num for (s, num) in data5]
#    y6 = [-num for (s, num) in data6]
#    y7 = [-num for (s, num) in data7]
#    y8 = [-num for (s, num) in data8]
#    y6F = [-num for (s, num) in data6F]
#    y7F = [-num for (s, num) in data7F]
#    y8F = [-num for (s, num) in data8F]
#    ax, fig = plt.subplots(figsize=(12,6))
#    plt.plot(x1, y1, 'ro-',label = 'CCSD(T)/aXZ')
#    plt.plot(x2, y2, 'bs-',label = 'CCSD(T**)-F12a/aXZ')
#    plt.plot(x3, y3, 'g^-',label = 'CCSD(T**)-F12b/aXZ')
#    plt.plot(x4, y4, 'm*-',label = 'CCSD(T**)-F12c/aXZ')
#    plt.plot(x2F, y2F, 'b+-',label = 'CCSD(T**)-F12a/XZ-F12')
#    plt.plot(x3F, y3F, 'gx-',label = 'CCSD(T**)-F12b/XZ-F12')
#    plt.plot(x4F, y4F, 'm.-',label = 'CCSD(T**)-F12c/XZ-F12')
#    plt.plot(x5, y5, 'ro-')
#    plt.plot(x6, y6, 'bs-')
#    plt.plot(x7, y7, 'g^-')
#    plt.plot(x8, y8, 'm*-')
#    plt.plot(x6F, y6F, 'b+-')
#    plt.plot(x7F, y7F, 'gx-')
#    plt.plot(x8F, y8F, 'm.-')
#    
#    labels = ['aDZ','aTZ','aQZ','a5Z','a6Z','aDTZ','aTQZ','aQ5Z','a56Z']
#    xlabel = np.append(x1, x5)
#    plt.axhline(-bench, color='black',linestyle='--')
#    plt.title('A24-'+str(rxn))
#    plt.xlabel('Level of Theory')
#    plt.ylabel('Interaction Energy (kcal)')
#    plt.legend(loc = 'bottom right')
#    plt.xticks(xlabel, labels)
#    plt.margins(0.2)
#    plt.show()
#
#for rxn in [1,2,3]:
#    title = 'A24-'+str(rxn)
#    bench = a24.hrxn[rxn].data[a24.hrxn[rxn].benchmark].value
#    data1 = []
#    for mc in mcs1:
#        try: data1.append((mc,a24.hrxn[rxn].data[mc].value))
#        except: KeyError
#        pass
#    data2 = []
#    for mc in mcs2:
#        try: data2.append((mc,a24.hrxn[rxn].data[mc].value))
#        except: KeyError
#        pass
#    data3 = []
#    for mc in mcs3:
#        try: data3.append((mc,a24.hrxn[rxn].data[mc].value))
#        except: KeyError
#        pass
#    data4 = []
#    for mc in mcs4:
#        try: data4.append((mc,a24.hrxn[rxn].data[mc].value))
#        except: KeyError
#        pass
#    data2F = []
#    for mc in mcs2f:
#        try: data2.append((mc,a24.hrxn[rxn].data[mc].value))
#        except: KeyError
#        pass
#    data3F = []
#    for mc in mcs3f:
#        try: data3.append((mc,a24.hrxn[rxn].data[mc].value))
#        except: KeyError
#        pass
#    data4F = []
#    for mc in mcs4f:
#        try: data4.append((mc,a24.hrxn[rxn].data[mc].value))
#        except: KeyError
#        pass
#    data5 = []
#    for mc in mcs5:
#        try: data5.append((mc,a24.hrxn[rxn].data[mc].value))
#        except: KeyError
#        pass
#    data6 = []
#    for mc in mcs6:
#        try: data6.append((mc,a24.hrxn[rxn].data[mc].value))
#        except: KeyError
#        pass
#    data7 = []
#    for mc in mcs7:
#        try: data7.append((mc,a24.hrxn[rxn].data[mc].value))
#        except: KeyError
#        pass
#    data8 = []
#    for mc in mcs8:
#        try: data8.append((mc,a24.hrxn[rxn].data[mc].value))
#        except: KeyError
#        pass
#    data6F = []
#    for mc in mcs6f:
#        try: data6.append((mc,a24.hrxn[rxn].data[mc].value))
#        except: KeyError
#        pass
#    data7F = []
#    for mc in mcs7f:
#        try: data7.append((mc,a24.hrxn[rxn].data[mc].value))
#        except: KeyError
#        pass
#    data8F = []
#    for mc in mcs8f:
#        try: data8.append((mc,a24.hrxn[rxn].data[mc].value))
#        except: KeyError
#        pass
#    plot_convergence(title, data1, data2, data3, data4, data2F, data3F, data4F, data5, data6, data7, data8, data6F, data7F, data8F)
#    print data1
#    print data2
#    print data3
#    print data4
#    print data2F
#    print data3F
#    print data4F
#    print data5
#    print data6
#    print data7
#    print data8
#    print data6F
#    print data7F
#    print data8F
