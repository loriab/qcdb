from __future__ import print_function
import os
import sys
import glob
import math
import collections
sys.path.append('C:\Users\Owner\Documents\GitHub\qcdb')
sys.path.append('C:\Users\Owner\Documents\GitHub\qcdb\databases')
import qcdb
from qcdb.psivarrosetta import useme2psivar
from qcdb.modelchems import Method, BasisSet, methods, bases
import pandas as pd
import numpy as np
pd.set_option('display.max_columns', 10)
pd.set_option('display.width', 200)


# <<< read usemefiles and convert to giant DataFrame >>>

dbse = 'A24'
dbobj = qcdb.Database(dbse)
path = r"""C:\Users\Owner\Documents\f12dilabiousemefiles\usemefiles"""
h2kc = qcdb.psi_hartree2kcalmol
    
rawdata = collections.defaultdict(dict)
for useme in glob.glob('%s/%s*useme*' % (path, dbse)):
    spl = os.path.basename(useme).split('.')
    #dbse = spl[0].split('-')[0]
    #ocalc = spl[0].split('-')[1]
    #optns = '_'.join(spl[0].split('-')[2:-1])
    basis = spl[0].split('-')[-1]
    piece = '.'.join(spl[1:])
    print(useme, basis, piece)
    tmp = pd.read_csv('%s' % (useme), index_col=0, sep='\s+', comment='#', na_values='None', names=['rxn', 'dimer', 'monoA-CP', 'monoB-CP'])
    print(useme, basis, piece)
    print(tmp.head(10))
    rawdata[basis][useme2psivar[piece]] = tmp.dropna(how='all')

datazip = {}
for key, val in rawdata.iteritems():
    datazip[key] = pd.concat(val)
df = pd.concat(datazip)
df.index.names = ['bstrt', 'psivar', 'rxn']


# <<< define utility functions >>>

def ie(rgts):
    return rgts['dimer'] - rgts['monoA-CP'] - rgts['monoB-CP']

def append_result_of_func_with_funcargs_to_master_DataFrame_atlevel_with_label(master, label, atlevel, func, funcargs):
    odr = {'psivar': [0, 2, 1],
           'bstrt': [2, 0, 1]}        
    data_rich_args = [master.xs(pv, level='psivar') if isinstance(pv, basestring) else pv for pv in funcargs]
    temp = func(data_rich_args)
    temp[atlevel] = label
    temp.set_index(atlevel, append=True, inplace=True)
    temp = temp.reorder_levels(odr[atlevel])
    return master.append(temp, verify_integrity=True)


# <<< define simple functions relating psi variables; follow args use in omega >>>

def xtpl_power(args):
    power, zHI, eHI, eLO = args
    return (eHI * zHI ** power - eLO * (zHI - 1) ** power) / (zHI ** power - (zHI - 1) ** power)

#def corl_xtpl_helgaker_2(zHI, eHI, eLO):
#    return xtpl_power(power=3.0, zHI=zHI, eHI=eHI, eLO=eLO)

def dispersion_weighting(args):
    omega, hb_mtd, dd_mtd = args
    return omega * hb_mtd + (1.0 - omega) * dd_mtd

def omega(args):
    alpha, beta, ratio = args
    return 0.5 * (1.0 + np.tanh(alpha + beta * ratio))
        
        
# <<< append to main DataFrame basic psivar equalities not explicit to useme structure >>>

lvl = 'psivar'
pv0 = collections.OrderedDict()
pv0['MP2-F12 TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'MP2-F12 CORRELATION ENERGY']}
pv0['CCSD-F12A TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD-F12A CORRELATION ENERGY']}
pv0['CCSD-F12B TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD-F12B CORRELATION ENERGY']}
pv0['CCSD-F12C TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD-F12C CORRELATION ENERGY']}
pv0['CCSD(T**)-F12A TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T**)-F12A CORRELATION ENERGY']}
pv0['CCSD(T**)-F12B TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T**)-F12B CORRELATION ENERGY']}
pv0['CCSD(T**)-F12C TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T**)-F12C CORRELATION ENERGY']}
pv0['CCSD(T)-F12A CORRELATION ENERGY'] = {'func': sum, 'args': ['CCSD-F12A CORRELATION ENERGY', '(T)-F12AB CORRECTION ENERGY']}
pv0['CCSD(T)-F12B CORRELATION ENERGY'] = {'func': sum, 'args': ['CCSD-F12B CORRELATION ENERGY', '(T)-F12AB CORRECTION ENERGY']}
pv0['CCSD(T)-F12C CORRELATION ENERGY'] = {'func': sum, 'args': ['CCSD-F12C CORRELATION ENERGY', '(T)-F12C CORRECTION ENERGY']}
pv0['CCSD(T)-F12A TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T)-F12A CORRELATION ENERGY']}
pv0['CCSD(T)-F12B TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T)-F12B CORRELATION ENERGY']}
pv0['CCSD(T)-F12C TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'CCSD(T)-F12C CORRELATION ENERGY']}
for pvar, action in pv0.iteritems():
    try:
        df = append_result_of_func_with_funcargs_to_master_DataFrame_atlevel_with_label(
            master=df, label=pvar, atlevel=lvl, func=action['func'], funcargs=action['args'])
    except KeyError as e:
        print("""Building df pv0: empty index '%s' because missing %s""" % (pvar, e))


# <<< miscellaneous pre-computing >>>

ratio_HFCABS_MP2F12 = ie(df.xs('HF-CABS TOTAL ENERGY', level='psivar')) / ie(df.xs('MP2-F12 TOTAL ENERGY', level='psivar'))
df_omega = pd.concat({
    'adz': omega([-1, 4, ratio_HFCABS_MP2F12['adz']]),
    'atz': omega([0.4, 0.6, ratio_HFCABS_MP2F12['atz']]),
    })
df_omega = pd.DataFrame(df_omega, columns=['dimer'])
df_omega['psivar'] = 'DW-CCSD(T)-F12 OMEGA'
df_omega['monoA-CP'] = df_omega['dimer']
df_omega['monoB-CP'] = df_omega['dimer']
df_omega.set_index('psivar', append=True, inplace=True)
df_omega = df_omega.reorder_levels([0, 2, 1])
df_omega.index.names = ['bstrt', 'psivar', 'rxn']


# <<< append to main DataFrame computable method quantities >>>

lvl = 'psivar'
pv1 = collections.OrderedDict()
pv1['DW-CCSD(T**)-F12 CORRELATION ENERGY'] = {'func': dispersion_weighting, 'args': [df_omega.xs('DW-CCSD(T)-F12 OMEGA', level=lvl), 'CCSD(T**)-F12A CORRELATION ENERGY', 'CCSD(T**)-F12B CORRELATION ENERGY']}
pv1['DW-CCSD(T**)-F12 TOTAL ENERGY'] = {'func': sum, 'args': ['HF-CABS TOTAL ENERGY', 'DW-CCSD(T**)-F12 CORRELATION ENERGY']}
for pvar, action in pv1.iteritems():
    df = append_result_of_func_with_funcargs_to_master_DataFrame_atlevel_with_label(
        master=df, label=pvar, atlevel=lvl, func=action['func'], funcargs=action['args'])

# <<< define extra Method and BasisSet objects >>>

# Note dict key must match object name; i.e., first two strings on line must be the same

bases['Hill1_adtz'] = BasisSet('Hill1_adtz', build=[['HillCC_adtz'], ['atz', 'HillCC_adtz']])  # TODO should have None or non-xtpl first element?
bases['Hill1_atqz'] = BasisSet('Hill1_atqz', build=[['HillCC_atqz'], ['aqz', 'HillCC_atqz']])
bases['Hill1_aq5z'] = BasisSet('Hill1_aq5z', build=[['HillCC_aq5z'], ['a5z', 'HillCC_aq5z']])
bases['Hill1_dtzf12'] = BasisSet('Hill1_dtzf12', build=[['HillCC_dtzf12'], ['tzf12', 'HillCC_dtzf12']])
bases['Hill1_tqzf12'] = BasisSet('Hill1_tqzf12', build=[['HillCC_tqzf12'], ['qzf12', 'HillCC_tqzf12']])
methods['CCSDTNSAF12'] = Method('CCSDTNSAF12', fullname='CCSD(T)-F12a')
methods['CCSDTNSBF12'] = Method('CCSDTNSBF12', fullname='CCSD(T)-F12b')
methods['CCSDTNSCF12'] = Method('CCSDTNSCF12', fullname='CCSD(T)-F12c')
bases['Hill2_dtzf12'] = BasisSet('Hill2_dtzf12', build=[None, None, ['tzf12', 'HillCC_dtzf12', 'HillT_dtzf12']])
bases['Hill2_tqzf12'] = BasisSet('Hill2_tqzf12', build=[None, None, ['qzf12', 'HillCC_tqzf12', 'HillT_tqzf12']])
bases['Hill2_adtz'] = BasisSet('Hill2_adtz', build=[None, None, ['atz', 'HillCC_adtz', 'HillT_adtz']])
bases['Hill2_atqz'] = BasisSet('Hill2_atqz', build=[None, None, ['aqz', 'HillCC_atqz', 'HillT_atqz']])
bases['Hill2_aq5z'] = BasisSet('Hill2_aq5z', build=[None, None, ['a5z', 'HillCC_aq5z', 'HillT_aq5z']])

# <<< append to main DataFrame computable basis treatment quantities >>>

lvl = 'bstrt'
pv2 = collections.OrderedDict()
pv2['adtz'] = {'func': xtpl_power, 'args': [3.0, 3, df.xs('atz', level=lvl), df.xs('adz', level=lvl)]}
pv2['atqz'] = {'func': xtpl_power, 'args': [3.0, 4, df.xs('aqz', level=lvl), df.xs('atz', level=lvl)]}
pv2['aq5z'] = {'func': xtpl_power, 'args': [3.0, 5, df.xs('a5z', level=lvl), df.xs('aqz', level=lvl)]}
pv2['dtzf12'] = {'func': xtpl_power, 'args': [3.0, 3, df.xs('tzf12', level=lvl), df.xs('dzf12', level=lvl)]}
pv2['tqzf12'] = {'func': xtpl_power, 'args': [3.0, 4, df.xs('qzf12', level=lvl), df.xs('tzf12', level=lvl)]}
# Hill xtpl for CCSD-F12b from Table X of JCP 131 194105 (2009)
# TODO should only be applied to CCF12, as per definition (literally, only CCF12B)
pv2['HillCC_adtz'] = {'func': xtpl_power, 'args': [2.483070, 3, df.xs('atz', level=lvl), df.xs('adz', level=lvl)]}
pv2['HillCC_atqz'] = {'func': xtpl_power, 'args': [4.255221, 4, df.xs('aqz', level=lvl), df.xs('atz', level=lvl)]}
pv2['HillCC_aq5z'] = {'func': xtpl_power, 'args': [4.910269, 5, df.xs('a5z', level=lvl), df.xs('aqz', level=lvl)]}
pv2['HillCC_dtzf12'] = {'func': xtpl_power, 'args': [3.144518, 3, df.xs('tzf12', level=lvl), df.xs('dzf12', level=lvl)]}
pv2['HillCC_tqzf12'] = {'func': xtpl_power, 'args': [4.595995, 4, df.xs('qzf12', level=lvl), df.xs('tzf12', level=lvl)]}
# Hill xtpl for unscaled (T)-F12 from Table XI of JCP 131 194105 (2009)
pv2['HillT_adtz'] = {'func': xtpl_power, 'args': [2.790300, 3, df.xs('atz', level=lvl), df.xs('adz', level=lvl)]}
pv2['HillT_dtzf12'] = {'func': xtpl_power, 'args': [2.615472, 3, df.xs('tzf12', level=lvl), df.xs('dzf12', level=lvl)]}
for pvar, action in pv2.iteritems():
    df = append_result_of_func_with_funcargs_to_master_DataFrame_atlevel_with_label(
        master=df, label=pvar, atlevel=lvl, func=action['func'], funcargs=action['args'])


# <<< define simple functions codifying cbs() piecing >>>

# TODO move these to be functions in Method or BasisSet

# min and max pieces for which basis treatment intent is clearly understood
# with the understanding that Method should be defining all those pieces that have uniform basis
#           desi max_mtd    min_bas    max_bas
#HF/adz     1    1          1           inf (1)
#MP2/adz    1    2          1           inf (1)
#CC/adz     1    3          1           inf (1)
#HF/adtz    1    1          1           2
#MP2/adtz   2    2          1           2
#CC/adtz    2    3          1           2
#CC/atqzadz 3    3          2           3 
#   stages = min(max_mtd, max_bas)

def compute_max_bas(bas):
    return len(bases[bas].build)
    
def generic_bas(Nstage, bas):
    return bases[bas].build[Nstage - 1]  # 1-indexed stage counting to 0-indexed array storage

def compute_max_mtd(stub):
    if stub in ['HFCABS', 'HF', 'SCF']:
        return 1
    elif 'MP2' in stub:
        return 2
    else:
        return 3

def generic_mtd(Nstage, stub):
    if Nstage == 1:
        return ['%s TOTAL ENERGY' % (stub)]
    elif Nstage == 2:
        return ['%s TOTAL ENERGY' % ('HF-CABS' if 'F12' in stub else 'SCF'), '%s CORRELATION ENERGY' % (stub)]
    elif Nstage == 3:
        return ['%s TOTAL ENERGY' % ('HF-CABS' if 'F12' in stub else 'SCF'), '%s CORRELATION ENERGY' % (stub), '%s CC CORRECTION ENERGY' % (stub)]

def build_from_lists(mtdlist, baslist):
    return ie(sum([df.loc[bas].loc[pcs] for pcs, bas in zip(mtdlist, baslist)]))

def build(method, basis):
    Nstage = min(compute_max_mtd(method), compute_max_bas(basis))
    baslist = generic_bas(Nstage, basis)
    mtdlist = generic_mtd(Nstage, methods[method].fullname.upper())
        #print '\n <<<', methods[method].fullname, '/', bases[basis].fullname, '>>>'
        #print 'stages:', 'M:', compute_max_mtd(method), 'B:', compute_max_bas(basis), 'U:', Nstage
        #print 'pcss:', mtdlist
        #print 'bass:', baslist
    return ie(sum([df.loc[bas].loc[pcs] for pcs, bas in zip(mtdlist, baslist)]))


# <<< assemble all model chemistries into columns of new DataFrame >>>

rxns = ['%s-%s' % (dbse, rxn) for rxn in dbobj.hrxn.keys()]
mine = pd.DataFrame({}, index=rxns)
mine.index.names = ['rxn']

for mtd in methods:
    if mtd == 'HFCABS' or 'F12' in mtd:
        if 'SC' not in mtd and 'MP2C' not in mtd and 'DWMP2' not in mtd:
            for bas in ['adz', 'atz', 'aqz', 'a5z', 'dzf12', 'tzf12', 'qzf12', 
                'adtz', 'atqz', 'aq5z', 'dtzf12', 'tqzf12', 
                'Hill1_adtz', 'Hill1_atqz', 'Hill1_aq5z', 'Hill1_dtzf12', 'Hill1_tqzf12']:
                mine['%s-%s-%s' % (mtd, 'CP', bas)] = build(mtd, bas)
                
#mine['CCSDTNSBF12-CP-Hill2_dtzf12'] = build_from_lists(['HF-CABS TOTAL ENERGY', 'CCSD-F12B CORRELATION ENERGY', '(T)-F12AB CORRECTION ENERGY'], ['tzf12', 'HillCC_dtzf12', 'HillT_dtzf12'])
mine['CCSDTNSBF12-CP-Hill2_adtz'] = build_from_lists(['HF-CABS TOTAL ENERGY', 'CCSD-F12B CORRELATION ENERGY', '(T)-F12AB CORRECTION ENERGY'], ['atz', 'HillCC_adtz', 'HillT_adtz'])

print('scf', h2kc * ie(df.loc['atz'].loc['HF-CABS TOTAL ENERGY'].loc['A24-9']))
print('cc-f12b', h2kc * ie(df.loc['HillCC_adtz'].loc['CCSD-F12B CORRELATION ENERGY'].loc['A24-9']))
print('(t)', h2kc * ie(df.loc['HillT_adtz'].loc['(T)-F12AB CORRECTION ENERGY'].loc['A24-9']))

# <<< test cases >>>

try:
    qcdb.compare_values(-2.4749, h2kc * mine['HFCABS-CP-adz'            ]['A24-9'], 4, 'HFCABS-CP-adz')
    qcdb.compare_values(-4.5972, h2kc * mine['CCSDTAF12-CP-adz'         ]['A24-9'], 4, 'CCSDTAF12-CP-adz')
    qcdb.compare_values(-4.4354, h2kc * mine['CCSDTBF12-CP-adz'         ]['A24-9'], 4, 'CCSDTBF12-CP-adz')
    qcdb.compare_values(-4.4726, h2kc * mine['CCSDTCF12-CP-adz'         ]['A24-9'], 4, 'CCSDTCF12-CP-adz')  # corr, translate error b/c tz scf
    qcdb.compare_values(-4.5830, h2kc * mine['DWCCSDTF12-CP-adz'        ]['A24-9'], 4, 'DWCCSDTF12-CP-adz')
    qcdb.compare_values(-4.5913, h2kc * mine['DWCCSDTF12-CP-atz'        ]['A24-9'], 4, 'DWCCSDTF12-CP-atz')  # corr, not separate atz def in reapsets
    qcdb.compare_values(-4.5669, h2kc * mine['MP2F12-CP-adz'            ]['A24-9'], 4, 'MP2F12-CP-adz')
    qcdb.compare_values(-4.5767, h2kc * mine['CCSDTAF12-CP-a5z'         ]['A24-9'], 4, 'CCSDTAF12-CP-a5z')
    qcdb.compare_values(-4.5855, h2kc * mine['CCSDTAF12-CP-aqz'         ]['A24-9'], 4, 'CCSDTAF12-CP-aqz')
    qcdb.compare_values(-4.6042, h2kc * mine['CCSDTAF12-CP-atz'         ]['A24-9'], 4, 'CCSDTAF12-CP-atz')
    qcdb.compare_values(-4.3559, h2kc * mine['CCSDTAF12-CP-dzf12'       ]['A24-9'], 4, 'CCSDTAF12-CP-dzf12')  # corr, translate error b/c tz scf
    qcdb.compare_values(-4.2422, h2kc * mine['CCSDTBF12-CP-dzf12'       ]['A24-9'], 4, 'CCSDTBF12-CP-dzf12')  # corr, translate error b/c tz scf
    qcdb.compare_values(-4.2744, h2kc * mine['CCSDTCF12-CP-dzf12'       ]['A24-9'], 4, 'CCSDTCF12-CP-dzf12')  # corr, translate error b/c tz scf
    qcdb.compare_values(-2.4582, h2kc * mine['HFCABS-CP-dzf12'          ]['A24-9'], 4, 'HFCABS-CP-dzf12')  # added
    #qcdb.compare_values(, h2kc * mine['HFCABS-CP-tzf12'          ]['A24-9'], 4, 'HFCABS-CP-tzf12')  # added
    #qcdb.compare_values(-2.4692, h2kc * mine['HFCABS-CP-qzf12'          ]['A24-9'], 4, 'HFCABS-CP-qzf12')  # added
    qcdb.compare_values(-4.5204, h2kc * mine['CCSDTAF12-CP-tzf12'       ]['A24-9'], 4, 'CCSDTAF12-CP-tzf12')
    qcdb.compare_values(-4.5692, h2kc * mine['CCSDTAF12-CP-qzf12'       ]['A24-9'], 4, 'CCSDTAF12-CP-qzf12')
    qcdb.compare_values(-4.5556, h2kc * mine['CCSDTBF12-CP-a5z'         ]['A24-9'], 4, 'CCSDTBF12-CP-a5z')
    qcdb.compare_values(-4.5495, h2kc * mine['CCSDTBF12-CP-aqz'         ]['A24-9'], 4, 'CCSDTBF12-CP-aqz')
    qcdb.compare_values(-4.5365, h2kc * mine['CCSDTBF12-CP-atz'         ]['A24-9'], 4, 'CCSDTBF12-CP-atz')
    qcdb.compare_values(-4.6084, h2kc * mine['CCSDTAF12-CP-adtz'        ]['A24-9'], 4, 'CCSDTAF12-CP-adtz')  # corr, working from older version of Dom's f12dilabio.py
    qcdb.compare_values(-4.5673, h2kc * mine['CCSDTAF12-CP-aq5z'        ]['A24-9'], 4, 'CCSDTAF12-CP-aq5z')  # corr, working from older version of Dom's f12dilabio.py
    qcdb.compare_values(-4.5705, h2kc * mine['CCSDTAF12-CP-atqz'        ]['A24-9'], 4, 'CCSDTAF12-CP-atqz')  # corr, working from older version of Dom's f12dilabio.py
    qcdb.compare_values(-4.5851, h2kc * mine['CCSDTAF12-CP-dtzf12'      ]['A24-9'], 4, 'CCSDTAF12-CP-dtzf12')  # corr, working from older version of Dom's f12dilabio.py
    qcdb.compare_values(-4.6017, h2kc * mine['CCSDTAF12-CP-tqzf12'      ]['A24-9'], 4, 'CCSDTAF12-CP-tqzf12')  # corr, working from older version of Dom's f12dilabio.py
    qcdb.compare_values(-4.5804, h2kc * mine['CCSDTBF12-CP-adtz'        ]['A24-9'], 4, 'CCSDTBF12-CP-adtz')  # corr, working from older version of Dom's f12dilabio.py
    qcdb.compare_values(-4.5619, h2kc * mine['CCSDTBF12-CP-aq5z'        ]['A24-9'], 4, 'CCSDTBF12-CP-aq5z')  # corr, working from older version of Dom's f12dilabio.py
    qcdb.compare_values(-4.5576, h2kc * mine['CCSDTBF12-CP-atqz'        ]['A24-9'], 4, 'CCSDTBF12-CP-atqz')  # corr, working from older version of Dom's f12dilabio.py
    qcdb.compare_values(-4.5547, h2kc * mine['CCSDTBF12-CP-dtzf12'      ]['A24-9'], 4, 'CCSDTBF12-CP-dtzf12')  # corr, working from older version of Dom's f12dilabio.py
    qcdb.compare_values(-4.5873, h2kc * mine['CCSDTBF12-CP-tqzf12'      ]['A24-9'], 4, 'CCSDTBF12-CP-tqzf12')  # corr, working from older version of Dom's f12dilabio.py 
    qcdb.compare_values(-4.5752, h2kc * mine['CCSDTCF12-CP-adtz'        ]['A24-9'], 4, 'CCSDTCF12-CP-adtz')  # corr, working from older version of Dom's f12dilabio.py
    qcdb.compare_values(-4.5447, h2kc * mine['CCSDTCF12-CP-atqz'        ]['A24-9'], 4, 'CCSDTCF12-CP-atqz')  # corr, working from older version of Dom's f12dilabio.py
    qcdb.compare_values(-4.5462, h2kc * mine['CCSDTCF12-CP-dtzf12'      ]['A24-9'], 4, 'CCSDTCF12-CP-dtzf12')  # corr, working from older version of Dom's f12dilabio.py
    qcdb.compare_values(-4.5756, h2kc * mine['CCSDTCF12-CP-tqzf12'      ]['A24-9'], 4, 'CCSDTCF12-CP-tqzf12')  # corr, working from older version of Dom's f12dilabio.py
    qcdb.compare_values(-4.5377, h2kc * mine['CCSDTBF12-CP-qzf12'       ]['A24-9'], 4, 'CCSDTBF12-CP-qzf12')
    qcdb.compare_values(-4.4653, h2kc * mine['CCSDTBF12-CP-tzf12'       ]['A24-9'], 4, 'CCSDTBF12-CP-tzf12')
    qcdb.compare_values(-4.5451, h2kc * mine['CCSDTCF12-CP-aqz'         ]['A24-9'], 4, 'CCSDTCF12-CP-aqz')
    qcdb.compare_values(-4.5439, h2kc * mine['CCSDTCF12-CP-atz'         ]['A24-9'], 4, 'CCSDTCF12-CP-atz')
    qcdb.compare_values(-4.5324, h2kc * mine['CCSDTCF12-CP-qzf12'       ]['A24-9'], 4, 'CCSDTCF12-CP-qzf12')
    qcdb.compare_values(-4.4688, h2kc * mine['CCSDTCF12-CP-tzf12'       ]['A24-9'], 4, 'CCSDTCF12-CP-tzf12')
    qcdb.compare_values(-2.4739, h2kc * mine['HFCABS-CP-a5z'            ]['A24-9'], 4, 'HFCABS-CP-a5z')
    qcdb.compare_values(-2.4737, h2kc * mine['HFCABS-CP-aqz'            ]['A24-9'], 4, 'HFCABS-CP-aqz')
    qcdb.compare_values(-2.4719, h2kc * mine['HFCABS-CP-atz'            ]['A24-9'], 4, 'HFCABS-CP-atz')
    qcdb.compare_values(-4.6070, h2kc * mine['MP2F12-CP-atz'            ]['A24-9'], 4, 'MP2F12-CP-atz')
    qcdb.compare_values(-4.5927, h2kc * mine['MP2F12-CP-aqz'            ]['A24-9'], 4, 'MP2F12-CP-aqz')
    qcdb.compare_values(-4.5954, h2kc * mine['MP2F12-CP-a5z'            ]['A24-9'], 4, 'MP2F12-CP-a5z')
    qcdb.compare_values(-4.6251, h2kc * mine['MP2F12-CP-adtz'           ]['A24-9'], 4, 'MP2F12-CP-adtz')  # corr, working from older version of Dom's f12dilabio.py
    qcdb.compare_values(-4.5980, h2kc * mine['MP2F12-CP-aq5z'           ]['A24-9'], 4, 'MP2F12-CP-aq5z')  # corr, working from older version of Dom's f12dilabio.py
    qcdb.compare_values(-4.5810, h2kc * mine['MP2F12-CP-atqz'           ]['A24-9'], 4, 'MP2F12-CP-atqz')  # corr, working from older version of Dom's f12dilabio.py
    qcdb.compare_values(-4.6099, h2kc * mine['CCSDTAF12-CP-Hill1_adtz'  ]['A24-9'], 4, 'CCSDTAF12-CP-Hill1_adtz')
    qcdb.compare_values(-4.5769, h2kc * mine['CCSDTAF12-CP-Hill1_atqz'  ]['A24-9'], 4, 'CCSDTAF12-CP-Hill1_atqz')
    qcdb.compare_values(-4.5722, h2kc * mine['CCSDTAF12-CP-Hill1_aq5z'  ]['A24-9'], 4, 'CCSDTAF12-CP-Hill1_aq5z')
    qcdb.compare_values(-4.5800, h2kc * mine['CCSDTAF12-CP-Hill1_dtzf12']['A24-9'], 4, 'CCSDTAF12-CP-Hill1_dtzf12')
    qcdb.compare_values(-4.5854, h2kc * mine['CCSDTAF12-CP-Hill1_tqzf12']['A24-9'], 4, 'CCSDTAF12-CP-Hill1_tqzf12')
    qcdb.compare_values(-4.5965, h2kc * mine['CCSDTBF12-CP-Hill1_adtz'  ]['A24-9'], 4, 'CCSDTBF12-CP-Hill1_adtz')
    qcdb.compare_values(-4.5541, h2kc * mine['CCSDTBF12-CP-Hill1_atqz'  ]['A24-9'], 4, 'CCSDTBF12-CP-Hill1_atqz')
    qcdb.compare_values(-4.5586, h2kc * mine['CCSDTBF12-CP-Hill1_aq5z'  ]['A24-9'], 4, 'CCSDTBF12-CP-Hill1_aq5z')
    qcdb.compare_values(-4.5476, h2kc * mine['CCSDTBF12-CP-Hill1_dtzf12']['A24-9'], 4, 'CCSDTBF12-CP-Hill1_dtzf12')
    qcdb.compare_values(-4.5624, h2kc * mine['CCSDTBF12-CP-Hill1_tqzf12']['A24-9'], 4, 'CCSDTBF12-CP-Hill1_tqzf12')
    qcdb.compare_values(-4.5867, h2kc * mine['CCSDTCF12-CP-Hill1_adtz'  ]['A24-9'], 4, 'CCSDTCF12-CP-Hill1_adtz')
    qcdb.compare_values(-4.5449, h2kc * mine['CCSDTCF12-CP-Hill1_atqz'  ]['A24-9'], 4, 'CCSDTCF12-CP-Hill1_atqz')
    qcdb.compare_values(-4.5400, h2kc * mine['CCSDTCF12-CP-Hill1_dtzf12']['A24-9'], 4, 'CCSDTCF12-CP-Hill1_dtzf12')
    qcdb.compare_values(-4.5539, h2kc * mine['CCSDTCF12-CP-Hill1_tqzf12']['A24-9'], 4, 'CCSDTCF12-CP-Hill1_tqzf12')
    qcdb.compare_values(-4.1343, h2kc * mine['CCSDAF12-CP-atz'          ]['A24-9'], 4, 'CCSDAF12-CP-atz')
    qcdb.compare_values(-4.1111, h2kc * mine['CCSDAF12-CP-a5z'          ]['A24-9'], 4, 'CCSDAF12-CP-a5z')
    qcdb.compare_values(-4.1183, h2kc * mine['CCSDAF12-CP-aqz'          ]['A24-9'], 4, 'CCSDAF12-CP-aqz')
    qcdb.compare_values(-4.5610, h2kc * mine['CCSDTNSAF12-CP-atz'       ]['A24-9'], 4, 'CCSDTNSAF12-CP-atz')
    qcdb.compare_values(-4.4933, h2kc * mine['CCSDTNSBF12-CP-atz'       ]['A24-9'], 4, 'CCSDTNSBF12-CP-atz')
    #qcdb.compare_values(, h2kc * mine['CCSDTNSCF12-CP-atz'       ]['A24-9'], 4, 'CCSDTNSCF12-CP-atz')
    qcdb.compare_values(-4.5892, h2kc * mine['CCSDTNSBF12-CP-Hill2_adtz']['A24-9'], 4, 'CCSDTNSBF12-CP-Hill2_adtz')

    # numbers in 6 checks below aren't right 
    #qcdb.compare_values(-4.5828, h2kc * mine['CCSDTNSAF12-CP-Hill2_dtzf12']['A24-9'], 4, 'CSDTNSAF12-CP-Hill2_dtzf12')
    #qcdb.compare_values(-4.5615, h2kc * mine['CCSDTNSBF12-CP-Hill2_dtzf12']['A24-9'], 4, 'CSDTNSBF12-CP-Hill2_dtzf12')
    #qcdb.compare_values(-4.5475, h2kc * mine['CCSDTNSCF12-CP-Hill2_dtzf12']['A24-9'], 4, 'CSDTNSCF12-CP-Hill2_dtzf12')
    #qcdb.compare_values(-4.5797, h2kc * mine['CCSDTNSAF12-CP-Hill2_tqzf12']['A24-9'], 4, 'CSDTNSAF12-CP-Hill2_tqzf12')
    #qcdb.compare_values(-4.5580, h2kc * mine['CCSDTNSBF12-CP-Hill2_tqzf12']['A24-9'], 4, 'CSDTNSBF12-CP-Hill2_tqzf12')
    #qcdb.compare_values(-4.5499, h2kc * mine['CCSDTNSCF12-CP-Hill2_tqzf12']['A24-9'], 4, 'CSDTNSCF12-CP-Hill2_tqzf12')
    
except KeyError as e:
    print(e)
    pass
else:
    print('End of Tests')


# <<< write qcdb data loader >>>

project = 'f12dilabio'
f1 = open('%s_%s.py' % (dbse, project), 'w')
f1.write('\ndef load_%s(dbinstance):\n\n' % (project))

for mc in mine.columns:
    method, bsse, basis = mc.split('-')  # TODO could be done better
    for rxn in rxns:
        database, reaction = rxn.split('-')
        value = h2kc * mine[mc][rxn]
        if pd.isnull(value):
            pass
        else:
            f1.write("""    dbinstance.add_ReactionDatum(dbse='%s', rxn=%s, method='%s', mode='%s', basis='%s', value=%.4f)\n""" %
                (database, reaction, method, bsse, basis, value))

f1.close()


# <<< collecting section >>>

#print 'indexval before', sorted(set([tup[1] for tup in df.index.values]))

# TODO write test case missing checker into the build loop

# if basis treatment isn't designed for SCF, i.e. adtz where Helgaker xtpl is meant for corl, should first item of build be None?

# form indices for mine from master df
#rxns = []
#[rxns.append(tup[2]) for tup in df.index.values if not rxns.count(tup[2])]

#print 'scf', df.loc['tzf12'].loc['HF-CABS TOTAL ENERGY']
#print 'cc', df.loc['HillCC_dtzf12'].loc['CCSD-F12B CORRELATION ENERGY']
#print '(t)', df.loc['HillT_dtzf12'].loc['(T)-F12AB CORRECTION ENERGY']
#print 'search', df.xs('(T)-F12AB CORRECTION ENERGY', level='psivar')
#print mine['CCSDTNSBF12-CP-Hill2_dtzf12']
#print mine.columns
