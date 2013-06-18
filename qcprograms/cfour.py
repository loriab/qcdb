import re

def cfour_harvest(outtext):
    """
    """
    psivar = {}

    # Many regexes
    NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"

    scf1 = re.compile(r'^\s+' + r'(E\(SCF\))' + r'\s+=\s+' + NUMBER + r'\s+a\.u\.\s*$',
        re.MULTILINE)

    scf2 = re.compile(r'^\s+' + r'(E\(SCF\)=)' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        re.MULTILINE)
    scf3 = re.compile(r'^\s+' + r'(SCF has converged.)' + r'\s*$' + 
                      r'(.*?)' + 
                      r'^\s+' + r'(\d+)' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        re.MULTILINE | re.DOTALL | re.IGNORECASE)

    mp2 = re.compile(r'\s+' + r'(E2\(TOT\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*',
        re.MULTILINE | re.DOTALL)

    # Process SCF
    if re.search(scf1, outtext):
        psivar['SCF TOTAL ENERGY'] = float(re.search(scf1, outtext).group(2))
    elif re.search(scf2, outtext):
        psivar['SCF TOTAL ENERGY'] = float(re.search(scf2, outtext).group(2))
    elif re.search(scf3, outtext):
        psivar['SCF TOTAL ENERGY'] = float(re.search(scf3, outtext).group(4))

    if psivar.has_key('SCF TOTAL ENERGY'):
        psivar['CURRENT REFERENCE ENERGY'] = psivar['SCF TOTAL ENERGY']
        psivar['CURRENT ENERGY'] = psivar['SCF TOTAL ENERGY']

    # Process MP2
    if re.search(mp2, outtext):
        psivar['MP2 TOTAL ENERGY'] = float(re.search(mp2, outtext).group(2))

    print psivar
    return psivar

