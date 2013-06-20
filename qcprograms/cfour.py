import re

def cfour_harvest(outtext):
    """
    """
    psivar = {}

    
    # Many regexes
    NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"

    nre = re.compile(r'^\s+' + r'(?:Nuclear repulsion energy :)' + r'\s+' + NUMBER + r'\s+a\.u\.\s*$',
        re.MULTILINE)

    scf1 = re.compile(r'^\s+' + r'(E\(SCF\))' + r'\s+=\s+' + NUMBER + r'\s+a\.u\.\s*$',
        re.MULTILINE)

    scf2 = re.compile(r'^\s+' + r'(E\(SCF\)=)' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        re.MULTILINE)

    scf3 = re.compile(r'^\s+' + r'(SCF has converged.)' + r'\s*$' + 
                      r'(.*?)' + 
                      r'^\s+' + r'(\d+)' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        re.MULTILINE | re.DOTALL)

    mp2r = re.compile(r'^\s+' + r'(E2\(AA\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' + 
                      r'^\s+' + r'(E2\(AB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' + 
                      r'^\s+' + r'(E2\(TOT\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' + 
                      r'^\s+' + r'(Total MP2 energy)' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*$', 
        re.MULTILINE)

    mp2u = re.compile(r'^\s+' + r'(E2\(AA\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' + 
                      r'^\s+' + r'(E2\(BB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' + 
                      r'^\s+' + r'(E2\(AB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' + 
                      r'^\s+' + r'(E2\(TOT\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' + 
                      r'^\s+' + r'(Total MP2 energy)' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*$', 
        re.MULTILINE)

    mp2ro = re.compile(r'^\s+' + r'(E2\(AA\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' + 
                       r'^\s+' + r'(E2\(BB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' + 
                       r'^\s+' + r'(E2\(AB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' + 
                       r'^\s+' + r'(E2\(SINGLE\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' + 
                       r'^\s+' + r'(E2\(TOT\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' + 
                       r'^\s+' + r'(Total MP2 energy)' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*$', 
        re.MULTILINE)

    mp3r = re.compile(r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
                      r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        re.MULTILINE | re.DOTALL)
    mp3ro = re.compile(r'^\s+' + r'(?:S-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
                       r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
                       r'^\s+' + r'(?:S-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
                       r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        re.MULTILINE | re.DOTALL)
    mp4sdtqr = re.compile(r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
                          r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
                          r'^\s+' + r'(?:D-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
                          r'^\s+' + r'(?:Q-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
                          r'^\s+' + r'(?:S-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
                          r'^\s+' + r'(?:T-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        re.MULTILINE | re.DOTALL)
    mp4sdtqro = re.compile(r'^\s+' + r'(?:S-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
                           r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
                           r'^\s+' + r'(?:S-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
                           r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
                           r'^\s+' + r'(?:L-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
                           r'^\s+' + r'(?:NL-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
                           r'^\s+' + r'(?:WT12-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
                           r'^\s+' + r'(?:T-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        re.MULTILINE | re.DOTALL)

    ccsd = re.compile(r'^\s+' + r'(CCSD|CCSD\(T\))' + r'\s+(energy will be calculated.)\s*' +
                      r'(.*?)' +
                      r'^\s+' + r'(\d+)' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+DIIS\s*' +
                      r'^\s*(-+)\s*' +
                      r'^\s*(A miracle (?:has come|come) to pass. The CC iterations have converged.)\s*$',
        re.MULTILINE | re.DOTALL)

    scsccsd = re.compile(r'^\s*' + r'(@CCENRG-I, Correlation energies.)' + r'\s+(ECCAA)\s+' + NUMBER + r'\s*' +
                         r'^\s+(ECCBB)\s+' + NUMBER + '\s*' +
                         r'^\s+(ECCAB)\s+' + NUMBER + '\s*' +
                         r'^\s+(Total)\s+' + NUMBER + '\s*',
        re.MULTILINE | re.DOTALL)

    scscc = re.compile(r'^\s+' + r'Amplitude equations converged in' + r'\s*\d+\s*' + r'iterations.\s*' +
                       r'^\s+' + r'The AA contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
                       r'^\s+' + r'The BB contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
                       r'^\s+' + r'The AB contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
                       r'^\s+' + r'The total correlation energy is\s+' + NUMBER + r'\s+a.u.\s*' +
                       r'^\s+' + r'The CC iterations have converged.' + r'\s*$',
        re.MULTILINE | re.DOTALL)

    ccsd_t_vcc = re.compile(r'^\s+' + r'(E\(SCF\))' + r'\s+=\s+' + NUMBER + r'\s+a\.u\.\s*' +
                            r'(.*?)' +
                            r'^\s+' + r'(E\(CCSD\))' + r'\s+=\s+' + NUMBER + r'\s*' +
                            r'(.*?)' +
                            r'^\s+' + r'(E\(CCSD\(T\)\))' + r'\s+=\s+' + NUMBER + r'\s*$',
        re.MULTILINE | re.DOTALL)

    ccsd_t_ecc = re.compile(r'^\s+' + r'(E\(SCF\))' + r'\s+=\s+' + NUMBER + r'\s+a\.u\.\s*' +
                            r'(.*?)' +
                            r'^\s+' + r'(CCSD energy)' + r'\s+' + NUMBER + r'\s*' +
                            r'(.*?)' +
                            r'^\s+' + r'(Total perturbative triples energy:)' + r'\s+' + NUMBER + r'\s*' +
                            r'^\s*(-+)\s*' +
                            r'^\s+' + r'(CCSD\(T\) energy)' + r'\s+' + NUMBER + r'\s*$',
        re.MULTILINE | re.DOTALL)

    cc3 = re.compile(r'^\s+' + r'(CC3)' + r'\s+(energy will be calculated.)\s*' +
                      r'(.*?)' +
                      r'^\s+' + r'(\d+)' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+DIIS\s*' +
                      r'^\s*(-+)\s*' +
                      r'^\s*(A miracle (?:has come|come) to pass. The CC iterations have converged.)\s*$',
        re.MULTILINE | re.DOTALL)

    ccsdt = re.compile(r'^\s+' + r'(CCSDT)' + r'\s+(energy will be calculated.)\s*' +
                       r'(.*?)' +
                       r'^\s+' + r'(\d+)' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+DIIS\s*' +
                       r'^\s*(-+)\s*' +
                       r'^\s*(A miracle (?:has come|come) to pass. The CC iterations have converged.)\s*$',
        re.MULTILINE | re.DOTALL)

#              E(SCF)                =  -76.062748460180 a.u.
#  CCSD energy                              -76.338453952539
#  E4T  to CCSD(T)                           -0.007532615477
#  E5ST to CCSD(T)                            0.000269017447
#  Total perturbative triples energy:        -0.007263598030
#--------------------------------------------------------------------------------
#  CCSD(T) energy                     -76.345717550569

# CCSD(T) energy will be calculated.
# ...
#     20     -0.2757054923585205   -76.338453952539   DIIS
# -----------------------------------------------------------
# A miracle come to pass. The CC iterations have converged.

#   CCSDT   energy will be calculated.
# ...
#           36        -0.218848569882     -55.808195538491  DIIS
#      -----------------------------------------------------------
#      A miracle has come to pass. The CC iterations have converged.

# @CCENRG-I, Correlation energies. ECCAA       -0.027981812388
#                                  ECCBB       -0.011925433159
#                                  ECCAB       -0.173390810142
#                                  Total       -0.213298055690


    # Process NRE
    if re.search(nre, outtext):
        print('matched nre')
        psivar['NUCLEAR REPULSION ENERGY'] = float(re.search(nre, outtext).group(1))

    # Process SCF
    if re.search(scf1, outtext):
        print('matched scf1')
        psivar['SCF TOTAL ENERGY'] = float(re.search(scf1, outtext).group(2))
    elif re.search(scf2, outtext):
        print('matched scf2')
        psivar['SCF TOTAL ENERGY'] = float(re.search(scf2, outtext).group(2))
    elif re.search(scf3, outtext):
        print('matched scf3')
        psivar['SCF TOTAL ENERGY'] = float(re.search(scf3, outtext).group(4))
    if psivar.has_key('SCF TOTAL ENERGY'):
        psivar['CURRENT REFERENCE ENERGY'] = psivar['SCF TOTAL ENERGY']
        psivar['CURRENT ENERGY'] = psivar['SCF TOTAL ENERGY']

    # Process MP2
    if re.search(mp2r, outtext):
        print('matched mp2r')
        mobj = re.search(mp2r, outtext)
        psivar['MP2 SAME-SPIN ENERGY'] = 2 * float(mobj.group(2))
        psivar['MP2 OPPOSITE-SPIN ENERGY'] = float(mobj.group(4))
        psivar['MP2 CORRELATION ENERGY'] = 2 * float(mobj.group(2)) + float(mobj.group(4))
        psivar['MP2 TOTAL ENERGY'] = float(mobj.group(8))
    elif re.search(mp2u, outtext):  # elif not really necessary since matches are exclusive, just saving searches
        print('matched mp2u')
        mobj = re.search(mp2u, outtext)
        psivar['MP2 SAME-SPIN ENERGY'] = float(mobj.group(2)) + float(mobj.group(4))
        psivar['MP2 OPPOSITE-SPIN ENERGY'] = float(mobj.group(6))
        psivar['MP2 CORRELATION ENERGY'] = float(mobj.group(2)) + float(mobj.group(4)) + float(mobj.group(6))
        psivar['MP2 TOTAL ENERGY'] = float(mobj.group(10))
    elif re.search(mp2ro, outtext):
        print('matched mp2ro')
        mobj = re.search(mp2ro, outtext)
#        print(float(mobj.group(8)))
        psivar['MP2 SAME-SPIN ENERGY'] = float(mobj.group(2)) + float(mobj.group(4))
        psivar['MP2 OPPOSITE-SPIN ENERGY'] = float(mobj.group(6))
        psivar['MP2 SINGLES ENERGY'] = float(mobj.group(8))
        psivar['MP2 CORRELATION ENERGY'] = float(mobj.group(2)) + float(mobj.group(4)) + float(mobj.group(6)) + float(mobj.group(8))
        psivar['MP2 TOTAL ENERGY'] = float(mobj.group(12))
    if psivar.has_key('MP2 TOTAL ENERGY') and psivar.has_key('MP2 CORRELATION ENERGY'):
        psivar['CURRENT CORRELATION ENERGY'] = psivar['MP2 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['MP2 TOTAL ENERGY']

    # Process MP3
    if re.search(mp3r, outtext):
        print('matched mp3r')
        mobj = re.search(mp3r, outtext)
        lmp2 = float(mobj.group(1))
        lmp3 = float(mobj.group(3))
        #print(float(mobj.group(2)))
        #print(float(mobj.group(3)))
        #print(mobj.groups(12))
        psivar['MP3 CORRELATION ENERGY'] = lmp2 + lmp3
        psivar['MP3 TOTAL ENERGY'] = float(mobj.group(4))
        psivar['MP2.5 CORRELATION ENERGY'] = lmp2 + 0.5 * lmp3
        psivar['MP2.5 TOTAL ENERGY'] = psivar['MP2.5 CORRELATION ENERGY'] + psivar['SCF TOTAL ENERGY']
    if re.search(mp3ro, outtext):
        print('matched mp3ro')
        mobj = re.search(mp3ro, outtext)
        lmp2 = float(mobj.group(1)) + float(mobj.group(3))
        lmp3 = float(mobj.group(5)) + float(mobj.group(7))
        psivar['MP3 CORRELATION ENERGY'] = lmp2 + lmp3
        psivar['MP3 TOTAL ENERGY'] = float(mobj.group(8))
        psivar['MP2.5 CORRELATION ENERGY'] = lmp2 + 0.5 * lmp3
        psivar['MP2.5 TOTAL ENERGY'] = psivar['MP2.5 CORRELATION ENERGY'] + psivar['SCF TOTAL ENERGY']
    if psivar.has_key('MP3 TOTAL ENERGY') and psivar.has_key('MP3 CORRELATION ENERGY'):
        psivar['CURRENT CORRELATION ENERGY'] = psivar['MP3 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['MP3 TOTAL ENERGY']

    # Process MP4
    if re.search(mp4sdtqr, outtext):
        print('matched mp4sdtqr')
        mobj = re.search(mp4sdtqr, outtext)
        lmp2 = float(mobj.group(1))
        lmp3 = float(mobj.group(3))
        lmp4sdq = float(mobj.group(5)) + float(mobj.group(7)) + float(mobj.group(9))
        lmp4t = float(mobj.group(11))
        #print(float(mobj.group(2)))
        #print(float(mobj.group(3)))
        print(mobj.groups(12))
        psivar['MP3 CORRELATION ENERGY'] = lmp2 + lmp3
        psivar['MP3 TOTAL ENERGY'] = float(mobj.group(4))
        psivar['MP2.5 CORRELATION ENERGY'] = lmp2 + 0.5 * lmp3
        psivar['MP2.5 TOTAL ENERGY'] = psivar['MP2.5 CORRELATION ENERGY'] + psivar['SCF TOTAL ENERGY']
        #print(lmp2 + lmp3 + psivar['SCF TOTAL ENERGY'], mobj.group(8))
        #print(lmp2 + lmp3 + lmp4sdq + psivar['SCF TOTAL ENERGY'], mobj.group(14))
        #print(lmp2 + lmp3 + lmp4sdq + lmp4t + psivar['SCF TOTAL ENERGY'], mobj.group(16))
        psivar['MP4(SDQ) CORRELATION ENERGY'] = lmp2 + lmp3 + lmp4sdq
        psivar['MP4(SDQ) TOTAL ENERGY'] = float(mobj.group(10))
        psivar['MP4(T) CORRECTION ENERGY'] = lmp4t
        psivar['MP4(SDTQ) CORRELATION ENERGY'] = lmp2 + lmp3 + lmp4sdq + lmp4t
        psivar['MP4(SDTQ) TOTAL ENERGY'] = float(mobj.group(12))
        psivar['MP4 CORRELATION ENERGY'] = lmp2 + lmp3 + lmp4sdq + lmp4t
        psivar['MP4 TOTAL ENERGY'] = float(mobj.group(12))
    if re.search(mp4sdtqro, outtext):
        print('matched mp4sdtqro')
        mobj = re.search(mp4sdtqro, outtext)
        lmp2 = float(mobj.group(1)) + float(mobj.group(3))
        lmp3 = float(mobj.group(5)) + float(mobj.group(7))
        lmp4sdq = float(mobj.group(9)) + float(mobj.group(11)) + float(mobj.group(13))
        lmp4t = float(mobj.group(15))
        #print(float(mobj.group(2)))
        #print(float(mobj.group(3)))
        #print(float(mobj.group(4)))
        psivar['MP3 CORRELATION ENERGY'] = lmp2 + lmp3
        psivar['MP3 TOTAL ENERGY'] = float(mobj.group(8))
        psivar['MP2.5 CORRELATION ENERGY'] = lmp2 + 0.5 * lmp3
        psivar['MP2.5 TOTAL ENERGY'] = psivar['MP2.5 CORRELATION ENERGY'] + psivar['SCF TOTAL ENERGY']
        #print(lmp2 + lmp3 + psivar['SCF TOTAL ENERGY'], mobj.group(8))
        #print(lmp2 + lmp3 + lmp4sdq + psivar['SCF TOTAL ENERGY'], mobj.group(14))
        #print(lmp2 + lmp3 + lmp4sdq + lmp4t + psivar['SCF TOTAL ENERGY'], mobj.group(16))
        psivar['MP4(SDQ) CORRELATION ENERGY'] = lmp2 + lmp3 + lmp4sdq
        psivar['MP4(SDQ) TOTAL ENERGY'] = float(mobj.group(14))
        psivar['MP4(T) CORRECTION ENERGY'] = lmp4t
        psivar['MP4(SDTQ) CORRELATION ENERGY'] = lmp2 + lmp3 + lmp4sdq + lmp4t
        psivar['MP4(SDTQ) TOTAL ENERGY'] = float(mobj.group(16))
        psivar['MP4 CORRELATION ENERGY'] = lmp2 + lmp3 + lmp4sdq + lmp4t
        psivar['MP4 TOTAL ENERGY'] = float(mobj.group(16))
    if psivar.has_key('MP4 TOTAL ENERGY') and psivar.has_key('MP4 CORRELATION ENERGY'):
        psivar['CURRENT CORRELATION ENERGY'] = psivar['MP4 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['MP4 TOTAL ENERGY']

    # Process CCSD
    if re.search(ccsd, outtext):
        print('matched ccsd')
        mobj = re.search(ccsd, outtext)
        psivar['CCSD CORRELATION ENERGY'] = float(mobj.group(5))
        psivar['CCSD TOTAL ENERGY'] = float(mobj.group(6))
        if re.search(scsccsd, outtext):  # PRINT=2 to get SCS-CCSD components
            print('matched scsccsd')
            mobj = re.search(scsccsd, outtext)
            psivar['CCSD SAME-SPIN ENERGY'] = float(mobj.group(3)) + float(mobj.group(5))
            psivar['CCSD OPPOSITE-SPIN ENERGY'] = float(mobj.group(7))
        elif re.search(scscc, outtext):  # PRINT=2 to get SCS components
            print('matched scscc')
            mobj = re.search(scscc, outtext)
#            print(mobj.group(1))
#            print(mobj.group(2))
#            print(mobj.group(3))
            psivar['CCSD SAME-SPIN ENERGY'] = float(mobj.group(1)) + float(mobj.group(2))
            psivar['CCSD OPPOSITE-SPIN ENERGY'] = float(mobj.group(3))
    if psivar.has_key('CCSD TOTAL ENERGY') and psivar.has_key('CCSD CORRELATION ENERGY'):
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSD CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSD TOTAL ENERGY']

    # Process CCSD(T)
    if re.search(ccsd_t_vcc, outtext):
        print('matched ccsd(t) vcc')
        mobj = re.search(ccsd_t_vcc, outtext)
        psivar['(T) CORRECTION ENERGY'] = float(mobj.group(8)) - float(mobj.group(5))
        psivar['CCSD(T) CORRELATION ENERGY'] = float(mobj.group(8)) - float(mobj.group(2))
        psivar['CCSD(T) TOTAL ENERGY'] = float(mobj.group(8))
    if re.search(ccsd_t_ecc, outtext):
        print('matched ccsd(t) ecc')
        mobj = re.search(ccsd_t_ecc, outtext)
        psivar['(T) CORRECTION ENERGY'] = float(mobj.group(8))
        psivar['CCSD(T) CORRELATION ENERGY'] = float(mobj.group(11)) - float(mobj.group(2))
        psivar['CCSD(T) TOTAL ENERGY'] = float(mobj.group(11))
    if psivar.has_key('CCSD(T) TOTAL ENERGY') and psivar.has_key('CCSD(T) CORRELATION ENERGY'):
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSD(T) CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSD(T) TOTAL ENERGY']

    # Process CC3
    if re.search(cc3, outtext):
        print('matched cc3')
        mobj = re.search(cc3, outtext)
        psivar['CC3 CORRELATION ENERGY'] = float(mobj.group(5))
        psivar['CC3 TOTAL ENERGY'] = float(mobj.group(6))
        if re.search(scscc, outtext):  # PRINT=2 to get SCS components
            print('matched scscc')
            mobj = re.search(scscc, outtext)
            print(mobj.group(1))
            print(mobj.group(2))
            print(mobj.group(3))
            psivar['CC3 SAME-SPIN ENERGY'] = float(mobj.group(1)) + float(mobj.group(2))
            psivar['CC3 OPPOSITE-SPIN ENERGY'] = float(mobj.group(3))
    if psivar.has_key('CC3 TOTAL ENERGY') and psivar.has_key('CC3 CORRELATION ENERGY'):
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CC3 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CC3 TOTAL ENERGY']

    # Process CCSDT
    if re.search(ccsdt, outtext):
        print('matched ccsdt')
        mobj = re.search(ccsdt, outtext)
        psivar['CCSDT CORRELATION ENERGY'] = float(mobj.group(5))
        psivar['CCSDT TOTAL ENERGY'] = float(mobj.group(6))
    if psivar.has_key('CCSDT TOTAL ENERGY') and psivar.has_key('CCSDT CORRELATION ENERGY'):
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSDT CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSDT TOTAL ENERGY']

    return psivar
