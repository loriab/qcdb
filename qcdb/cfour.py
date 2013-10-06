import re
import collections
from decimal import Decimal
from pdict import PreservingDict
from periodictable import *
from molecule import Molecule
from orient import OrientMols


def harvest_output(outtext):
    """
    Function to separate portions of a CFOUR output file *outtest*,
    divided by xjoda.
    """
    pass_psivar = []
    pass_coord = []
    pass_grad = []

    for outpass in re.split(r'--invoking executable xjoda', outtext, re.MULTILINE):
        psivar, c4coord, c4grad = harvest_outfile_pass(outpass)
        pass_psivar.append(psivar)
        pass_coord.append(c4coord)
        pass_grad.append(c4grad)

        #print '\n\nXXXXXXXXXXXXXXXXXXXXXXXXXX\n\n'
        #print outpass
        #print psivar, c4coord, c4grad
        #print psivar, c4grad
        #print '\n\nxxxxxxxxxxxxxxxxxxxxxxxxxx\n\n'

    retindx = -1 if pass_coord[-1] else -2

#    print '    <<<  C4 PSIVAR  >>>'
#    for item in pass_psivar[retindx]:
#        print('       %30s %16.8f' % (item, pass_psivar[retindx][item]))
#    print '    <<<  C4 COORD   >>>'
#    for item in pass_coord[retindx]:
#        print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
#    print '    <<<   C4 GRAD   >>>'
#    for item in pass_grad[retindx]:
#        print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

    return pass_psivar[retindx], pass_coord[retindx], pass_grad[retindx]


def harvest_outfile_pass(outtext):
    """
    Function to read CFOUR output file *outtext* and parse important quantum chemical information from it in
    """
    psivar = PreservingDict()
    psivar_coord = None
    psivar_grad = None

#    TODO: BCC
#          CI
#          QCISD(T)
#          other ROHF tests
#          vcc/ecc

    NUMBER = "((?:[-+]?\\d*\\.\\d+(?:[DdEe][-+]?\\d+)?)|(?:[-+]?\\d+\\.\\d*(?:[DdEe][-+]?\\d+)?))"

    # Process NRE
    mobj = re.search(r'^\s+' + r'(?:Nuclear repulsion energy :)' + r'\s+' + NUMBER + r'\s+a\.u\.\s*$',
        outtext, re.MULTILINE)

    if mobj:
        print('matched nre')
        psivar['NUCLEAR REPULSION ENERGY'] = mobj.group(1)
#        # Coordinates for atom-only calc
#        if float(mobj.group(1)) < 1.0E-6:
#            print('matched atom')
#            # Dinky Molecule and the element's invented TODO
#            molxyz = '1 bohr\n\n%s 0.0 0.0 0.0\n' % ('H')
#            psivar_coord = Molecule.init_with_xyz(molxyz, no_com=True, no_reorient=True, contentsNotFilename=True)

    # Process SCF
    mobj = re.search(
        r'^\s+' + r'(?:E\(SCF\))' + r'\s+=\s+' + NUMBER + r'\s+a\.u\.\s*$',
        outtext, re.MULTILINE)
    if mobj:
        print('matched scf1')
        psivar['SCF TOTAL ENERGY'] = mobj.group(1)

    mobj = re.search(
        r'^\s+' + r'(?:E\(SCF\)=)' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE)
    if mobj:
        print('matched scf2')
        psivar['SCF TOTAL ENERGY'] = mobj.group(1)

    mobj = re.search(
        r'^\s+' + r'(?:SCF has converged.)' + r'\s*$' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:\d+)' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched scf3')
        psivar['SCF TOTAL ENERGY'] = mobj.group(1)

    # Process MP2
    mobj = re.search(
        r'^\s+' + r'(?:E2\(AA\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(AB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(TOT\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:Total MP2 energy)' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*$',
        outtext, re.MULTILINE)
    if mobj:
        print('matched mp2r')
        psivar['MP2 SAME-SPIN ENERGY'] = 2 * Decimal(mobj.group(1))
        psivar['MP2 OPPOSITE-SPIN ENERGY'] = mobj.group(2)
        psivar['MP2 CORRELATION ENERGY'] = 2 * Decimal(mobj.group(1)) + Decimal(mobj.group(2))
        psivar['MP2 TOTAL ENERGY'] = mobj.group(4)

    mobj = re.search(
        r'^\s+' + r'(?:E2\(AA\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(BB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(AB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(TOT\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:Total MP2 energy)' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*$',
        outtext, re.MULTILINE)
    if mobj:
        print('matched mp2u')
        psivar['MP2 SAME-SPIN ENERGY'] = Decimal(mobj.group(1)) + Decimal(mobj.group(2))
        psivar['MP2 OPPOSITE-SPIN ENERGY'] = mobj.group(3)
        psivar['MP2 CORRELATION ENERGY'] = Decimal(mobj.group(1)) + \
            Decimal(mobj.group(2)) + Decimal(mobj.group(3))
        psivar['MP2 TOTAL ENERGY'] = mobj.group(5)

    mobj = re.search(
        r'^\s+' + r'(?:E2\(AA\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(BB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(AB\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(SINGLE\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:E2\(TOT\))' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'(?:Total MP2 energy)' + r'\s+=\s+' + NUMBER + r'\s+a.u.\s*$',
        outtext, re.MULTILINE)
    if mobj:
        print('matched mp2ro')
        psivar['MP2 SAME-SPIN ENERGY'] = Decimal(mobj.group(1)) + Decimal(mobj.group(2))
        psivar['MP2 OPPOSITE-SPIN ENERGY'] = mobj.group(3)
        psivar['MP2 SINGLES ENERGY'] = mobj.group(4)
        psivar['MP2 CORRELATION ENERGY'] = Decimal(mobj.group(1)) + \
            Decimal(mobj.group(2)) + Decimal(mobj.group(3)) + Decimal(mobj.group(4))
        psivar['MP2 TOTAL ENERGY'] = mobj.group(6)

    # Process MP3
    mobj = re.search(
        r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched mp3r')
        dmp2 = Decimal(mobj.group(1))
        dmp3 = Decimal(mobj.group(3))
        psivar['MP2 CORRELATION ENERGY'] = dmp2
        psivar['MP2 TOTAL ENERGY'] = mobj.group(2)
        psivar['MP3 CORRELATION ENERGY'] = dmp2 + dmp3
        psivar['MP3 TOTAL ENERGY'] = mobj.group(4)
        psivar['MP2.5 CORRELATION ENERGY'] = dmp2 + Decimal('0.500000000000') * dmp3
        psivar['MP2.5 TOTAL ENERGY'] = psivar['MP2.5 CORRELATION ENERGY'] + psivar['SCF TOTAL ENERGY']

    mobj = re.search(
        r'^\s+' + r'(?:S-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:S-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched mp3ro')
        dmp2 = Decimal(mobj.group(1)) + Decimal(mobj.group(3))
        dmp3 = Decimal(mobj.group(5)) + Decimal(mobj.group(7))
        psivar['MP3 CORRELATION ENERGY'] = dmp2 + dmp3
        psivar['MP3 TOTAL ENERGY'] = mobj.group(8)
        psivar['MP2.5 CORRELATION ENERGY'] = dmp2 + Decimal('0.500000000000') * dmp3
        psivar['MP2.5 TOTAL ENERGY'] = psivar['MP2.5 CORRELATION ENERGY'] + psivar['SCF TOTAL ENERGY']

    # Process MP4
    mobj = re.search(
        r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:Q-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:S-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched mp4r')
        dmp2 = Decimal(mobj.group(1))
        dmp3 = Decimal(mobj.group(3))
        dmp4sdq = Decimal(mobj.group(5)) + Decimal(mobj.group(7)) + Decimal(mobj.group(9))
        psivar['MP2 CORRELATION ENERGY'] = dmp2
        psivar['MP2 TOTAL ENERGY'] = mobj.group(2)
        psivar['MP3 CORRELATION ENERGY'] = dmp2 + dmp3
        psivar['MP3 TOTAL ENERGY'] = mobj.group(4)
        psivar['MP2.5 CORRELATION ENERGY'] = dmp2 + Decimal('0.500000000000') * dmp3
        psivar['MP2.5 TOTAL ENERGY'] = psivar['MP2.5 CORRELATION ENERGY'] + psivar['SCF TOTAL ENERGY']
        psivar['MP4(SDQ) CORRELATION ENERGY'] = dmp2 + dmp3 + dmp4sdq
        psivar['MP4(SDQ) TOTAL ENERGY'] = mobj.group(10)

    mobj = re.search(
        r'^\s+' + r'(?:S-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(2\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:S-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:D-MBPT\(3\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:L-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:NL-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched mp4ro')
        dmp2 = Decimal(mobj.group(1)) + Decimal(mobj.group(3))
        dmp3 = Decimal(mobj.group(5)) + Decimal(mobj.group(7))
        dmp4sdq = Decimal(mobj.group(9)) + Decimal(mobj.group(11))
        psivar['MP2 CORRELATION ENERGY'] = dmp2
        psivar['MP2 TOTAL ENERGY'] = mobj.group(4)
        psivar['MP3 CORRELATION ENERGY'] = dmp2 + dmp3
        psivar['MP3 TOTAL ENERGY'] = mobj.group(8)
        psivar['MP2.5 CORRELATION ENERGY'] = dmp2 + Decimal('0.500000000000') * dmp3
        psivar['MP2.5 TOTAL ENERGY'] = psivar['MP2.5 CORRELATION ENERGY'] + psivar['SCF TOTAL ENERGY']
        psivar['MP4(SDQ) CORRELATION ENERGY'] = dmp2 + dmp3 + dmp4sdq
        psivar['MP4(SDQ) TOTAL ENERGY'] = mobj.group(12)

    mobj = re.search(
        r'^\s+' + r'(?:D-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:Q-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:S-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:T-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched mp4tr')
        dmp4sdq = Decimal(mobj.group(1)) + Decimal(mobj.group(3)) + Decimal(mobj.group(5))
        dmp4t = Decimal(mobj.group(7))
        psivar['MP4(SDQ) CORRELATION ENERGY'] = psivar['MP3 CORRELATION ENERGY'] + dmp4sdq
        psivar['MP4(SDQ) TOTAL ENERGY'] = mobj.group(6)
        psivar['MP4(T) CORRECTION ENERGY'] = dmp4t
        psivar['MP4(SDTQ) CORRELATION ENERGY'] = psivar['MP3 CORRELATION ENERGY'] + dmp4sdq + dmp4t
        psivar['MP4(SDTQ) TOTAL ENERGY'] = mobj.group(8)
        psivar['MP4 CORRELATION ENERGY'] = psivar['MP4(SDTQ) CORRELATION ENERGY']
        psivar['MP4 TOTAL ENERGY'] = psivar['MP4(SDTQ) TOTAL ENERGY']

    mobj = re.search(
        r'^\s+' + r'(?:L-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:NL-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:WT12-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*' +
        r'^\s+' + r'(?:T-MBPT\(4\))' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched mp4tro')
        dmp4sdq = Decimal(mobj.group(1)) + Decimal(mobj.group(3))
        dmp4t = Decimal(mobj.group(5)) + Decimal(mobj.group(7))  # TODO: WT12 with T, not SDQ?
        psivar['MP4(SDQ) CORRELATION ENERGY'] = psivar['MP3 CORRELATION ENERGY'] + dmp4sdq
        psivar['MP4(SDQ) TOTAL ENERGY'] = mobj.group(4)
        psivar['MP4(T) CORRECTION ENERGY'] = dmp4t
        psivar['MP4(SDTQ) CORRELATION ENERGY'] = psivar['MP3 CORRELATION ENERGY'] + dmp4sdq + dmp4t
        psivar['MP4(SDTQ) TOTAL ENERGY'] = mobj.group(8)
        psivar['MP4 CORRELATION ENERGY'] = psivar['MP4(SDTQ) CORRELATION ENERGY']
        psivar['MP4 TOTAL ENERGY'] = psivar['MP4(SDTQ) TOTAL ENERGY']

    # Process CC Iterations
    mobj = re.search(
        r'^\s+' + r'(?P<fullCC>(?P<iterCC>CC(?:\w+))(?:\(T\))?)' + r'\s+(?:energy will be calculated.)\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:\d+)' + r'\s+' + NUMBER + r'\s+' + NUMBER + r'\s+DIIS\s*' +
        r'^\s*(?:-+)\s*' +
        r'^\s*(?:A miracle (?:has come|come) to pass. The CC iterations have converged.)\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched cc with full %s iterating %s' % (mobj.group('fullCC'), mobj.group('iterCC')))
        #print mobj.group(1), mobj.group(2), mobj.group(3), mobj.group(4)
        psivar['%s CORRELATION ENERGY' % (mobj.group('iterCC'))] = mobj.group(3)
        psivar['%s TOTAL ENERGY' % (mobj.group('iterCC'))] = mobj.group(4)

    # Process CC(T)
    mobj = re.search(
        r'^\s+' + r'(?:E\(SCF\))' + r'\s+=\s+' + NUMBER + r'\s+a\.u\.\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:E\(CCSD\))' + r'\s+=\s+' + NUMBER + r'\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:E\(CCSD\(T\)\))' + r'\s+=\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched ccsd(t) vcc')
        psivar['SCF TOTAL ENERGY'] = mobj.group(1)
        psivar['CCSD TOTAL ENERGY'] = mobj.group(2)
        psivar['(T) CORRECTION ENERGY'] = Decimal(mobj.group(3)) - Decimal(mobj.group(2))
        psivar['CCSD(T) CORRELATION ENERGY'] = Decimal(mobj.group(3)) - Decimal(mobj.group(1))
        psivar['CCSD(T) TOTAL ENERGY'] = mobj.group(3)

    mobj = re.search(
        r'^\s+' + r'(?:E\(SCF\))' + r'\s+=\s+' + NUMBER + r'\s+a\.u\.\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:CCSD energy)' + r'\s+' + NUMBER + r'\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'(?:Total perturbative triples energy:)' + r'\s+' + NUMBER + r'\s*' +
        r'^\s*(?:-+)\s*' +
        r'^\s+' + r'(?:CCSD\(T\) energy)' + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched ccsd(t) ecc')
        psivar['SCF TOTAL ENERGY'] = mobj.group(1)
        psivar['CCSD TOTAL ENERGY'] = mobj.group(2)
        psivar['(T) CORRECTION ENERGY'] = mobj.group(3)
        psivar['CCSD(T) CORRELATION ENERGY'] = Decimal(mobj.group(4)) - Decimal(mobj.group(2))
        psivar['CCSD(T) TOTAL ENERGY'] = mobj.group(4)

    mobj = re.search(
        r'^\s+' + r'(?:CCSD energy)' + r'\s+' + NUMBER + r'\s*' +
        r'^\s*(?:-+)\s*' +
        r'^\s+' + r'(?:CCSD\(T\) energy)' + r'\s+' + NUMBER + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:
        print('matched ccsd(t) lamb')
        #print mobj.group(1), mobj.group(2)
        psivar['CCSD TOTAL ENERGY'] = mobj.group(1)
        psivar['(T) CORRECTION ENERGY'] = Decimal(mobj.group(2)) - Decimal(mobj.group(1))
        psivar['CCSD(T) CORRELATION ENERGY'] = Decimal(mobj.group(2)) - psivar['SCF TOTAL ENERGY']
        psivar['CCSD(T) TOTAL ENERGY'] = mobj.group(2)

    # Process SCS-CC
    mobj = re.search(
        r'^\s+' + r'(?P<fullCC>(?P<iterCC>CC(?:\w+))(?:\(T\))?)' + r'\s+(?:energy will be calculated.)\s*' +
        r'(?:.*?)' +
        r'^\s*' + r'(?:@CCENRG-I, Correlation energies.)' + r'\s+(?:ECCAA)\s+' + NUMBER + r'\s*' +
        r'^\s+(?:ECCBB)\s+' + NUMBER + '\s*' +
        r'^\s+(?:ECCAB)\s+' + NUMBER + '\s*' +
        r'^\s+(?:Total)\s+' + NUMBER + '\s*',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:  # PRINT=2 to get SCS-CC components
        print('matched scscc')
        psivar['%s SAME-SPIN ENERGY' % (mobj.group('iterCC'))] = Decimal(mobj.group(3)) + Decimal(mobj.group(4))
        psivar['%s OPPOSITE-SPIN ENERGY' % (mobj.group('iterCC'))] = mobj.group(5)
        psivar['%s CORRELATION ENERGY' % (mobj.group('iterCC'))] = mobj.group(6)

    mobj = re.search(
        r'^\s+' + r'(?P<fullCC>(?P<iterCC>CC(?:\w+))(?:\(T\))?)' + r'\s+(?:energy will be calculated.)\s*' +
        r'(?:.*?)' +
        r'^\s+' + r'Amplitude equations converged in' + r'\s*\d+\s*' + r'iterations.\s*' +
        r'^\s+' + r'The AA contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'The BB contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'The AB contribution to the correlation energy is:\s+' + NUMBER + r'\s+a.u.\s*' +
        r'^\s+' + r'The total correlation energy is\s+' + NUMBER + r'\s+a.u.\s*' +
        r'(?:.*?)' +
        #r'^\s+' + r'The CC iterations have converged.' + r'\s*$',
        r'^\s+' + r'(?:A miracle come to pass. )?' + r'The CC iterations have converged.' + r'\s*$',
        outtext, re.MULTILINE | re.DOTALL)
    if mobj:  # PRINT=2 to get SCS components
        print('matched scscc2')
        psivar['%s SAME-SPIN ENERGY' % (mobj.group('iterCC'))] = Decimal(mobj.group(3)) + Decimal(mobj.group(4))
        psivar['%s OPPOSITE-SPIN ENERGY' % (mobj.group('iterCC'))] = mobj.group(5)
        psivar['%s CORRELATION ENERGY' % (mobj.group('iterCC'))] = mobj.group(6)

    # Process gradient
    mobj = re.search(
        r'\s+' + r'Molecular gradient' + r'\s*' +
        r'\s+' + r'------------------' + r'\s*' +
        r'\s+' + r'\n' +
        r'(?:(?:\s+[A-Z]+\s*#\d+\s+[xyz]\s+[-+]?\d+\.\d+\s*\n)+)' +  # optional, it seems
        r'\n\n' +  # optional, it seems
        r'((?:\s+[A-Z]+\s*#\d+\s+\d?\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)' +
        r'\n\n' +
        r'\s+' + 'Molecular gradient norm',
        outtext, re.MULTILINE)
    if mobj:
        print('matched molgrad')
        atoms = []
        psivar_grad = []
        for line in mobj.group(1).splitlines():
            lline = line.split()
            atoms.append(lline[0])
            #psivar_gradient.append([Decimal(lline[-3]), Decimal(lline[-2]), Decimal(lline[-1])])
            psivar_grad.append([float(lline[-3]), float(lline[-2]), float(lline[-1])])

    # Process geometry
    mobj = re.search(
#        r'\s+(?:-+)\s*' +
#        r'^\s+' + r'Z-matrix   Atomic            Coordinates (in bohr)' + r'\s*' +
        r'^\s+' + r'Symbol    Number           X              Y              Z' + r'\s*' +
        r'^\s+(?:-+)\s*' +
        r'((?:\s+[A-Z]+\s+[0-9]+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s+[-+]?\d+\.\d+\s*\n)+)' +
        r'^\s+(?:-+)\s*',
        outtext, re.MULTILINE)
    if mobj:
        print('matched geom')
        molxyz = '%d bohr\n\n' % len(mobj.group(1).splitlines())
        for line in mobj.group(1).splitlines():
            lline = line.split()
            molxyz += '%s %16s %16s %16s\n' % (lline[0], lline[-3], lline[-2], lline[-1])
        # Rather a dinky Molecule as no ghost, charge, or multiplicity
        psivar_coord = Molecule.init_with_xyz(molxyz, no_com=True, no_reorient=True, contentsNotFilename=True)
        #print molxyz

    # Process atom geometry
    mobj = re.search(
        r'^\s+' + r'@GETXYZ-I,     1 atoms read from ZMAT.' + r'\s*' +
        r'^\s+' + r'[0-9]+\s+([A-Z]+)\s+[0-9]+\s+' + NUMBER + r'\s*',
        outtext, re.MULTILINE)
    if mobj:
        print('matched atom')
        # Dinky Molecule
        molxyz = '1 bohr\n\n%s 0.0 0.0 0.0\n' % (mobj.group(1))
        psivar_coord = Molecule.init_with_xyz(molxyz, no_com=True, no_reorient=True, contentsNotFilename=True)

    # Process error codes
    mobj = re.search(
        r'^\s*' + r'--executable ' + r'(\w+)' + r' finished with status' + r'\s+' + r'([1-9][0-9]*)',
        outtext, re.MULTILINE)
    if mobj:
        print('matched error')
        print mobj.group(1), mobj.group(2)
        psivar['CFOUR ERROR CODE'] = mobj.group(2)

    # Process CURRENT energies (TODO: needs better way)
    if 'SCF TOTAL ENERGY' in psivar:
        psivar['CURRENT REFERENCE ENERGY'] = psivar['SCF TOTAL ENERGY']
        psivar['CURRENT ENERGY'] = psivar['SCF TOTAL ENERGY']

    if 'MP2 TOTAL ENERGY' in psivar and 'MP2 CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['MP2 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['MP2 TOTAL ENERGY']

    if 'MP3 TOTAL ENERGY' in psivar and 'MP3 CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['MP3 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['MP3 TOTAL ENERGY']

    if 'MP4 TOTAL ENERGY' in psivar and 'MP4 CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['MP4 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['MP4 TOTAL ENERGY']

#    if ('%s TOTAL ENERGY' % (mobj.group('fullCC')) in psivar) and \
#       ('%s CORRELATION ENERGY' % (mobj.group('fullCC')) in psivar):
#        psivar['CURRENT CORRELATION ENERGY'] = psivar['%s CORRELATION ENERGY' % (mobj.group('fullCC')]
#        psivar['CURRENT ENERGY'] = psivar['%s TOTAL ENERGY' % (mobj.group('fullCC')]

    if 'CCSD TOTAL ENERGY' in psivar and 'CCSD CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSD CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSD TOTAL ENERGY']

    if 'CCSD(T) TOTAL ENERGY' in psivar and 'CCSD(T) CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSD(T) CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSD(T) TOTAL ENERGY']

    if 'CC3 TOTAL ENERGY' in psivar and 'CC3 CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CC3 CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CC3 TOTAL ENERGY']

    if 'CCSDT TOTAL ENERGY' in psivar and 'CCSDT CORRELATION ENERGY' in psivar:
        psivar['CURRENT CORRELATION ENERGY'] = psivar['CCSDT CORRELATION ENERGY']
        psivar['CURRENT ENERGY'] = psivar['CCSDT TOTAL ENERGY']

    return psivar, psivar_coord, psivar_grad


def harvest(p4Mol, c4out, **largs):
    """Parses all the pieces of output from Cfour: the stdout in
    *c4out* and the contents of various scratch files like GRD stored
    in their namesake keys in *largs*. Since all Cfour output uses
    its own orientation and atom ordering for the given molecule,
    a qcdb.Molecule *p4Mol*, if supplied, is used to transform the
    Cfour output back into consistency with *p4Mol*.

    """
    # Collect results from output file and subsidiary files
    outPsivar, outMol, outGrad = harvest_output(c4out)

    if 'GRD' in largs:
        grdMol, grdGrad = harvest_GRD(largs['GRD'])
    else:
        grdMol, grdGrad = None, None

    if 'FCMFINAL' in largs:
        fcmHess = harvest_FCM(largs['FCMFINAL'])
    else:
        fcmHess = None

    # Reconcile the coordinate information: several cases
    #   Case                            p4Mol   GRD      Check consistency           Apply orientation?
    #   sp with mol thru cfour {}       None    None              outMol             N.C.
    #   opt with mol thru cfour {}      None    grdMol            outMol && grdMol   N.C.
    #   sp with mol thru molecule {}    p4Mol   None     p4Mol && outMol             p4Mol <-- outMol
    #   opt with mol thru molecule {}   p4Mol   grdMol   p4Mol && outMol && grdMol   p4Mol <-- grdMol

    if outMol:
#        try:
        if grdMol:
#            print 'grdMol', grdMol
            if abs(outMol.nuclear_repulsion_energy() - grdMol.nuclear_repulsion_energy()) > 1.0e-3:
                raise ValidationError("""Cfour outfile (NRE: %f) inconsistent with Cfour GRD (NRE: %f).""" % \
                        (outMol.nuclear_repulsion_energy(), grdMol.nuclear_repulsion_energy()))
#        except:
#            print 'ack'
        if p4Mol:
            if abs(outMol.nuclear_repulsion_energy() - p4Mol.nuclear_repulsion_energy()) > 1.0e-3:
                raise ValidationError("""Cfour outfile (NRE: %f) inconsistent with Psi4 input (NRE: %f).""" % \
                    (outMol.nuclear_repulsion_energy(), p4Mol.nuclear_repulsion_energy()))
    else:
        raise ValidationError("""No coordinate information extracted from Cfour output.""")

#    try:
    print '    <<<   [1] P4-MOL   >>>'
    if p4Mol:
        p4Mol.print_out_in_bohr()
    print '    <<<   [2] C4-OUT-MOL   >>>'
    if outMol:
        outMol.print_out_in_bohr()
    print '    <<<   [3] C4-GRD-MOL   >>>'
    if grdMol:
        grdMol.print_out_in_bohr()
#    except UnboundLocalError:
#        pass

    # Set up array reorientation object
    if p4Mol and grdMol:
        p4c4 = OrientMols(p4Mol, grdMol)
        print p4c4
        oriCoord = p4c4.transform_coordinates2(grdMol)
        oriGrad = p4c4.transform_gradient(grdGrad)
    elif p4Mol and outMol:
        p4c4 = OrientMols(p4Mol, outMol)
        print p4c4
        oriCoord = p4c4.transform_coordinates2(outMol)
        oriGrad = None
    elif outMol:
        oriCoord = None
        oriGrad = None

    print '    <<<   [4] C4-ORI-MOL   >>>'
    if oriCoord is not None:
        for item in oriCoord:
            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))

    print '    <<<   [1] C4-GRD-GRAD   >>>'
    if grdGrad is not None:
        for item in grdGrad:
            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))
    print '    <<<   [2] C4-ORI-GRAD   >>>'
    if oriGrad is not None:
        for item in oriGrad:
            print('       %16.8f %16.8f %16.8f' % (item[0], item[1], item[2]))




#    p4c4 = OrientMols(p4Mol, grdMol)
#    print p4c4
#    print p4c4.transform_coordinates2(grdMol)
#    print 'pre grad'
#    oriGrad = p4c4.transform_gradient(grdGrad)
#    print 'post grad'







    #p4mol                             c4mol
    #mol.atom_map_and_orient_from_cfour(newmol)



#    if not p4Mol and not grdMol:
#        pass
#    elif not p4Mol and grdMol:
#        if abs(outMol.nuclear_repulsion_energy() - grdMol.nuclear_repulsion_energy()) > 1.0e-3:
#            raise ValidationError("""Cfour outfile (NRE: %f) inconsistent with Cfour GRD (NRE: %f).""" % \
#                (outMol.nuclear_repulsion_energy(), grdMol.nuclear_repulsion_energy()))
#        pass
#    elif p4Mol and not grdMol:
#        if abs(outMol.nuclear_repulsion_energy() - p4Mol.nuclear_repulsion_energy()) > 1.0e-3:
#            raise ValidationError("""Cfour outfile (NRE: %f) inconsistent with Psi4 input (NRE: %f).""" % \
#                (outMol.nuclear_repulsion_energy(), p4Mol.nuclear_repulsion_energy()))
#        pass
#    elif p4Mol and grdMol:
#        if (abs(outMol.nuclear_repulsion_energy() - grdMol.nuclear_repulsion_energy()) > 1.0e-3) or \
#           (abs(outMol.nuclear_repulsion_energy() - p4Mol.nuclear_repulsion_energy()) > 1.0e-3):
#            raise ValidationError("""Cfour outfile (NRE: %f) inconsistent with Cfour GRD (NRE: %f) or Psi4 input (NRE: %f).""" % \
#                (outMol.nuclear_repulsion_energy(), grdMol.nuclear_repulsion_energy(), p4Mol.nuclear_repulsion_energy()))
#        pass
#    else:
#        raise ValidationError("""Inexplicable pattern of Psi4 molecule and Cfour GRD information.""")




    #p4mol                             c4mol
    #mol.atom_map_and_orient_from_cfour(newmol)

    if oriGrad:
        retGrad = oriGrad
    elif grdGrad:
        retGrad = grdGrad
    else:
        retGrad = None

    return outPsivar, retGrad


def harvest_GRD(grd):
    """Parses the contents *grd* of the Cfour GRD file into the gradient
    array and coordinate information. The coordinate info is converted
    into a rather dinky Molecule (no charge, multiplicity, or fragment),
    but this is these coordinates that govern the reading of molecule
    orientation by Cfour. Return qcdb.Molecule and gradient array.

    """
    grd = grd.splitlines()
    Nat = int(grd[0].split()[0])
    molxyz = '%d bohr\n\n' % (Nat)

    grad = []
    for at in range(Nat):
        mline = grd[at + 1].split()
        el = 'GH' if int(float(mline[0])) == 0 else z2el[int(float(mline[0]))]
        molxyz += '%s %16s %16s %16s\n' % (el, mline[-3], mline[-2], mline[-1])
        lline = grd[at + 1 + Nat].split()
        grad.append([float(lline[-3]), float(lline[-2]), float(lline[-1])])
    mol = Molecule.init_with_xyz(molxyz, no_com=True, no_reorient=True, contentsNotFilename=True)

    return mol, grad
    

def harvest_FCM(fcm):
    """Parses the contents *fcm* of the Cfour FCMFINAL file into a hessian array.

    """
    fcm = fcm.splitlines()
    Nat = int(fcm[0].split()[0])
    Ndof = int(fcm[0].split()[1])

    empty = True
    hess = []
    for df in range(Ndof):
        for at in range(Nat):
            lline = fcm[Ndof * at + at + 1].split()
            if empty:
                if (abs(float(lline[0])) > 1.0e-8) or \
                   (abs(float(lline[1])) > 1.0e-8) or \
                   (abs(float(lline[2])) > 1.0e-8):
                    empty = False
            fcm.append([float(lline[0]), float(lline[1]), float(lline[2])])

    return None if empty else hess

#def cfour_harvest_files(mol, grd=None, fcmfinal=None):
#    """
#    grad = []
#
#    if grd:
#        grd = grd.splitlines()
#        Nat = int(grd[0].split()[0])
#        molxyz = '%d bohr\n\n' % (Nat)
#        for at in range(Nat):
#            molxyz += grd[at + 1] + '\n'
#            lline = grd[at + 1 + Nat].split()
#            grad.append([float(lline[-3]), float(lline[-2]), float(lline[-1])])
#
#        newmol = Molecule.init_with_xyz(molxyz, no_com=True, no_reorient=True, contentsNotFilename=True)
#        newmol.print_out()
#
#        mol.atom_map_and_orient_from_cfour(newmol)
#
#def cfour_harvest_GRD(outtext):
#    """Function to
#
#    """
#    atoms = []
#    coord = []
#    grad = []
#
#    outtext = outtext.splitlines()
#    Nat = int(outtext[0].split()[0])
#    xyzstring = '%d bohr\n\n' % (Nat)
#    for at in range(Nat):
#        lline = outtext[at + 1].split()
#        atoms.append(float(lline[0]))
#        coord.append([float(lline[-3]), float(lline[-2]), float(lline[-1])])
#        xyzstring += outtext[at + 1] + '\n'
#    for at in range(Nat):
#        lline = outtext[at + 1 + Nat].split()
#        if abs(float(lline[0]) - atoms[at]) > 1.0E-6:
#            raise qcdb.exceptions.ValidationError("""Inconsistent CFOUR GRD file.""")
#        grad.append([float(lline[-3]), float(lline[-2]), float(lline[-1])])
#
#    newmol = Molecule.init_with_xyz(xyzstring, no_com=True, no_reorient=True, contentsNotFilename=True)
#    newmol.print_out()


def cfour_memory(mem):
    """Transform input *mem* in MB into psi4-type options for cfour.

    """
    text = ''

    # prepare memory keywords to be set as c-side keywords
    options = collections.defaultdict(lambda: collections.defaultdict(dict))
    options['CFOUR']['CFOUR_MEMORY_SIZE']['value'] = int(mem)
    options['CFOUR']['CFOUR_MEM_UNIT']['value'] = 'MB'

    return text, options


def cfour_calclevel(dertype):
    """
    """
    text = ''
    options = collections.defaultdict(lambda: collections.defaultdict(dict))

    if dertype == 0:
        pass
    elif dertype == 1:
        options['CFOUR']['CFOUR_DERIV_LEVEL']['value'] = 'FIRST'
    else:
        print('bad dertype')
        exit()

    return text, options


def cfour_method(name):
    """Function to

    """
    lowername = name.lower()
    text = ''

    options = collections.defaultdict(lambda: collections.defaultdict(dict))

    if lowername == 'c4-scf':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'SCF'

    elif lowername == 'c4-mp2':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'MP2'

    elif lowername == 'c4-mp3':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'MP3'

    elif lowername == 'c4-mp4(sdq)':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'SDQ-MP4'

    elif lowername == 'c4-mp4':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'MP4'

    elif lowername == 'c4-ccsd':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CCSD'
        options['CFOUR']['CFOUR_CC_PROGRAM']['value'] = 'ECC'

    elif lowername == 'c4-cc3':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CC3'

    elif lowername == 'c4-ccsd(t)':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CCSD(T)'
        options['CFOUR']['CFOUR_CC_PROGRAM']['value'] = 'ECC'

    elif lowername == 'c4-ccsdt':
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CCSDT'

    return text, options


def cfour_list():
    """Return an array of Cfour methods with energies. Appended
    to procedures['energy'].

    """
    val = []
    val.append('cfour')
    val.append('c4-scf')
    val.append('c4-mp2')
    val.append('c4-mp3')
    val.append('c4-mp4(sdq)')
    val.append('c4-mp4')
    val.append('c4-ccsd')
    val.append('c4-cc3')
    val.append('c4-ccsd(t)')
    val.append('c4-ccsdt')
    return val


def cfour_gradient_list():
    """Return an array of Cfour methods with analytical gradients.
    Appended to procedures['gradient'].

    """
    val = []
    val.append('cfour')
    val.append('c4-scf')
    val.append('c4-mp2')
    val.append('c4-mp3')
    val.append('c4-mp4(sdq)')
    val.append('c4-mp4')
    val.append('c4-ccsd')
    val.append('c4-cc3')
    val.append('c4-ccsd(t)')
    val.append('c4-ccsdt')
    return val


def cfour_psivar_list():
    """Return a dict with keys of most Cfour methods and values of dicts
    with the PSI Variables returned by those methods. Used by cbs()
    wrapper to avoid unnecessary computations in compound methods.
    Result is appended to ``VARH``.

    """
    VARH = {}
    VARH['c4-scf'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY'}
    VARH['c4-mp2'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY',
                        'c4-mp2corl': 'MP2 CORRELATION ENERGY'}
    VARH['c4-mp3'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY',
                        'c4-mp2corl': 'MP2 CORRELATION ENERGY',
                      'c4-mp2.5corl': 'MP2.5 CORRELATION ENERGY',
                        'c4-mp3corl': 'MP3 CORRELATION ENERGY'}
    VARH['c4-mp4(sdq)'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY',
                        'c4-mp2corl': 'MP2 CORRELATION ENERGY',
                      'c4-mp2.5corl': 'MP2.5 CORRELATION ENERGY',
                        'c4-mp3corl': 'MP3 CORRELATION ENERGY',
                   'c4-mp4(sdq)corl': 'MP4(SDQ) CORRELATION ENERGY'}
    VARH['c4-mp4'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY',
                        'c4-mp2corl': 'MP2 CORRELATION ENERGY',
                      'c4-mp2.5corl': 'MP2.5 CORRELATION ENERGY',
                        'c4-mp3corl': 'MP3 CORRELATION ENERGY',
                   'c4-mp4(sdq)corl': 'MP4(SDQ) CORRELATION ENERGY',
                        'c4-mp4corl': 'MP4(SDTQ) CORRELATION ENERGY'}
    VARH['c4-ccsd'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY',
                        'c4-mp2corl': 'MP2 CORRELATION ENERGY',
                       'c4-ccsdcorl': 'CCSD CORRELATION ENERGY'}
    VARH['c4-cc3'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY',
                        'c4-mp2corl': 'MP2 CORRELATION ENERGY',
                        'c4-cc3corl': 'CC3 CORRELATION ENERGY'}
    VARH['c4-ccsd(t)'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY',
                        'c4-mp2corl': 'MP2 CORRELATION ENERGY',
                       'c4-ccsdcorl': 'CCSD CORRELATION ENERGY',
                    'c4-ccsd(t)corl': 'CCSD(T) CORRELATION ENERGY'}
    VARH['c4-ccsdt'] = {
                         'c4-scftot': 'SCF TOTAL ENERGY',
                        'c4-mp2corl': 'MP2 CORRELATION ENERGY',
                       'c4-ccsdcorl': 'CCSD CORRELATION ENERGY',
                      'c4-ccsdtcorl': 'CCSDT CORRELATION ENERGY'}

    return VARH
