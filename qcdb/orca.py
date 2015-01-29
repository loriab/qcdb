from collections import defaultdict
from pdict import PreservingDict
from molecule import Molecule


def harvest(p4Mol, orca_out, **largs):
    """Harvest variables, gradient, and the molecule from the output and other files"""

    # Split into lines as it is much easier to find what is needed
    out_lines = orca_out.split('\n')

    mol = harvest_molecule_from_outfile(out_lines)

    file_name = "NONE"
    grad = harvest_engrad(file_name)

    psivar = PreservingDict()
    dipole, magnitude = harvest_dipole(out_lines, psivar)

    harvest_scf_from_outfile(out_lines, psivar)

    # harvest energi(es)

    return psivar, grad, mol


def muster_memory(mem):
    """Transform input *mem* in MB into psi4-type options for orca.

    """
    text = ''

    # prepare memory keywords to be set as c-side keywords
    options = defaultdict(lambda: defaultdict(dict))
    options['ORCA']['ORCA_MAXCORE']['value'] = int(mem)
    #options['CFOUR']['CFOUR_MEM_UNIT']['value'] = 'MB' # orca expects always mb

    for item in options['ORCA']:
        options['ORCA'][item]['clobber'] = True
    return text, options


def muster_modelchem(name, dertype):
    """Transform calculation method *name* and derivative level *dertype*
    into options for orca. While deliberately requested pieces,
    generally |cfour__cfour_deriv_level| and |cfour__cfour_calc_level|,
    are set to complain if contradicted ('clobber' set to True), other
    'recommended' settings, like |cfour__cfour_cc_program|, can be
    countermanded by keywords in input file ('clobber' set to False).
    Occasionally, want these pieces to actually overcome keywords in
    input file ('superclobber' set to True).

    """
    text = ''
    lowername = name.lower()
    options = defaultdict(lambda: defaultdict(dict))

    if dertype == 0:
        options['ORCA']['ORCA_RUNTYP']['value'] = 'ENERGY'
    elif dertype == 1:
        options['ORCA']['ORCA_RUNTYP']['value'] = 'ENGRAD'
    #elif dertype == 2:
    #    options['ORCA']['ORCA_RUNTYP']['value'] = 'SECOND'
    else:
        raise ValidationError("""Requested Cfour dertype %d is not available.""" % (dertype))

    #if lowername == 'cfour':
    #    pass
    if lowername == 'oc-b3lyp':
        options['ORCA']['ORCA_FUNCTIONAL']['value'] = 'B3LYP_G'

    elif lowername == 'c4-ccsdt':
        pass
        options['CFOUR']['CFOUR_CALC_LEVEL']['value'] = 'CCSDT'
        options['CFOUR']['CFOUR_CC_PROGRAM']['value'] = 'ECC'

    else:
        raise ValidationError("""Requested Orca computational methods %d is not available.""" % (lowername))

    # Set clobbering
    if 'ORCA_RUNTYP' in options['ORCA']:
        options['ORCA']['ORCA_RUNTYP']['clobber'] = True
        options['ORCA']['ORCA_RUNTYP']['superclobber'] = True
    if 'ORCA_FUNCTIONAL' in options['ORCA']:
        options['ORCA']['ORCA_FUNCTIONAL']['clobber'] = True
        options['ORCA']['ORCA_FUNCTIONAL']['superclobber'] = True
    #if 'CFOUR_CC_PROGRAM' in options['CFOUR']:
    #    options['CFOUR']['CFOUR_CC_PROGRAM']['clobber'] = False

    return text, options


def orca_list():
    """Return an array of Orca methods with energies. Appended
    to procedures['energy'].

    """
    val = []
    val.append('oc-b3lyp')
    return val


def orca_gradient_list():
    """Return an array of Orca methods with analytical gradients.
    Appended to procedures['gradient'].

    """
    val = []
    val.append('oc-b3lyp')
    return val


def harvest_molecule_from_outfile(lines):
    """Return a molecule of the last geometry"""
    #Sample molecule block

    #----------------------------
    #CARTESIAN COORDINATES (A.U.)
    #----------------------------
    #  NO LB      ZA    FRAG    MASS        X           Y           Z
    #   0 O     8.0000    0    15.999         -0.043407801307192         -0.055556028344352          0.000000000000000
    #   1 H     1.0000    0     1.008          1.780497256508764         -0.017018089151928          0.000000000000000
    #   2 H     1.0000    0     1.008         -0.462170608038134          1.719154625261312          0.000000000000000
    #

    geom_start = find_start(lines, 'CARTESIAN COORDINATES (A.U.)')
    if geom_start == -1:
        return Molecule()

    # Geometry starts 3 lines after header and ends with a blank line
    geom_start += 3
    end = ''
    mol_str = ''
    for i, line in enumerate(lines[geom_start:], start=geom_start):
        if line == end:
            break
        num, atom, z, frag, mass, x, y, z = line.split()
        mol_str += '{} {} {} {}\n'.format(atom, x, y, z)

    return Molecule.init_with_xyz(mol_str)


def harvest_scf_from_outfile(lines, psivar):
    """Harvest SCF results from the SCF section of the output file"""
    #Sample SCF results block

    #----------------
    #TOTAL SCF ENERGY
    #----------------
    #
    #Total Energy       :          -76.02602169 Eh           -2068.77322 eV
    #
    #Components:
    #Nuclear Repulsion  :            9.12509697 Eh             248.30651 eV
    #Electronic Energy  :          -85.15111867 Eh           -2317.07974 eV
    #
    #One Electron Energy:         -123.01434123 Eh           -3347.39040 eV
    #Two Electron Energy:           37.86322256 Eh            1030.31067 eV
    #
    #Virial components:
    #Potential Energy   :         -151.99262033 Eh           -4135.92947 eV
    #Kinetic Energy     :           75.96659864 Eh            2067.15624 eV
    #Virial Ratio       :            2.00078223
    #
    #

    scf_start = find_start(lines, 'TOTAL SCF ENERGY')
    if scf_start == -1:
        return ''

    # Energies in SCF block
    psivar['SCF TOTAL ENERGY'] = float(lines[scf_start + 3].split()[3])
    psivar['NUCLEAR REPULSION ENERGY'] = float(lines[scf_start + 6].split()[3])


def harvest_dipole(lines):
    """Harvest the dipole, and return as a tuple (x, y, z)
    Multiple different dipole moments are output if post-HF calculations are
    run, resulting in highly similar blocks. It by default collects the last

    TODO: collect all the different types of dipole moments
    """
    #Sample dipole moment results block

    #-------------
    #DIPOLE MOMENT
    #-------------
    #                                X             Y             Z
    #Electronic contribution:     -0.11359      -0.14669      -0.00000
    #Nuclear contribution   :      0.61892       0.79867       0.00000
    #                        -----------------------------------------
    #Total Dipole Moment    :      0.50533       0.65198      -0.00000
    #                        -----------------------------------------
    #Magnitude (a.u.)       :      0.82489
    #Magnitude (Debye)      :      2.09670
    #
    #

    dipole_start = find_start(lines, 'DIPOLE MOMENT')

    # Dipole x, y, z are the last items 6 lines down in the dipole block
    dipole = tuple(map(float, lines[dipole_start + 6].split()[-3:]))
    # Dipole magnitude is 8 line down in the dipole block
    magnitude = float(lines[dipole_start + 8][-1])

    return dipole, magnitude


def harvest_mp2(lines, psivar):
    """Harvest the MP2 results"""
    # Sample MP2 energy line (works for both MP2 and RI-MP2)

    #---------------------------------------
    #MP2 TOTAL ENERGY:      -76.226803665 Eh
    #---------------------------------------

    for line in lines:
        if line[:16] == 'MP2 TOTAL ENERGY':
            #psivar['MP2 ENERGY'] = float(line.split()[-2])
            break


def harvest_coupled_cluster(lines, psivar):
    """Harvest the coupled cluster results
    WARNING: Canonical and DLPNO print out the coupled cluster results differently
    """
    # Sample (canonical) coupled cluster results block

    #----------------------
    #COUPLED CLUSTER ENERGY
    #----------------------
    #
    #E(0)                                       ...    -76.063720080
    #E(CORR)                                    ...     -0.288938791
    #E(TOT)                                     ...    -76.352658871
    #Singles Norm <S|S>**1/2                    ...      0.021106262
    #T1 diagnostic                              ...      0.007462191
    #

    #Sample DLPNO coupled cluster block (CCSD)

    #----------------------
    #COUPLED CLUSTER ENERGY
    #----------------------
    #
    #E(0)                                       ...    -76.026019996
    #E(CORR)(strong-pairs)                      ...     -0.211953159
    #E(CORR)(weak-pairs)                        ...     -0.000007244
    #E(CORR)(corrected)                         ...     -0.211960403
    #E(TOT)                                     ...    -76.237980399
    #Singles Norm <S|S>**1/2                    ...      0.014443573
    #T1 diagnostic                              ...      0.005106574
    #

    cc_start = find_start(lines, 'COUPLED CLUSTER ENERGY')
    if cc_start == -1:
        return

    cc_reference = float(lines[cc_start + 3].split()[-1])


def harvest_engrad(engrad):
    """Parse the engrad file for the gradient"""
    mol = Molecule('')
    grad = []
    return mol, grad


def find_start(lines, start_str, reverse=True):
    """Find the start of a block, iterate backwards by default,
    Usually the last one is wanted
    If not found, return -1
    """
    start = -1
    # Iterate backwards until the last value is found
    if reverse:
        for i, line in reversed(list(enumerate(lines))):
            if start_str == line:
                return i
    else:
        for i, line in enumerate(lines):
            if start_str == line:
                return i
    return start
