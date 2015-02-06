from collections import defaultdict

def muster_memory(mem):
    """Transform input *mem* in MB into psi4-type options for orca.

    """
    text = """%MaxCore {}\n""".format(int(mem))

    # prepare memory keywords to be set as c-side keywords
    options = defaultdict(lambda: defaultdict(dict))
    #options['ORCA']['ORCA_MAXCORE']['value'] = int(mem)
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

