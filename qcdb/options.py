import math
from exceptions import *


def format_option_for_cfour(opt, val):
    """Function to reformat value *val* for option *opt* from python
    into cfour-speak. Arrays are the primary target.

    """
    text = ''

    # Transform list from [[3, 0, 1, 1], [2, 0, 1, 0]] --> 3-0-1-1/2-0-1-0
    if isinstance(val, list):
        if type(val[0]).__name__ == 'list':
            if type(val[0][0]).__name__ == 'list':
                raise ValidationError('Option has level of array nesting inconsistent with CFOUR.')
            else:
                # option is 2D array
                for no in range(len(val)):
                    for ni in range(len(val[no])):
                        text += str(val[no][ni])
                        if ni < (len(val[no]) - 1):
                            text += '-'
                    if no < (len(val) - 1):
                        text += '/'
        else:
            # option is plain 1D array
            for n in range(len(val)):
                text += str(val[n])
                if n < (len(val) - 1):
                    text += '-'

    # Transform the basis sets that *must* be lowercase (dratted c4 input)
    elif (opt == 'CFOUR_BASIS') and (val.upper() in ['SVP', 'DZP', 'TZP', 'TZP2P', 'QZ2P', 'PZ3D2F', '13S9P4D3F']):
        text += str(val.lower())

    # Transform convergence options to order-of-magnitude
#    elif (opt.upper() in ['CFOUR_BRUCK_CONV', 'CFOUR_CC_CONV', 'CFOUR_CIS_CONV', 
#        'CFOUR_CONVERGENCE', 'CFOUR_CPHF_CONVER', 'CFOUR_ESTATE_CONV', 
#        'CFOUR_GEO_CONV', 'CFOUR_LINEQ_CONV', 'CFOUR_SCF_CONV']):
#        print val, math.log10(val)
#        text += str(-1 * math.log10(val))
#        print 'asdf', val, text

    # No Transform
    else:
        text += str(val)

    return opt[6:], text


def prepare_options_for_cfour(options):
    """Function to take the full snapshot of the liboptions object
    encoded in dictionary *options*, find the options directable toward
    Cfour (options['CFOUR']['CFOUR_*']) that aren't default, then write
    a *CFOUR deck with those options.

    """
    text = ''

    for opt, val in options['CFOUR'].items():
        if opt.startswith('CFOUR_'):
            if val['has_changed']:
                if not text:
                    text += """*CFOUR("""
                text += """%s=%s\n""" % (format_option_for_cfour(opt, val['value']))
    if text:
        text = text[:-1] + ')\n\n'

    return text


def reconcile_options(full, partial):
    """Function to take the full snapshot of the liboptions object
    encoded in dictionary *full* and reconcile it with proposed options
    value changes in *partial*. Overwrites *full* with *partial* if
    option untouched, touches *full* if *full* and *partial* are in
    agreement, balks if *full* and *partial* conflict. Returns *full*.

    """
    for okey, oval in partial.items():
        for ikey, ival in oval.items():
            if full[okey][ikey]['has_changed']:
                if full[okey][ikey]['value'] != ival['value']:
                    raise ValidationError("""Option %s value `%s` set by options block incompatible with value `%s` in memory/molecule/psi4options block.""" %
                        (ikey, full[okey][ikey]['value'], ival['value']))
                else:
                    # kw in full is touched, but in agreement with value in partial, no change
                    full[okey][ikey]['has_changed'] = True
            else:
                # If kw in full is untouched, overwrite it with value in partial
                full[okey][ikey]['value'] = ival['value']
                full[okey][ikey]['has_changed'] = True
                print 'Overwriting %s with %s' % (ikey, ival['value'])

    return full
