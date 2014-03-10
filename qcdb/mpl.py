"""Module with matplotlib plotting routines. These are not hooked up to
any particular qcdb data structures but can be called with basic
arguments.

"""


def thread(data, labels, color=None, title='', xlimit=4.0, mae=None, mape=None):
    """Generates a tiered slat diagram between model chemistries with
    errors in list *data*, which is supplied as part of the dictionary for
    each participating reaction, along with *dbse* and *rxn* keys in
    argument *data*. The plot is labeled with *title* and each tier with
    an element of *labels* and plotted at *xlimit* from the zero-line. If
    *color* is None, slats are black, if 'sapt', colors are taken from
    sapt_colors module. Summary statistics *mae* are plotted on the
    overbound side and relative statistics *mape* on the underbound side.

    """
    from random import random
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    #import sapt_colors

    Nweft = len(labels)
    positions = range(-1, -1 * Nweft - 1, -1)

    # initialize plot
    fig, ax = plt.subplots(figsize=(12, 7))
    plt.xlim([-xlimit, xlimit])
    plt.ylim([-1 * Nweft - 1, 0])
    plt.yticks([])

    # label plot and tiers
    ax.text(-xlimit + 0.25, -0.25, title,
        verticalalignment='bottom', horizontalalignment='left',
        family='Times New Roman', weight='bold', fontsize=12)
    for weft in labels:
        ax.text(-xlimit + 0.25, -(1.2 + labels.index(weft)), weft,
            verticalalignment='bottom', horizontalalignment='left',
            family='Times New Roman', weight='bold', fontsize=18)

    # plot reaction errors and threads
    for rxn in data:
        xvals = rxn['data']
        toplblposn = next(item for item in xvals if item is not None)
        botlblposn = next(item for item in reversed(xvals) if item is not None)
        
#       global_color   rxn_color
##       None           None             grey
##                      colorname        colorname
##                      0-1 float        0-1 float
##       colorname       None           glob colorname
##                       colorname       glob colorname
##                       0-1 float       glob colorname
#       sapt           None            grey
#                       colorname       grey
#                       0-1 float       mpl.cm.jet
#       rgb             None           grey
#                       colorname       grey
#                       0-1 float       <.33 blue, >.67 red, else green

        # validate any sapt color
        try:
            rxn['color'] *= 1.0
        except (KeyError, TypeError):
            saptcolor = None
        else:
            if rxn['color'] >= 0.0 and rxn['color'] <= 1.0:
                saptcolor = rxn['color']
            else:
                saptcolor = None

        if color is None:
            # no color argument, so take from rxn
            if rxn['color'] is None:
                clr = 'grey'
            elif saptcolor is not None:
                clr = mpl.cm.jet(saptcolor)
            else:
                clr = rxn['color']
        elif color == 'sapt':
            # sapt color from rxn if available
            if saptcolor is not None:
                clr = mpl.cm.jet(saptcolor)
            else:
                clr = 'grey'
        elif color == 'rgb':
            # HB/MX/DD sapt color from rxn if available
            if saptcolor is not None:
                if saptcolor < 0.333:
                    clr = 'blue'
                elif saptcolor < 0.667:
                    clr = 'green'
                else:
                    clr = 'red'
            else:
                clr = 'grey'
        else:
            # color argument is name of mpl color
            clr = color

#        try:
#            clr = 'black' if color is None else mpl.cm.jet(sapt_colors.sapt_colors[rxn['dbse']][rxn['sys']])
#        except KeyError:
#            clr = 'black'
            
        ax.plot(xvals, positions, '-', color=clr)
        ax.plot(xvals, positions, '|', color=clr, markersize=10.0)

        ax.text(toplblposn, -0.75 + 0.6 * random(), rxn['sys'],
            verticalalignment='bottom', horizontalalignment='center',
            family='Times New Roman', fontsize=8)
        ax.text(botlblposn, -1 * Nweft - 0.75 + 0.6 * random(), rxn['sys'],
            verticalalignment='bottom', horizontalalignment='center',
            family='Times New Roman', fontsize=8)

    # plot trimmings
    if mae is not None:
        ax.plot([-x for x in mae], positions, 's', color='black')
    if mape is not None:  # equivalent to MAE for a 10 kcal/mol interaction energy
        ax.plot([0.025 * x for x in mape], positions, 'o', color='black')

    plt.axvline(0, color='black')
    plt.show()


if __name__ == "__main__":

    merge_dats = [
    {'dbse':'HSG', 'sys':'1', 'data':[0.3508, 0.1234, 0.0364, 0.0731, 0.0388]},
    {'dbse':'HSG', 'sys':'3', 'data':[0.2036, -0.0736, -0.1650, -0.1380, -0.1806]},
    {'dbse':'S22', 'sys':'14', 'data':[None, -3.2144, None, None, None]},
    {'dbse':'S22', 'sys':'15', 'data':[-1.5090, -2.5263, -2.9452, -2.8633, -3.1059]},
    {'dbse':'S22', 'sys':'22', 'data':[0.3046, -0.2632, -0.5070, -0.4925, -0.6359]}]

    thread(merge_dats, labels=['d', 't', 'dt', 'q', 'tq'], color='sapt',
        title='MP2-CPa[]z', mae=[0.25, 0.5, 0.5, 0.3, 1.0], mape=[20.1, 25, 15, 5.5, 3.6])
