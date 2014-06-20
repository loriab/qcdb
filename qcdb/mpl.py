"""Module with matplotlib plotting routines. These are not hooked up to
any particular qcdb data structures but can be called with basic
arguments.

"""


def segment_color(argcolor, saptcolor):
    """Find appropriate color expression between overall color directive
    *argcolor* and particular color availibility *rxncolor*.

    """
    import matplotlib as mpl

    # validate any sapt color
    if saptcolor is not None:
        if saptcolor < 0.0 or saptcolor > 1.0:
            saptcolor = None

    if argcolor is None:
        # no color argument, so take from rxn
        if rxncolor is None:
            clr = 'grey'
        elif saptcolor is not None:
            clr = mpl.cm.jet(saptcolor)
        else:
            clr = rxncolor
    elif argcolor == 'sapt':
        # sapt color from rxn if available
        if saptcolor is not None:
            clr = mpl.cm.jet(saptcolor)
        else:
            clr = 'grey'
    elif argcolor == 'rgb':
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
        clr = argcolor

    return clr


def oldthread(data, labels, color=None, title='', xlimit=4.0, mae=None, mape=None):
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


def bar(data, title=''):
    #"""Generates a tiered slat diagram between model chemistries with
    #errors in list *data*, which is supplied as part of the dictionary for
    #each participating reaction, along with *dbse* and *rxn* keys in
    #argument *data*. The plot is labeled with *title* and each tier with
    #an element of *labels* and plotted at *xlimit* from the zero-line. If
    #*color* is None, slats are black, if 'sapt', colors are taken from
    #sapt_colors module. Summary statistics *mae* are plotted on the
    #overbound side and relative statistics *mape* on the underbound side.

    #"""
    from random import random
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    # initialize plot, fix dimensions for consistent Illustrator import
    fig, ax = plt.subplots(figsize=(12, 7))
    plt.ylim([0, 4.86])
    plt.xlim([0, 6])
    plt.xticks([])

    # label plot and tiers
    ax.text(0.4, 3.8, title,
        verticalalignment='bottom', horizontalalignment='left',
        family='Times New Roman', weight='bold', fontsize=12)

    widths = [0.15, 0.02, 0.02, 0.02]  # TT, HB, MX, DD
    xval = 0.1  # starting posn along x-axis

    # plot bar sets
    for rxn in data:
        lefts = [xval, xval + 0.025, xval + 0.065, xval + 0.105]

        rect = ax.bar(lefts, rxn['data'], widths, linewidth=0)
        rect[0].set_color('grey')
        rect[1].set_color('red')
        rect[2].set_color('green')
        rect[3].set_color('blue')

        ax.text(xval + .08, 3.5, rxn['mc'],
            verticalalignment='center', horizontalalignment='right', rotation='vertical',
            family='Times New Roman', fontsize=8)
        xval += 0.20

    plt.show()
    plt.savefig('mplbar_' + title + '.pdf', bbox_inches='tight', transparent=True, format='PDF')


def flat(data, color=None, title='', xlimit=4.0, mae=None, mape=None, view=True):
    """Generates a slat diagram between model chemistries with errors in
    single-item list *data*, which is supplied as part of the dictionary
    for each participating reaction, along with *dbse* and *rxn* keys in
    argument *data*. Limits of plot are *xlimit* from the zero-line. If
    *color* is None, slats are black, if 'sapt', colors are taken from
    sapt_colors module. Summary statistic *mae* is plotted on the
    overbound side and relative statistic *mape* on the underbound side.
    Saves a file with name *title* and plots to screen if *view*.

    """
    from random import random
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    Nweft = 1
    positions = range(-1, -1 * Nweft - 1, -1)

    # initialize plot
    fig, ax = plt.subplots(figsize=(12, 0.33))
    plt.xlim([-xlimit, xlimit])
    plt.ylim([-1 * Nweft - 1, 0])
    plt.yticks([])
    plt.xticks([])
    fig.patch.set_visible(False)
    ax.patch.set_visible(False)
    ax.axis('off')

    plt.axvline(-1.0, color='grey', linewidth=4)
    plt.axvline(-0.3, color='grey', linewidth=4)
    plt.axvline(0.0, color='grey', linewidth=4)
    plt.axvline(0.3, color='grey', linewidth=4)
    plt.axvline(1.0, color='grey', linewidth=4)

    # plot reaction errors and threads
    for rxn in data:
        xvals = rxn['data']

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

        ax.plot(xvals, positions, '|', color=clr, markersize=13.0)

    # plot trimmings
    if mae is not None:
        #ax.plot(-1 * mae, positions, 's', color='black', markersize=18.0)
        plt.axvline(-1 * mae, color='black', linewidth=12)
    if mape is not None:  # equivalent to MAE for a 10 kcal/mol interaction energy
        ax.plot(0.025 * mape, positions, 'o', color='black', markersize=15.0)

    plt.show()
    plt.savefig('mplflat_' + title + '.pdf', bbox_inches='tight', transparent=True, format='PDF')
    if not view:
        plt.close()


#def mpl_distslat_multiplot_files(pltfile, dbid, dbname, xmin, xmax, mcdats, labels, titles):
#    """Saves a plot with basename *pltfile* with a slat representation
#    of the modelchems errors in *mcdat*. Plot is in PNG, PDF, & EPS
#    and suitable for download, no mouseover properties. Both labeled
#    and labelless (for pub) figures are constructed.
#
#    """
#    import matplotlib as mpl
#    from matplotlib.axes import Subplot
#    import sapt_colors
#    from matplotlib.figure import Figure
#
#    nplots = len(mcdats)
#    fht = nplots * 0.8
#    fig, axt = plt.subplots(figsize=(12.0, fht))
#    plt.subplots_adjust(left=0.01, right=0.99, hspace=0.3)
#
#    axt.set_xticks([])
#    axt.set_yticks([])
#    plt.axis('off')
#
#    for item in range(nplots):
#        mcdat = mcdats[item]
#        label = labels[item]
#        title = titles[item]
#
#        erdat = np.array(mcdat)
#        yvals = np.ones(len(mcdat))
#        y = np.array([sapt_colors.sapt_colors[dbname][i] for i in label])
#
#        ax = Subplot(fig, nplots, 1, item + 1)
#        fig.add_subplot(ax)
#        sc = ax.scatter(erdat, yvals, c=y, s=3000, marker="|", cmap=mpl.cm.jet, vmin=0, vmax=1)
#
#        ax.set_yticks([])
#        ax.set_xticks([])
#        ax.set_frame_on(False)
#        ax.set_xlim([xmin, xmax])
#
#    # Write files with only slats
#    plt.savefig('scratch/' + pltfile + '_plain' + '.png', transparent=True, format='PNG')
#    plt.savefig('scratch/' + pltfile + '_plain' + '.pdf', transparent=True, format='PDF')
#    plt.savefig('scratch/' + pltfile + '_plain' + '.eps', transparent=True, format='EPS')
#
#    # Rewrite files with guides and labels
#    for item in range(nplots):
#        ax_again = fig.add_subplot(nplots, 1, item + 1)
#        ax_again.set_title(titles[item], fontsize=8)
#        ax_again.text(xmin + 0.3, 1.0, stats(np.array(mcdats[item])), fontsize=7, family='monospace', verticalalignment='center')
#        ax_again.plot([0, 0], [0.9, 1.1], color='#cccc00', lw=2)
#        ax_again.set_frame_on(False)
#        ax_again.set_yticks([])
#        ax_again.set_xticks([-12.0, -8.0, -4.0, -2.0, -1.0, 0.0, 1.0, 2.0, 4.0, 8.0, 12.0])
#        ax_again.tick_params(axis='both', which='major', labelbottom='off', bottom='off')
#    ax_again.set_xticks([-12.0, -8.0, -4.0, -2.0, -1.0, 0.0, 1.0, 2.0, 4.0, 8.0, 12.0])
#    ax_again.tick_params(axis='both', which='major', labelbottom='on', bottom='off')
#
#    plt.savefig('scratch/' + pltfile + '_trimd' + '.png', transparent=True, format='PNG')
#    plt.savefig('scratch/' + pltfile + '_trimd' + '.pdf', transparent=True, format='PDF')
#    plt.savefig('scratch/' + pltfile + '_trimd' + '.eps', transparent=True, format='EPS')


def thread(data, labels, color=None, title='', xlimit=4.0, mae=None, mape=None):
    """Generates a tiered slat diagram between model chemistries with
    errors (or simply values) in list *data*, which is supplied as part of the
    dictionary for each participating reaction, along with *dbse* and *rxn* keys
    in argument *data*. The plot is labeled with *title* and each tier with
    an element of *labels* and plotted at *xlimit* from the zero-line. If
    *color* is None, slats are black, if 'sapt', colors are taken from *color*
    key in *data* [0, 1]. Summary statistics *mae* are plotted on the
    overbound side and relative statistics *mape* on the underbound side.

    """
    from random import random
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    # initialize tiers/wefts
    Nweft = len(labels)
    lenS = 0.2
    gapT = 0.04
    positions = range(-1, -1 * Nweft - 1, -1)
    posnS = []
    for weft in range(Nweft):
        posnS.extend([positions[weft] + lenS, positions[weft] - lenS, None])
    posnT = []
    for weft in range(Nweft - 1):
        posnT.extend([positions[weft] - lenS - gapT, positions[weft + 1] + lenS + gapT, None])

    # initialize plot
    fht = Nweft * 0.8
    fig, ax = plt.subplots(figsize=(12, fht))
    plt.subplots_adjust(left=0.01, right=0.99, hspace=0.3)
    plt.xlim([-xlimit, xlimit])
    plt.ylim([-1 * Nweft - 1, 0])
    plt.yticks([])

    # label plot and tiers
    ax.text(-0.9 * xlimit, -0.25, title,
        verticalalignment='bottom', horizontalalignment='left',
        family='Times New Roman', weight='bold', fontsize=12)
    for weft in labels:
        ax.text(-0.9 * xlimit, -(1.2 + labels.index(weft)), weft,
            verticalalignment='bottom', horizontalalignment='left',
            family='Times New Roman', weight='bold', fontsize=18)

    # plot reaction errors and threads
    for rxn in data:
        
        # preparation
        xvals = rxn['data']
        clr = segment_color(color, rxn['color'] if 'color' in rxn else None)
        slat = []
        for weft in range(Nweft):
            slat.extend([xvals[weft], xvals[weft], None])
        thread = []
        for weft in range(Nweft - 1):
            thread.extend([xvals[weft], xvals[weft + 1], None])

        # plotting
        ax.plot(slat, posnS, color=clr, linewidth=1.0, solid_capstyle='round')
        ax.plot(thread, posnT, color=clr, linewidth=0.5, solid_capstyle='round',
            alpha=0.3)

        # labeling
        try:
            toplblposn = next(item for item in xvals if item is not None)
            botlblposn = next(item for item in reversed(xvals) if item is not None)
        except StopIteration:
            pass
        else:
            ax.text(toplblposn, -0.75 + 0.6 * random(), rxn['sys'],
                verticalalignment='bottom', horizontalalignment='center',
                family='Times New Roman', fontsize=8)
            ax.text(botlblposn, -1 * Nweft - 0.75 + 0.6 * random(), rxn['sys'],
                verticalalignment='bottom', horizontalalignment='center',
                family='Times New Roman', fontsize=8)

    # plot trimmings
    if mae is not None:
        ax.plot([-x for x in mae], positions, 's', color='black')
    if mape is not None:  # equivalent to MAE for a 10 kcal/mol IE
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

    more_dats = [
    {'mc':'MP2-CP-adz', 'data':[1.0, 0.8, 1.4, 1.6]},
    {'mc':'MP2-CP-adtz', 'data':[0.6, 0.2, 0.4, 0.6]}]

    bar(more_dats, title='asdf')

    single_dats = [
    {'dbse':'HSG', 'sys':'1', 'data':[0.3508]},
    {'dbse':'HSG', 'sys':'3', 'data':[0.2036]},
    {'dbse':'S22', 'sys':'14', 'data':[None]},
    {'dbse':'S22', 'sys':'15', 'data':[-1.5090]},
    {'dbse':'S22', 'sys':'22', 'data':[0.3046]}]

    #flat(single_dats, color='sapt', title='fg_MP2_adz', mae=0.25, mape=20.1)

    flat([{'sys': '1', 'color': 0.6933450559423702, 'data': [0.45730000000000004]}, {'sys': '2', 'color': 0.7627027688599753, 'data': [0.6231999999999998]}, {'sys': '3', 'color': 0.7579958735528617, 'data': [2.7624999999999993]}, {'sys': '4', 'color': 0.7560883254421639, 'data': [2.108600000000001]}, {'sys': '5', 'color': 0.7515161912065955, 'data': [2.2304999999999993]}, {'sys': '6', 'color': 0.7235223893438876, 'data': [1.3782000000000014]}, {'sys': '7', 'color': 0.7120099024225569, 'data': [1.9519000000000002]}, {'sys': '8', 'color': 0.13721565059144678, 'data': [0.13670000000000004]}, {'sys': '9', 'color': 0.3087395095814767, 'data': [0.2966]}, {'sys': '10', 'color': 0.25493207637105103, 'data': [-0.020199999999999996]}, {'sys': '11', 'color': 0.24093814608979347, 'data': [-1.5949999999999998]}, {'sys': '12', 'color': 0.3304746631959777, 'data': [-1.7422000000000004]}, {'sys': '13', 'color': 0.4156050644764822, 'data': [0.0011999999999989797]}, {'sys': '14', 'color': 0.2667207259626991, 'data': [-2.6083999999999996]}, {'sys': '15', 'color': 0.3767053567641695, 'data': [-1.5090000000000003]}, {'sys': '16', 'color': 0.5572641509433963, 'data': [0.10749999999999993]}, {'sys': '17', 'color': 0.4788598239641578, 'data': [0.29669999999999996]}, {'sys': '18', 'color': 0.3799031371351281, 'data': [0.10209999999999964]}, {'sys': '19', 'color': 0.5053227185999078, 'data': [0.16610000000000014]}, {'sys': '20', 'color': 0.2967660584483015, 'data': [-0.37739999999999974]}, {'sys': '21', 'color': 0.38836460733750316, 'data': [-0.4712000000000005]}, {'sys': '22', 'color': 0.5585849893078809, 'data': [0.30460000000000065]}, {'sys': 'BzBz_PD36-1.8', 'color': 0.1383351040559965, 'data': [-1.1921]}, {'sys': 'BzBz_PD34-2.0', 'color': 0.23086034843049832, 'data': [-1.367]}, {'sys': 'BzBz_T-5.2', 'color': 0.254318060864096, 'data': [-0.32230000000000025]}, {'sys': 'BzBz_T-5.1', 'color': 0.26598486566733337, 'data': [-0.3428]}, {'sys': 'BzBz_T-5.0', 'color': 0.28011258347610224, 'data': [-0.36060000000000025]}, {'sys': 'PyPy_S2-3.9', 'color': 0.14520332101084785, 'data': [-0.9853000000000001]}, {'sys': 'PyPy_S2-3.8', 'color': 0.1690757103699542, 'data': [-1.0932]}, {'sys': 'PyPy_S2-3.5', 'color': 0.25615734567417053, 'data': [-1.4617]}, {'sys': 'PyPy_S2-3.7', 'color': 0.19566550224566906, 'data': [-1.2103999999999995]}, {'sys': 'PyPy_S2-3.6', 'color': 0.22476748600170826, 'data': [-1.3333]}, {'sys': 'BzBz_PD32-2.0', 'color': 0.31605681987208084, 'data': [-1.6637]}, {'sys': 'BzBz_T-4.8', 'color': 0.31533827331543723, 'data': [-0.38759999999999994]}, {'sys': 'BzBz_T-4.9', 'color': 0.2966146678069063, 'data': [-0.3759999999999999]}, {'sys': 'BzH2S-3.6', 'color': 0.38284814928043304, 'data': [-0.1886000000000001]}, {'sys': 'BzBz_PD32-1.7', 'color': 0.3128835191478639, 'data': [-1.8703999999999998]}, {'sys': 'BzMe-3.8', 'color': 0.24117892478245323, 'data': [-0.034399999999999986]}, {'sys': 'BzMe-3.9', 'color': 0.22230903086047088, 'data': [-0.046499999999999986]}, {'sys': 'BzH2S-3.7', 'color': 0.36724255203373696, 'data': [-0.21039999999999992]}, {'sys': 'BzMe-3.6', 'color': 0.284901522674611, 'data': [0.007099999999999884]}, {'sys': 'BzMe-3.7', 'color': 0.2621086166558813, 'data': [-0.01770000000000005]}, {'sys': 'BzBz_PD32-1.9', 'color': 0.314711251903219, 'data': [-1.7353999999999998]}, {'sys': 'BzBz_PD32-1.8', 'color': 0.3136181753200793, 'data': [-1.8039999999999998]}, {'sys': 'BzH2S-3.8', 'color': 0.3542001591399945, 'data': [-0.22230000000000016]}, {'sys': 'BzBz_PD36-1.9', 'color': 0.14128552184232473, 'data': [-1.1517]}, {'sys': 'BzBz_S-3.7', 'color': 0.08862098445220466, 'data': [-1.3414]}, {'sys': 'BzH2S-4.0', 'color': 0.33637540012259076, 'data': [-0.2265999999999999]}, {'sys': 'BzBz_PD36-1.5', 'color': 0.13203548045236127, 'data': [-1.3035]}, {'sys': 'BzBz_S-3.8', 'color': 0.0335358832178858, 'data': [-1.2022]}, {'sys': 'BzBz_S-3.9', 'color': 0.021704594689389095, 'data': [-1.0747]}, {'sys': 'PyPy_T3-5.1', 'color': 0.3207725129126432, 'data': [-0.2958000000000003]}, {'sys': 'PyPy_T3-5.0', 'color': 0.3254925304351165, 'data': [-0.30710000000000015]}, {'sys': 'BzBz_PD36-1.7', 'color': 0.13577087141986593, 'data': [-1.2333000000000003]}, {'sys': 'PyPy_T3-4.8', 'color': 0.3443704059902452, 'data': [-0.32010000000000005]}, {'sys': 'PyPy_T3-4.9', 'color': 0.3333442013628509, 'data': [-0.3158999999999996]}, {'sys': 'PyPy_T3-4.7', 'color': 0.35854000505665756, 'data': [-0.31530000000000014]}, {'sys': 'BzBz_PD36-1.6', 'color': 0.13364651314909243, 'data': [-1.2705000000000002]}, {'sys': 'BzMe-4.0', 'color': 0.20560117919562013, 'data': [-0.05389999999999984]}, {'sys': 'MeMe-3.6', 'color': 0.16934865900383142, 'data': [0.18420000000000003]}, {'sys': 'MeMe-3.7', 'color': 0.1422332591197123, 'data': [0.14680000000000004]}, {'sys': 'MeMe-3.4', 'color': 0.23032794290360467, 'data': [0.29279999999999995]}, {'sys': 'MeMe-3.5', 'color': 0.19879551978386897, 'data': [0.23260000000000003]}, {'sys': 'MeMe-3.8', 'color': 0.11744404936205816, 'data': [0.11680000000000001]}, {'sys': 'BzBz_PD34-1.7', 'color': 0.22537382457222138, 'data': [-1.5286999999999997]}, {'sys': 'BzBz_PD34-1.6', 'color': 0.22434088042760192, 'data': [-1.5754000000000001]}, {'sys': 'BzBz_PD32-2.2', 'color': 0.3189891685300601, 'data': [-1.5093999999999999]}, {'sys': 'BzBz_S-4.1', 'color': 0.10884135031532088, 'data': [-0.8547000000000002]}, {'sys': 'BzBz_S-4.0', 'color': 0.06911476296747143, 'data': [-0.9590000000000001]}, {'sys': 'BzBz_PD34-1.8', 'color': 0.22685419834431494, 'data': [-1.476]}, {'sys': 'BzBz_PD34-1.9', 'color': 0.2287079261672095, 'data': [-1.4223999999999997]}, {'sys': 'BzH2S-3.9', 'color': 0.3439077006047999, 'data': [-0.22739999999999982]}, {'sys': 'FaNNFaNN-4.1', 'color': 0.7512716174974567, 'data': [1.7188999999999997]}, {'sys': 'FaNNFaNN-4.0', 'color': 0.7531388297328865, 'data': [1.9555000000000007]}, {'sys': 'FaNNFaNN-4.3', 'color': 0.7478064149182957, 'data': [1.2514000000000003]}, {'sys': 'FaNNFaNN-4.2', 'color': 0.7493794908838113, 'data': [1.4758000000000013]}, {'sys': 'FaOOFaON-4.0', 'color': 0.7589275618320565, 'data': [2.0586]}, {'sys': 'FaOOFaON-3.7', 'color': 0.7619465815742713, 'data': [3.3492999999999995]}, {'sys': 'FaOOFaON-3.9', 'color': 0.7593958895631474, 'data': [2.4471000000000007]}, {'sys': 'FaOOFaON-3.8', 'color': 0.7605108059280967, 'data': [2.8793999999999986]}, {'sys': 'FaONFaON-4.1', 'color': 0.7577459277014137, 'data': [1.8697999999999997]}, {'sys': 'FaOOFaON-3.6', 'color': 0.7633298028299997, 'data': [3.847599999999998]}, {'sys': 'FaNNFaNN-3.9', 'color': 0.7548200901251662, 'data': [2.2089]}, {'sys': 'FaONFaON-3.8', 'color': 0.7582294603551467, 'data': [2.967699999999999]}, {'sys': 'FaONFaON-3.9', 'color': 0.7575285282217349, 'data': [2.578900000000001]}, {'sys': 'FaONFaON-4.2', 'color': 0.7594549221042256, 'data': [1.5579999999999998]}, {'sys': 'FaOOFaNN-3.6', 'color': 0.7661655616885379, 'data': [3.701599999999999]}, {'sys': 'FaOOFaNN-3.7', 'color': 0.7671068376007428, 'data': [3.156500000000001]}, {'sys': 'FaOOFaNN-3.8', 'color': 0.766947626251711, 'data': [2.720700000000001]}, {'sys': 'FaONFaNN-3.9', 'color': 0.7569836601896789, 'data': [2.4281000000000006]}, {'sys': 'FaONFaNN-3.8', 'color': 0.758024548462959, 'data': [2.7561999999999998]}, {'sys': 'FaOOFaOO-3.6', 'color': 0.7623422640217077, 'data': [3.851800000000001]}, {'sys': 'FaOOFaOO-3.7', 'color': 0.7597430792159379, 'data': [3.2754999999999974]}, {'sys': 'FaOOFaOO-3.4', 'color': 0.7672554950739594, 'data': [5.193299999999999]}, {'sys': 'FaOOFaOO-3.5', 'color': 0.764908813123865, 'data': [4.491900000000001]}, {'sys': 'FaONFaNN-4.2', 'color': 0.7549212942233738, 'data': [1.534699999999999]}, {'sys': 'FaONFaNN-4.0', 'color': 0.7559404310956357, 'data': [2.1133000000000024]}, {'sys': 'FaONFaNN-4.1', 'color': 0.7551574698775625, 'data': [1.813900000000002]}, {'sys': 'FaONFaON-4.0', 'color': 0.7572064604483282, 'data': [2.2113999999999994]}, {'sys': 'FaOOFaOO-3.8', 'color': 0.7573810956831686, 'data': [2.7634000000000007]}, {'sys': '1', 'color': 0.2784121805328983, 'data': [0.3508]}, {'sys': '2', 'color': 0.22013842798900166, 'data': [-0.034600000000000186]}, {'sys': '3', 'color': 0.12832496088281312, 'data': [0.20360000000000023]}, {'sys': '4', 'color': 0.6993695033529733, 'data': [1.9092000000000002]}, {'sys': '5', 'color': 0.7371192790053749, 'data': [1.656600000000001]}, {'sys': '6', 'color': 0.5367033190796172, 'data': [0.27970000000000006]}, {'sys': '7', 'color': 0.3014220615964802, 'data': [0.32289999999999974]}, {'sys': '8', 'color': 0.01605867807629261, 'data': [0.12199999999999994]}, {'sys': '9', 'color': 0.6106300539083558, 'data': [0.3075999999999999]}, {'sys': '10', 'color': 0.6146680031333968, 'data': [0.6436000000000002]}, {'sys': '11', 'color': 0.6139747851721759, 'data': [0.4551999999999996]}, {'sys': '12', 'color': 0.32122739401126593, 'data': [0.44260000000000005]}, {'sys': '13', 'color': 0.24678148099136055, 'data': [-0.11789999999999967]}, {'sys': '14', 'color': 0.23700950710597016, 'data': [0.42689999999999995]}, {'sys': '15', 'color': 0.23103396678138563, 'data': [0.3266]}, {'sys': '16', 'color': 0.1922070769654413, 'data': [0.0696000000000001]}, {'sys': '17', 'color': 0.19082151944747366, 'data': [0.11159999999999992]}, {'sys': '18', 'color': 0.2886200282444196, 'data': [0.4114]}, {'sys': '19', 'color': 0.23560171133945224, 'data': [-0.1392]}, {'sys': '20', 'color': 0.3268270751294533, 'data': [0.5593]}, {'sys': '21', 'color': 0.7324460869158442, 'data': [0.6806000000000001]}],
    color='sapt', title='MP2-CP-adz', mae=1.21356003247, mape=24.6665886087, xlimit=4.0)
