from __future__ import absolute_import
import sys
sys.path.append('/theoryfs2/ds/cdsgroup/qcdb')
sys.path.append('/theoryfs2/ds/cdsgroup/qcdb/qcdb')
import qcdb

import numpy as np
from . import mpl
import matplotlib.pyplot as plt

def plot_convergence(rxn, dbse, proj, loadfrompickle, std_mtd, std_bsse, std_bas_list, std_curvestyle, mtd_list, bsse, bas_list, curvestyle,
                    failoninc = True, assigned_title='default', xlabel='default', ylabel='default', legend=True, 
                    markersize='default', linewidth='default',
                    standalone=True, plotpath='default', filename='default', filetype='default'):

    """Prepares convergence plots for all members of array *rxn* of
    reactions from database *dbse* in project(s) *proj* of modelchemistry
    defined by all combinations of the method list, *mtd_list*,
    counterpoise scheme, *bsse*, and basis set list, *bas_list*, arrays.

    Convergence curves are plotted relative to reference value, with
    style defined in array *curvestyle* for each of the curves to be
    created.  Curve styles are three-character strings giving marker
    style, curve color, and line style.

    Plot title *title* will default to a standard string identifier
    based on *rxn* and *dbse*, unless otherwise specified. Axis labels
    *xlabel* and *ylabel* may be specified, but have default values
    "Basis Set" and "Interaction Energy" otherwise. Ticks on the x axis
    will be basis sets defined in *bas_list*.

    Inclusion of legend may be turned on or off by the truth value
    of *legend*.  Marker and line size may be changed by adjusting
    *markersize* and *linewidth* arguments, both accepting float values.
    Default linewidth is 1, default marker size is 5.5.

    Plots are opened if *standalone* is true, in .pdf format by
    default. File type can be changed by *filetype*, and will be saved in
    location designated by *plotpath* as *dbse*_*rxn*_convplot.*filetype*
    if *filename* is left as 'default'.

    """
    # load data from corresponding function arguments
    asdf = qcdb.Database([dbse], loadfrompickle=loadfrompickle)
    for pr in proj:
            asdf.load_qcdata_byproject(pr)
    
    # Build 2D modelchems array for novel methods
    mcs = []
    for mtd in mtd_list:
        for b in bsse:
            temp_mcs_row = []
            for bas in bas_list:
                temp_bas = []
                temp_bas.append(mtd + '-' + b + '-' + bas)
                temp_mcs_row.append(temp_bas[0])
        mcs.append(temp_mcs_row)

    # Build 1D standard modelchem array
    std_mcs = []
    for std_bas in std_bas_list:
        std_mcs.append(std_mtd + '-' + std_bsse + '-' + std_bas)

    # Get benchmark values for all rxns/modelchems in database
    bench = asdf.benchmark # Where to look for the benchmark values
    temp_ref_array = []
    for lmc, lbm, orxn in asdf.get_reactions(modelchem=bench): 
        temp_ref_array.append(orxn.data[bench].value) # Get all benchmark values for database
    
    # TODO: Form rxn array based on more intuitive string argument?

    # Iterate over array *rxn* and create plot for each one
    for r in rxn:
        # Assign plot title
        
        if assigned_title == 'default':
            working_title = dbse + '-' + str(r) + ' ConvPlot'
        else:
            working_title = assigned_title
        # Get benchmark value for r
        ref = temp_ref_array[r-1]
        # Build rank 1 *std_dbdata* array containing the data for the standard convergence curve
        std_dbdata = []
        for std_mc in std_mcs:
            all_std_mc_data = []
            if std_mc in asdf.mcs.keys():
                for lmc, lbm, orxn in asdf.get_reactions(modelchem=std_mc, failoninc=failoninc):
                    try:
                       all_std_mc_data.append(orxn.data[std_mc].value)
                    except KeyError, e:
                        if failoninc:
                            raise e
                        else:
                            all_std_mc_data.append(np.NaN)
            else:
                print "Standard model chemistry not in mcs.keys(). Try again."
            std_dbdata.append(all_std_mc_data[r-1])
        std_curve_label=std_mtd + '-' + std_bsse
        std_x_vals = np.arange(1, len(std_dbdata) + 1)[:]
        # Build rank 2 *dbdata* array containing all data for rxn *r* to be plotted 
        # Build rank 1 *curve_labels* array for use in legend
        # Build rank 1 *x_vals* array against which to plot rows of *dbdata*
        dbdata = []
        x_vals = []
        curve_labels = []
        for row in mcs: # Each row contains different combo of *mcs* and *bsse*
            temp_dbdata_row = []
            for mc in row: # Iterate over different values of *bas* for same *mcs* & *bsse*
                all_dbse_mc_data = []
                if mc in asdf.mcs.keys(): # If the modelchem exists at all, regardless of some missing data
                    for lmc, lbm, orxn in asdf.get_reactions(modelchem=mc,failoninc=failoninc): # Gets rxndata dictionary 
                        try:
                            all_dbse_mc_data.append(orxn.data[mc].value) # Appends all *dbse*/*mc* rxn datum to temporary array
                        except KeyError, e: # Protects against missing data
                            if failoninc: 
                                raise e # Raises offending mc label
                            else:
                                all_dbse_mc_data.append(np.NaN) # Appends None to dbdata array, so nothing gets plotted for no data
                else:
                    for i in range(24): # Need to append 24 NaN's to simulate totally empty array
                        all_dbse_mc_data.append(np.NaN) # Append None to dbdata array for an entire mc if it doesn't exist
                temp_dbdata_row.append(all_dbse_mc_data[r-1]) # Appends rxn *r*/*mc* datum to array for later append to *dbdata*
                # Split modelchem *mc* name *mtd*-*bsse*-*bas* by '-' delimiter, save in temp array 
                temp_curve_label = mc.split("-")
            # Concatenate *mtd* + '-' + *bsse* to form unique curve_label, add to *curve_labels* array
            curve_labels.append(temp_curve_label[0] + '-' + temp_curve_label[1])
            dbdata.append(temp_dbdata_row) # Add row to *dbdata* array
        x_vals = np.arange(1, len(bas_list) + 1)[:] # Evelnly spaced x values for every basis set in *bas*
        # Set default *linewidth* and *markersize*
        if linewidth=='default':
            linewidth = 1
        if markersize=='default':
            markersize = 5.5 
        # Regardless of this decision, assume from here on that *curvestyle* is an array.
        # Iterate over curves (by each label in *curve_labels*) and plot x_vals, dbdata[row], curvestyle, curve_labels
        for l in range(0, len(curve_labels)):
            plt.plot(x_vals, dbdata[l], curvestyle[l], label= curve_labels[l],linewidth=linewidth, markersize=markersize)
        # Plot standard curve
        plt.plot(std_x_vals, std_dbdata, std_curvestyle, label=std_curve_label, linewidth=linewidth, markersize=markersize)
        # Set other plot options
        plt.axhline(ref, color='black', linestyle='--')
        plt.title(working_title)
        if xlabel=='default':
            plt.xlabel('Basis Set')
        else: 
            plt.xlabel(xlabel)
        if ylabel=='default':
            plt.ylabel('Interaction Energy')
        else:
            plt.ylabel(ylabel)
        if legend:
            plt.legend(loc = 'upper right') 
        # Make basis labels array
        bas_labels = []
        if len(bas_list) < len(std_bas_list):
            for bas in bas_list:
                bas_labels.append(bas)
            for std_bas in range(len(bas_list), len(std_bas_list)):
                bas_labels.append(std_bas_list[std_bas])
        else:
                bas_labels = bas_list
        plt.xticks(np.arange(1, len(std_bas_list) + 1), bas_labels)
        plt.margins(0.2)
        # Save or display plot
        if standalone:  
            plt.show()
        elif plotpath == 'default' and filename == 'default' and filetype == 'default':
            # Form a unique filename, save in working directory as .pdf
            plt.savefig(dbse + '-' + str(r) + '-convplot.pdf')
        elif plotpath == 'default' and filename == 'defulat':
            plt.savefig(dbse + '-' + str(r) + '-convplot' + filetype)
        elif plotpath == 'defulat':
            plt.savefig(filename + '.' + filetype)
        else:
            plt.savefig( plotpath + '/' + filename + '.' + filetype)
        del working_title
        plt.close()
        # Delete the values from curve_lables and ref so it doesn't re-plot
        del curve_labels
        del ref
        del dbdata
        del std_dbdata

plot_convergence(
rxn = range(1,25),
dbse = 'A24',
proj = ['f12dilabio','dilabio'],
loadfrompickle=True,
std_mtd = 'CCSDT',
std_bsse = 'CP',
std_bas_list = ['adz','atz','aqz','a5z','a6z'],
std_curvestyle = 'ro-',
mtd_list = ['CCSDTAF12','CCSDTBF12','CCSDTCF12'],
bsse = ['CP'],
bas_list = ['adz','atz','aqz','a5z'],
curvestyle=['bs-','g^-','m*-'],
failoninc=False,
assigned_title='default',
xlabel='default',
ylabel='default',
legend=True,
markersize=7,
linewidth='default',
standalone=False,
plotpath='default',
filename='default',
filetype='default'
)
