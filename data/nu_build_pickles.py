import os
import sys
import time
#try:
import cPickle as pickle
#except ImportError:
#    import pickle
sys.path.append('/Users/loriab/linux/qcdb')
sys.path.append('/Users/loriab/linux/qcdb/databases')
#import qcdb
import qcdb.dbwrap

homewrite = '.'
dbnet = {}
#dbnet['S22'] = ['dft', 'saptone', 'pt2']
#dbnet['NBC10'] = ['dft', 'saptone', 'pt2']
#dbnet['HBC6'] = ['dft', 'saptone', 'pt2']
#dbnet['HSG'] = ['dft', 'saptone', 'pt2', 'bfdbmm']
#dbnet['SSI'] = ['saptmisc', 'bfdbmm', 'dfit', 'bfdbdft', 'pt2misc', 'ccmisc', 'dfitm', 'bfdbdftm', 'dftmisc']  #, 'efp']
#dbnet['BBI'] = ['saptmisc', 'bfdbmm', 'dfit', 'bfdbdft', 'pt2misc', 'ccmisc', 'dfitm', 'bfdbdftm', 'dftmisc']
#dbnet['PCONF'] = ['dfit', 'dfitm']
#dbnet['SCONF'] = ['dfit', 'dfitm']
#dbnet['ACONF'] = ['dfit', 'dfitm']
#dbnet['CYCONF'] = ['dfit', 'dfitm']
#dbnet['NBC10ext'] = ['saptmisc', 'dfit']
#dbnet['ACHC'] = ['saptmisc', 'dfit']
#dbnet['UBQ'] = ['saptmisc', 'bfdbmm']
#dbnet['S22by7'] = ['saptmisc']
#dbnet['S66'] = ['saptmisc']
#dbnet['A24'] = ['saptmisc', 'dilabio']
#dbnet['JSCH'] = ['saptmisc']
#dbnet[''] = []

# [1-4 Aug 2017] LAB
#   * added merz3 & 1hsg citations to HSG_bfdbmm.py
#   * added merz3 & 1ubq citations to UBQ_bfdbmm.py
#   * anon's allowed to remain b/c not dist. in qcdb conda: JSCH, S66, S22by7
#   * updated dfit, 1hsg, merz3, 1ubq, bfdbefp Citations
#     * so regen ACHC, *CONF, NBC10ext just in case
#   * added merz3 citations to BBI_bfdbdftm.py, BBI_bfdbdft.py (was bfdbdft), BBI_bfdbmm.py (some were bfdbmm), BBI_dftmisc.py, BBI_ccmisc.py, BBI_pt2misc.py, BBI_saptmisc.py
#   * added merz3 citations to SSI_bfdbdftm.py, SSI_bfdbdft.py (was bfdbdft), BBI_bfdbmm.py (some were bfdbmm), SSI_dftmisc.py, SSI_ccmisc.py, SSI_pt2misc.py, SSI_saptmisc.py
#   * added bfdbefp citation to SSI_efp.py
#   * added merz3 citation to UBQ_saptmisc.py as well as updating it to alpha0
#   * added cdsgroup (internal) citation to S66_saptmisc.py and JSCH_saptmisc.py
#   * added lots more citations & tagl so pretty much regen all
#   * added dfit citation to S22by7_saptmisc.py


#dbse = 'SSI'  # UNCOMMENT for local WDb

for db, lproj in dbnet.iteritems():
    print '\n<<< %s >>>' % (db)
    t0 = time.time()
    asdf = qcdb.dbwrap.WrappedDatabase(db)  # COMMENT for local WDb
    dbse = asdf.dbse  # COMMENT for local WDb
    t1 = time.time()
    print '%-70s %8.1f' % ('database.py --> WrappedDatabase', t1-t0)

    WDbfilename = homewrite + '/' + db + '_WDb.pickle'
    with open(WDbfilename, 'wb') as handle:  # COMMENT for local WDb
        pickle.dump(asdf, handle, pickle.HIGHEST_PROTOCOL)  # COMMENT for local WDb
    t2 = time.time()
    print '%-70s %8.1f' % ('* WrappedDatabase --> database.pickle', t2-t1)

    for pj in lproj:
        print '  * ' + pj
        t3 = time.time()
        with open(WDbfilename, 'rb') as handle:
            qwer = pickle.load(handle)    
        t4 = time.time()
        qwer.load_qcdata_byproject(pj)
        #qwer.load_qcdata_hdf5_trusted(project=pj)
        print '%-70s %8.1f' % ('    * database.pickle --> WrappedDatabase', t4-t3)
        t5 = time.time()
        hrxnfilename = homewrite + '/' + dbse + '_hrxn_' + pj + '.pickle'
        print '%-70s %8.1f' % ('    * WrappedDatabase --> WrappedDatabase + project', t5-t4)
        guts = {}
        for rxn, orxn in qwer.hrxn.iteritems():
            guts[rxn] = orxn.data
        with open(hrxnfilename, 'wb') as handle:
            #pickle.dump(qwer.hrxn, handle, pickle.HIGHEST_PROTOCOL)
            pickle.dump(guts, handle, pickle.HIGHEST_PROTOCOL)
        t6 = time.time()
        print '%-70s %8.1f' % ('    * WrappedDatabase+Project --> database_hrxn_project.pickle', t6-t5)

    print '  * Simulation'
    t7 = time.time()
    zxcv = qcdb.Database(db, loadfrompickle=True, path=homewrite)
    t8 = time.time()
    print '%-70s %8.1f' % ('    * database.pickle --> Database', t8-t7)
    for pj in lproj:
        t9 = time.time()
        zxcv.load_qcdata_hrxn_byproject(pj, path=homewrite)
        t10 = time.time()
        print '%-70s %8.1f' % ('    * ' + pj + 'database_hrxn_project.pickle --> Database + project', t10-t9)
    t11 = time.time()
    print '%-70s %8.1f' % ('  * pickles --> Database + all projects', t11-t7)
    nmc = len(zxcv.fancy_mcs().keys())
    t12 = time.time()
    print '%-70s %8.1f' % ('  * Access all projects ' + str(nmc), t12-t7)



