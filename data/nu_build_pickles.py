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

# S22, A24, S22by7, BBI, NBC10, NBC10ext, S66, HSG, HBC6, SSI
homewrite = '.'
dbnet = {}
dbnet['S22'] = ['dft', 'saptone', 'pt2']
dbnet['NBC10'] = ['dft', 'saptone', 'pt2']
dbnet['HBC6'] = ['dft', 'saptone', 'pt2']
dbnet['HSG'] = ['dft', 'saptone', 'pt2']
dbnet['SSI'] = ['bfdbmm', 'merz3', 'dfit']
dbnet['PCONF'] = ['dfit']

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



