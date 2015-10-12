import os
import sys
import time
#try:
import cPickle as pickle
#except ImportError:
#    import pickle
sys.path.append('/Users/loriab/linux/qcdb')
sys.path.append('/Users/loriab/linux/qcdb/databases')
import qcdb
#import qcdb.dbwrap

homewrite = '.'
dbnet = {}
#dbnet['S22'] = ['dft', 'saptone', 'pt2']
#dbnet['NBC10'] = ['dft', 'saptone', 'pt2']
#dbnet['HBC6'] = ['dft', 'saptone', 'pt2']
#dbnet['HSG'] = ['dft', 'saptone', 'pt2']
#dbnet['SSI'] = ['saptmisc', 'bfdbmm', 'merz3', 'dfit', 'bfdbdft']
#dbnet['BBI'] = ['saptmisc', 'bfdbmm', 'merz3', 'dfit', 'bfdbdft']
#dbnet['PCONF'] = ['dfit']
#dbnet['SCONF'] = ['dfit']
#dbnet['ACONF'] = ['dfit']
#dbnet['CYCONF'] = ['dfit']
#dbnet['DB4'] = ['dft', 'saptone', 'pt2']
#dbnet['NBC10ext'] = ['dfit']
#dbnet['ACHC'] = ['saptmisc', 'dfit']
#dbnet['UBQ'] = ['saptmisc', 'bfdbmm']
#dbnet['S22by7'] = ['saptmisc']
#dbnet['S66'] = ['saptmisc']
#dbnet['A24'] = ['dilabio']
#dbnet['JSCH'] = ['saptmisc']
#dbnet[''] = []

for db, lproj in dbnet.iteritems():
    print '\n<<< %s >>>' % (db)
    t0 = time.time()
    if db == 'DB4':
        asdf = qcdb.DB4(loadfrompickle=True, path=homewrite)
    else:
        asdf = qcdb.Database(db, loadfrompickle=True, path=homewrite)
    dbse = asdf.dbse
    t1 = time.time()
    print '%-70s %8.1f' % ('database.py --> Database', t1-t0)
    print 'Benchmark', asdf.benchmark

    Dbfilename = homewrite + '/' + db + '_Db.pickle'
    with open(Dbfilename, 'wb') as handle:
        pickle.dump(asdf, handle, pickle.HIGHEST_PROTOCOL)
    t2 = time.time()
    print '%-70s %8.1f' % ('* Database --> database_Db.pickle', t2-t1)

    for pj in lproj:
        t3 = time.time()
        asdf.load_qcdata_byproject(pj)
        t4 = time.time()
        print '%-70s %8.1f' % ('    * Database --> Database + project %s' % (pj), t4-t3)

    t5 = time.time()
    DbLfilename = homewrite + '/' + db + '_DbLoaded.pickle'
    with open(DbLfilename, 'wb') as handle:
        pickle.dump(asdf, handle, pickle.HIGHEST_PROTOCOL)
    t6 = time.time()
    print '%-70s %8.1f' % ('* Database + projects --> database_DbLoaded.pickle', t6-t5)


    print '  * Simulation'
    t7 = time.time()
    with open(Dbfilename, 'rb') as handle:
        zxcv = pickle.load(handle)
    t8 = time.time()
    print '%-70s %8.1f' % ('    * database_Db.pickle --> Database', t8-t7)

    t11 = time.time()
    with open(DbLfilename, 'rb') as handle:
        zxcvL = pickle.load(handle)
    t12 = time.time()
    print '%-70s %8.1f' % ('  * database_DbLoaded.pickle --> Database + all projects', t12-t11)


    nmc = len(zxcv.fancy_mcs().keys())
    print zxcv.sset.keys()
    print '%-70s %8.1f' % ('  * Access Database ' + str(nmc), t8-t7)
    nmc = len(zxcvL.fancy_mcs().keys())
    print zxcvL.sset.keys()
    print '%-70s %8.1f' % ('  * Access all projects ' + str(nmc), t12-t11)



