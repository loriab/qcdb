import os
import sys
import time
try:
    import cPickle as pickle
except ImportError:
    import pickle
sys.path.append('/Users/loriab/linux/qcdb')
sys.path.append('/Users/loriab/linux/qcdb/databases')
#import qcdb
import qcdb.dbwrap

databases = ['S22', 'NBC10', 'HBC6', 'HSG']
for db in databases:
    print '\n<<< %s >>>' % (db)
    t0 = time.time()
    #asdf = qcdb.Database(db)
    asdf = qcdb.dbwrap.WrappedDatabase(db)
    t1 = time.time()
    print asdf.sset.keys()

    filename = db + '.pickle'
    asdf.load_qcdata_byproject('dft')
    asdf.load_qcdata_byproject('saptone')
    asdf.load_qcdata_byproject('pt2')
    #asdf.load_qcdata_byproject('sflow')
    #asdf.load_qcdata_hdf5_trusted(project='dft') #, path='/Users/loriab/linux/qcdb/sandbox/bfdb')
    t4 = time.time()
    
    with open(filename, 'wb') as handle:
      pickle.dump(asdf, handle)
    t5 = time.time()
    
    with open(filename, 'rb') as handle:
      asdf3 = pickle.load(handle)
    t6 = time.time()
    
    print '%-40s %8.1f' % ('form Database from db.py', t1-t0)
    #print '%-40s %8.1f' % ('load qcdata from hdf5', t4-t1)
    print '%-40s %8.1f' % ('load qcdata from ReactionDatums', t4-t1)
    print '%-40s %8.1f' % ('write pickled loaded Db', t5-t4)
    print '%-40s %8.1f' % ('read pickled loaded Db', t6-t5)

