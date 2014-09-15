
Obtaining |qcdb| Locally
========================

There is a |qcdb| repository on Lori's GitHub account, but the active
repository is stored at https://github.com/cdsgroup/qcdb , since it contains unpublished research
results.  An ordinary git clone operation on Mac or Windows (with GitHub
for Windows) suffices to acquire the repository. ::

    git clone https://github.com/cdsgroup/qcdb.git

The above is sufficient to run |qcdb| in Canopy or to run |qcdb| scripts.
If you want to be able to simply ``import qcdb`` in a python file or to
use the QC interface capabilities in |PSIfour| that need an advanced
version of the |qcdb| module (an old, stable version is included in the
|PSIfour| repository), the :envvar:`PYTHONPATH` needs to find |qcdb|.
Note the directory into which you placed the |qcdb| repository. Assuming
you have a directory structure such that the following is a file
``/Users/johndoe/linux/qcdb/qcdb/__init__.py``, execute something like the
following. ::

    export PYTHONPATH=/Users/johndoe/linux/qcdb:$PYTHONPATH

.. note :: Be aware that having |qcdb| in your :envvar:`PYTHONPATH` can 
   interfere with |PSIfour| since the latter imports an old version of 
   |qcdb| as part of the driver workings. In this case, it is recommended 
   that whatever script you have to set :envvar:`PATH`, :envvar:`PSIDATADIR` 
   in preparation to running |PSIfour| be modified to use a restricted 
   :envvar:`PYTHONPATH`.

Initialization in Canopy
========================

A blank Canopy session needs to be instructed where to find |qcdb| and to
build an object containing all the data you want to analyze. It is easiest
to just keep around a little script that defines these initialization
steps. ::

    import sys
    sys.path.append('/Users/johndoe/linux/qcdb')
    sys.path.append('/Users/johndoe/linux/qcdb/databases')
    import qcdb

    # to work with S22, NBC10, HBC6, HSG together
    asdf = qcdb.FourDatabases()
    asdf.load_qcdata_byproject('saptone')

    # to work with a single database
    asdf = qcdb.Database('S22')
    asdf.load_qcdata_byproject('dft')

.. note:: Any time you alter the qcdb module code, the python kernal in 
   Canopy must be restarted and |qcdb| re-imported for those changes to 
   become apparent.

