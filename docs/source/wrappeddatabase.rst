
WrappedDatabase as a Class
==========================

Member Data
-----------

.. autoclass:: qcdb.dbwrap.WrappedDatabase
   :members: dbse, hrxn, hrgt, sset

..    def load_pickled(dbname, path=None):
..    def available_modelchems(self, union=True):
..    def benchmark(self):

Loading Subsets
---------------

.. autoclass:: qcdb.dbwrap.WrappedDatabase
   :members: add_Subset, load_subsets

Loading QC Data 
---------------

.. autoclass:: qcdb.dbwrap.WrappedDatabase
   :members: add_ReactionDatum, load_qcdata, load_qcdata_byproject, load_qcdata_hdf5_trusted

Statistics
----------

.. autoclass:: qcdb.dbwrap.WrappedDatabase
   :members: compute_errors, compute_statistics


.. comment .. autoclass:: qcdb.dbwrap.Database
.. comment    :members: load_qcdata_byproject, plot_bars
.. comment 
.. comment    .. method:: load_qcdata_byproject
.. comment 
.. comment       Boil the noodle *time* minutes.
.. comment 
.. comment    .. attribute:: hrxn
.. comment 
.. comment       atusldfkjsld

