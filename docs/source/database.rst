
Interactive Database Objects
============================

Member Data
-----------

.. autoclass:: qcdb.dbwrap.Database
   :members: dbse, dbdict, sset, mcs

..    def intersect_subsets(self):
..    def intersect_modelchems(self):


Loading Subsets
---------------

.. autoclass:: qcdb.dbwrap.Database
   :members: load_subsets

Loading QC Data 
---------------

.. autoclass:: qcdb.dbwrap.Database
   :members: load_qcdata, load_qcdata_byproject, load_qcdata_hdf5_trusted

Statistics
----------

.. autoclass:: qcdb.dbwrap.Database
   :members: compute_errors, compute_statistics, analyze_modelchems

Convenience Functions
---------------------

.. autoclass:: qcdb.dbwrap.Database
   :members: fancy_mcs

Plots
-----

.. autoclass:: qcdb.dbwrap.Database
   :members: plot_bars, plot_disthist, plot_flat, plot_all_flats, plot_modelchems, plot_iowa

Tables
------

.. autoclass:: qcdb.dbwrap.Database
   :members: table_generic, table_simple1, table_merge_abbr, table_merge_suppmat

