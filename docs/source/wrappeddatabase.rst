
Database as a Class
===================

Member Data
-----------

.. autoclass:: qcdb.dbwrap.Database
   :members: dbse, hrxn, hrgt, sset

Loading Subsets
---------------

.. autoclass:: qcdb.dbwrap.Database
   :members: add_Subset, load_subsets

Loading QC Data 
---------------

.. autoclass:: qcdb.dbwrap.Database
   :members: add_ReactionDatum, load_qcdata, load_qcdata_byproject

Statistics
----------

.. autoclass:: qcdb.dbwrap.Database
   :members: compute_errors, compute_statistics, analyze_modelchems

Plots
-----

.. autoclass:: qcdb.dbwrap.Database
   :members: plot_flat, plot_modelchems, plot_bars, plot_iowa

Tables
------

.. autoclass:: qcdb.dbwrap.Database
   :members: table_generic, table_simple1

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

