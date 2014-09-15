
Databases, Meeting and Acquiring
================================

The fundamental analysis unit of |qcdb| is the database, which is simply a
collection of geometries (reagents; rgt), their algebraic arrangement into
chemically interesting quantities (reactions; rxn), and a reference value
for the reactions (optional). Read
http://sirius.chem.vt.edu/psi4manual/master/quickadddatabase.html for
further details.

If you have a project with xyz files and no database, see
:source:`bin/ixyz2database.py` for a script to build one. Place it
in ``qcdb/databases``. No further registry required. And **do** fill in
the docstring at the top with description, citations, clues, etc. about
the origin of the database content; you and others will later be glad of
it.

Technical Details
-----------------

A database is a very plain file format, so as to not intimidate those
unfamiliar with python. Nevertheless, it is a proper python module and
importable as such (a one-file python module, so no ``__init__.py``) with
``import S22`` or the like. Below is code to import a generic database
*db_name*. ::

    try:
        database = __import__(db_name)
    except ImportError:
        print('\nPython module for database %s failed to load\n\n' % (db_name))
        print('\nSearch path that was tried:\n')
        print(", ".join(map(str, sys.path)))
        raise ValidationError("Python module loading problem for database " + str(db_name))
    else:
        dbse = database.dbse
        HRXN = database.HRXN
        ACTV = getattr(database, mode)
        RXNM = database.RXNM
        BIND = database.BIND
        TAGL = database.TAGL
        GEOS = database.GEOS
        try:
            DATA = database.DATA
        except AttributeError:
            DATA = {}

As is apparent, once imported, the structures are available to you.

.. py:data:: dbse 

   String of the internal name of the database. Usually the same as the file
   name, but can be different as for :source:`databases/HBC6.py` with ``HBC1``. 
   It is with this name that :source:`bin/herd-DB.py` will seek output 
   files for the database.

.. py:data:: GEOS

   Dictionary of database geometries. Contains a key-value pair for every 
   reagent where the key is dbse-rgt and the value is a :py:class:`qcdb.Molecule` 
   object. With the power of the |PSIfour| Molecule class (translated to python),
   reagents can be defined in terms of one another.

.. py:data:: HRXN

   Array of indices of all reactions in the database. This defines the 
   ordering of reactions in the database. For numerical indices, *e.g.* 
   S22-5, choose integers or strings and stick with it.

.. py:data:: ACTV

   Dictionary defining the ordering of reagents within a reaction. Contains 
   a key-value pair for every reaction where the key is dbse-rxn and the 
   value is an array with the ordered indices of the reagents in that reaction.

.. py:data:: RXNM

   Dictionary defining stoichiometric reagent contributions within a 
   reaction. Contains a key-value pair for every reaction where the key 
   is dbse-rxn and the value is a dictionary whose keys are the contributing 
   reagent indices and whose values are the stoichiometric indices.

.. py:data:: BIND

   Dictionary defining reference values for the chemical quantities defined 
   by the reactions. Contains a key-value pair for every reaction where 
   the key is dbse-rxn and the value is a float. If not natively in kcal/mol, 
   enter in native units, then convert whole dictionary at once.

.. py:data:: TAGL

   Dictionary of comment lines to label input files or calculation headers. 
   Contains a key-value pair for every reagent and reaction where the key 
   is dbse-rgt or dbse-rxn and the value is a string.

.. py:data:: DATA

   Optional dictionary of additional data pertaining to reactions or reagents. 
   In limited use at present for nuclear repulsion energies (handy for 
   checking integrity of GEOS objects) and SAPT data (for coloring plots). 
   Likely to be deprecated.

.. todo:: Because the database files are so plain in format (and because 
   this was the first thing I ever did in python), all the content code 
   gets executed upon import, which is unpythonic and can make imports 
   expensive (perhaps creating thousands of qcdb.Molecule objects).
   The plain-ness needs to be retained, but I'm open to defining GEOS
   as strings (to be fed into qcdb.Molecule constructor, *not* exec(),
   if definition of monomers from dimer could be preserved) or to importing
   instead a wrapped database, where S22 is an instance of class database,
   or to loading once then pickling, shelving, or HTF5-ing, with pieces
   to be served up upon request.

.. todo:: The distinction of what goes in the database file and what is 
   "data" belonging to the database is not rigid. I'm inclined to minimize 
   the former simply because I know how rapidly the latter can accumulate.



