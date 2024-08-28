DihedralsMeasure
----------------

The DihedralsMeasure step calculates and records the dihedral angles for specified sets of four particles over time.

Output format:

.. code-block::

   # Step 1
   dihedral1 dihedral2 dihedral3 ...
   # Step 2
   dihedral1 dihedral2 dihedral3 ...
   # Step 3
   dihedral1 dihedral2 dihedral3 ...
   ...

----

* **type**: ``ParticlesListMeasure``, ``DihedralsMeasure``
* **parameters**:

  * ``outputFilePath``: ``string``: Path to the output file

* **data**:

  * ``id_i``: ``int``: ID of the first particle in the dihedral
  * ``id_j``: ``int``: ID of the second particle in the dihedral
  * ``id_k``: ``int``: ID of the third particle in the dihedral
  * ``id_l``: ``int``: ID of the fourth particle in the dihedral

Example:

.. code-block::

   "dihedralsMeasure":{
     "type":["ParticlesListMeasure","DihedralsMeasure"],
     "parameters":{
       "intervalStep": 10000,
       "outputFilePath": "dihedrals.dat"
     },
     "labels":["id_i", "id_j", "id_k", "id_l"],
     "data":[
       [1, 2, 3, 4],
       [2, 3, 4, 5],
       [3, 4, 5, 6]
     ]
   }

.. note::
   The output file will contain the step number and the dihedral angle for each specified set of four particles.

.. tip::
   This measure is particularly useful for analyzing the conformational changes of molecules, such as protein backbone dihedrals or torsion angles in polymers.
