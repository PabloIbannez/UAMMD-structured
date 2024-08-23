DihedralsMeasure
----------------

The DihedralsMeasure step calculates and records the dihedral angles for specified sets of four particles over time.

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
