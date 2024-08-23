NVT
----

Selecting this ensembale will define a global temperature and volume.
These quantities are used for the other components of the simulation.
For example, integrators can use the temperature to fix the dynamic properties of
the particles in order to simulate a system at a given temperature
(which commonly means reproducing the canonical Boltzmann distribution).

The volume is given as a 3D vector, which is used to define the size of the simulation box.

* **type**: ``Ensemble``, ``NVT``.
* **parameters**: ``None``.
* **data**:

  .. list-table::
     :widths: 25 25
     :header-rows: 1
     :align: center

     * - temperature
       - box
     * - ``float``
       - [``float``, ``float``, ``float``]

----

Example:

.. code-block:: json

   "entryName": {
     "type": ["Ensemble", "NVT"],
     "labels": ["box", "temperature"],
     "data": [
        [[10.0, 10.0, 10.0], 1.0]
      ]
   }
