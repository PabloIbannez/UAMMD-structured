##########
GlobalData
##########

global data intro

*****
Units
*****

units intro

None
====

Data entry description:

* **type**: ``Units``, ``None``.
* **parameters**: ``None``.
* **data**: ``None``.

Example:

.. code-block:: json

   "entryName": {
      "type": ["Units", "None"],
   }

KcalMol_A
=========

Data entry description:

* **type**: ``Units``, ``KcalMol_A``.
* **parameters**: ``None``.
* **data**: ``None``.

Example:

.. code-block:: json

   "entryName": {
      "type": ["Units", "KcalMol_A"],
   }

*****
Types
*****

types intro

Basic
=====

Data entry description:

* **type**: ``Types``, ``Basic``.
* **parameters**: ``None``.
* **data**:

   .. list-table::
      :widths: 25 25 25 25
      :header-rows: 1
      :align: center

      * - name
        - mass
        - radius
        - charge
      * - ``string``
        - ``float``
        - ``float``
        - ``float``

Example:

.. code-block:: json

    "entryName":{
      "type":["Types","Basic"]
      "labels":["name","mass","radius","charge"],
      "data":[
         ["A",1.5,1.0, 1.0],
         ["B",1.0,0.5, 0.0],
         ["C",2.0,0.5,-1.0]
      ],
    }


***********
Fundamental
***********

fundamental intro

Time
====

Data entry description:

* **type**: ``Fundamental``, ``Time``.
* **parameters**:

  * ``currentStep`` : ``unsigned long long int``, *optional*, default: 0.

  * ``simulationTime`` : ``float``, *optional*, default: 0.0.

  * ``timeStep`` : ``float``, *optional*.

* **data**: ``None``.

Example:

.. code-block:: json

   "entryName": {
     "type": ["Fundamental", "Time"],
     "parameters": {
       "currentStep": 0,
       "timeStep": 0.001,
       "simulationTime": 0.0
     }
   }

DynamicallyBondedPatchyParticles
================================

Data entry description:

* **type**: ``Fundamental``, ``DynamicallyBondedPatchyParticles``.
* **parameters**:

  * ``energyThreshold`` : ``float``, *optional*, default: 0.0.

* **data**: ``None``.

Example:

.. code-block:: json

   "entryName": {
     "type": ["Fundamental", "DynamicallyBondedPatchyParticles"],
     "parameters": {
       "energyThreshold": -1.0
     }
   }

None
====

Data entry description:

* **type**: ``Fundamental``, ``None``.
* **parameters**: ``None``.
* **data**: ``None``.

Example:

.. code-block:: json

   "entryName": {
     "type": ["Fundamental", "None"]
   }

********
Ensemble
********

NVT
===

Data entry description:

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

Example:

.. code-block:: json

   "entryName": {
     "type": ["Ensemble", "NVT"],
     "labels": ["box", "temperature"],
     "data": [
        [[10.0, 10.0, 10.0], 1.0]
      ]
   }
