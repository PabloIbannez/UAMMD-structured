NVTlambda
----------

This ensemble extends the `NVT ensemble <NVT.html>`_ by adding a lambda parameter. This
parameter is used to implement alchemical free energy calculations. This parameter is used
to turn on and off some interactions in the system. 

* **type**: ``Ensemble``, ``NVTlambda``
* **parameters**: ``None``.
* **data**:

  .. list-table::
     :widths: 25 25 25
     :header-rows: 1
     :align: center

     * - temperature
       - box
       - :math:`\lambda`
     * - ``float``
       - [``float``, ``float``, ``float``]
       - ``float``

----

Example:

.. code-block:: json

   "entryName": {
     "type": ["Ensemble", "NVT"],
     "labels": ["box", "temperature","lambda"],
     "data": [
        [[10.0, 10.0, 10.0], 1.0, 1.0]
      ]
   }

.. note::
   Not all interactions are compatible :math:`\lambda` mechanism. Those compatible are explicitly indicated in the documentation.

.. warning::
   The :math:`\lambda` parameter is a float between 0 and 1. Otherwise, an error will be raised.

