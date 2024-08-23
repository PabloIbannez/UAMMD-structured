Basic
------

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
