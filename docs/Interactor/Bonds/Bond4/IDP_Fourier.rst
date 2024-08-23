IDP_Fourier
-----------

The IDP_Fourier potential is a specialized dihedral potential designed for modeling Intrinsically Disordered Proteins (IDPs) [Smith2015]_ .

.. math::

    U = \sum_{s=1}^4 [A_s \cos(s\phi) + B_s \sin(s\phi)]

where :math:`A_s` and :math:`B_s` are coefficients determined from simulations of IDPs.

----

* **type**: ``Bond4``, ``IDP_Fourier``
* **parameters**: ``None``
* **data**:

  * ``id_i``: ``int``: Id of the first particle
  * ``id_j``: ``int``: Id of the second particle
  * ``id_k``: ``int``: Id of the third particle
  * ``id_l``: ``int``: Id of the fourth particle

Example:

.. code-block::

   "idpFourierBonds":{
     "type":["Bond4","IDP_Fourier"],
     "parameters":{},
     "labels":["id_i", "id_j", "id_k", "id_l"],
     "data":[[0, 1, 2, 3],
             [1, 2, 3, 4],
             [2, 3, 4, 5]]
   }

.. warning::
    The IDP_Fourier potential uses predefined coefficients and does not require additional parameters per bond [Smith2015]_.
    The IDP_Fourier potential is only available in the KcalMol_A unit system.If other units are used, an error will be raised.

.. [Smith2015] Smith, W. Wendell, Po-Yi Ho, and Corey S. O'Hern. "Calibrated Langevin-dynamics simulations of intrinsically disordered proteins." Physical Review E 90, no. 4 (2014): 042709.
