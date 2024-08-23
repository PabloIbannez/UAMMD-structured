Bond4
=====

Bond4 potentials in **UAMMD-structured** are interactions involving four particles. These are primarily used to model dihedral angles in molecular simulations. A dihedral angle is defined by four particles (i, j, k, l).

The dihedral angle :math:`\phi` is the angle between the plane containing particles i, j, k and the plane containing particles j, k, l:

.. image:: /img/dihedral.png
    :align: center
    :scale: 50%

The energy and forces are typically functions of this dihedral angle :math: `\phi`, allowing for the modeling of various molecular conformations and rotational barriers.

.. code-block:: yaml

   topology:
     forceField:
       # ...
       bond4_example:
         type: ["Bond4", "Dihedral"]
         labels: ["id_i", "id_j", "id_k", "id_l", "n", "K", "phi0"]
         data:
           - [16, 17, 18, 19, 1, 43.2, 0.1]
       # ...

----

Available dihedral (:math:`\phi` dependent) potentials:

.. toctree::
   :maxdepth: 1

   Dihedral
   IDP_Fourier
