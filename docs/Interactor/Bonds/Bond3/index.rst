Bond3
=====

Bond3 potentials in **UAMMD-structured** are interactions involving three particles. 
While these can be used for various three-body interactions, 
they are primarily employed to define angular bonds in molecular simulations. 
In the context of angular bonds, the three particles (i, j, k) form an angle, with j being the central particle.

The angle :math:`\theta` is defined as the angle formed between the vectors j-i and j-k:

.. image:: /img/bond_angle.png
    :align: center
    :scale: 15%

The energy and forces are typically functions of this angle θ, :math:`r_{ji} = r_i - r_j` and :math:`r_{jk} = r_j - r_k`, 
where :math:`cos \theta = (\vec{r}_{ji} \cdot \vec{r}_{jk}) / (|r_{ji}| |r_{jk}|)`. 

Hence, we can define potentials such as :math:`U(\theta) = 1/2 K(\theta - \theta_0)^2`, 
which we refer to as "HarmonicAngular" in the following example:

.. code-block:: yaml

   topology:
     forceField:
       # ...
       bond3_example:
         type: ["Bond3", "HarmonicAngular"]
         labels: ["id_i", "id_j", "id_k", "K", "theta0"]
         data:
           - [12, 13, 14, 1.32, 1.15]
           - [13, 14, 15, 3.13, 0.53]
       # ...

----

Available angular (:math:`\theta` dependent) potentials:

.. toctree::
   :maxdepth: 1

   KratkyPorod
   HarmonicAngular
   BestChenHummerAngular
