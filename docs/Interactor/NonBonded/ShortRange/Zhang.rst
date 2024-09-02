Zhang
-----

The Zhang potential is a coarse-grained, one-particle-thick model for biological and biomimetic fluid membranes developed by Yuan et al. [Yuan2010]_ . It features an anisotropic pair potential that allows particles to self-assemble into fluid membranes with biologically relevant properties.

The potential consists of two parts:

1. A distance-dependent function :math:`U(r)`
2. An orientation-dependent function :math:`\Omega(r̂_{ij}, n_{i}, n_{j})`

The total potential is given by:

.. math::
    U(\mathbf{r}_{ij}, \mathbf{n}_i, \mathbf{n}_j) =
    \begin{cases}
        U_R(r) + \varepsilon[1 - \Omega(\hat{\mathbf{r}}_{ij}, \mathbf{n}_i, \mathbf{n}_j)], & r \leq r_{\text{min}} \\
        U_A(r)\Omega(\hat{\mathbf{r}}_{ij}, \mathbf{n}_i, \mathbf{n}_j), & r_{\text{min}} < r \leq r_c
    \end{cases}

The distance-dependent functions:

.. math::
    U_R(r) = \varepsilon \left[\left(\frac{r_{\text{min}}}{r}\right)^4 - 2\left(\frac{r_{\text{min}}}{r}\right)^2\right]

.. math::
    U_A(r) = -\varepsilon \cos^2\left[\frac{\pi}{2}\frac{r - r_{\text{min}}}{r_c - r_{\text{min}}}\right]

The orientation-dependent function:

.. math::
    \Omega = 1 + \mu[(\mathbf{n}_i \cdot \hat{\mathbf{r}}_{ij})(\mathbf{n}_j \cdot \hat{\mathbf{r}}_{ij}) + \sin\theta_0(\mathbf{n}_j - \mathbf{n}_i) \cdot \hat{\mathbf{r}}_{ij} - \sin^2\theta_0]

where :math:`r_{ij}` is the distance between particles, :math:`n_i` and :math:`n_j` are unit vectors representing particle orientations, :math:`r_{\text{min}}` is the distance of minimum energy, and :math:`r_c` is the cutoff distance.

----

* **type**: ``NonBonded``, ``Zhang``
* **parameters**:

  * ``alpha``: ``real``: Controls the slope of the attractive potential branch
  * ``beta``: ``real``: Related to membrane bending rigidity
  * ``gamma``: ``real``: Related to membrane spontaneous curvature
  * ``r_min``: ``real``: Distance of minimum energy [distance]
  * ``r_c``: ``real``: Cutoff distance [distance]
  * ``axis``: ``real3``: Reference axis for particle orientations (default: {0, 0, 1})

* **data**:

  * ``name_i``: ``string``: Type of particle i
  * ``name_j``: ``string``: Type of particle j
  * ``epsilon``: ``real``: Energy scale [energy]

Example:

.. code-block::

   "zhang":{
     "type":["NonBonded","Zhang"],
     "parameters":{
       "alpha": 4.0,
       "beta": 3.0,
       "gamma": 0.0,
       "r_min": 1.1225,
       "r_c": 2.6,
       "axis":[0, 0, 1]
     },
     "labels":["name_i", "name_j", "epsilon"],
     "data":[
       ["A", "A", 1.0],
       ["A", "B", 1.0],
       ["B", "B", 1.0]
     ]
   }

.. note::
   This potential requires orientation information for each particle. Ensure that the orientation (director) is properly set for each particle before using this potential.

.. [Yuan2010] Yuan, H., Huang, C., Li, J., Lykotrafitis, G., & Zhang, S. (2010). One-particle-thick, solvent-free, coarse-grained model for biological and biomimetic fluid membranes. Physical Review E, 82(1), 011905.
