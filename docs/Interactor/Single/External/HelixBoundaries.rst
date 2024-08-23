HelixBoundaries
---------------

The HelixBoundaries potential confines particles within a helical boundary, 
typically used for simulating systems with helical symmetry or chiral confinement.

This potentential is an extension of `ExternalTabulated <ExternalTabulated.html>`_ potential.
ExternalTabulated potential takes as argument the value of the energy and force at a set of grid points.
In the case of HelixBoundaries these values are computed to represent a helical boundary.

The helix boundary is defined as:

.. math::
   U(r_{\text{h}}) = 
   \begin{cases}
   \begin{align*}
    0 \;\;\;& r_{\text{h}} \leq R_{\text{helixInnerRadius}} \\
    \frac{1}{2} K (r_{\text{h}} - R_{\text{helixInnerRadius}})^2 \;\;\;& r_{\text{h}} > R_{\text{helixInnerRadius}} \\
   \end{align*}
   \end{cases}

where :math:`r_{\text{h}}` is the distance from the helix axis, :math:`R_{\text{helixInnerRadius}}` is the radius of the tube that follows the helix
(the helix itself has a radius of :math:`R_{\text{helix}}`), and :math:`K` is the strength of the confining potential.

We compute the previous potential in every point of a grid that is defined by the parameters ``nx``, ``ny``, and ``nz``, and then we interpolate the values of the potential.

Since the value of :math:`r_{\text{h}}` is not easily computed, we proceed as follows:

1. We generate :math:`n_{\text{pointsHelix}}` points along the helix, with a pitch of :math:`R_{\text{helixPitch}}` and a radius of :math:`R_{\text{helix}}`.
2. For each grid point, we compute the distance to that point to every point in the helix axis, the lowest distance is considered as :math:`r_{\text{h}}`.
3. We compute the value of the potential and its gradient at each grid point.

The second step is done in the GPU. This step is the most computationally expensive, since it requires :math:`O(n_{\text{pointsHelix}} \cdot nx \cdot ny \cdot nz)` operations.

----

* **type**: ``External``, ``HelixBoundaries``
* **parameters**:

  * ``helixPitch``: ``real``: Pitch of the helix :math:`[distance]`
  * ``helixRadius``: ``real``: Radius of the helix :math:`[distance]`
  * ``eps``: ``real``: Handedness of the helix (1.0 for right-handed, -1.0 for left-handed, default: 1.0)
  * ``helixInnerRadius``: ``real``: Inner radius of the confining potential :math:`[distance]`
  * ``K``: ``real``: Strength of the confining potential :math:`[energy/distance^2]`
  * ``nTurns``: ``int``: Number of turns in the helix
  * ``nPointsHelix``: ``int``: Number of discretization points along the helix
  * ``nx``, ``ny``, ``nz``: ``int``: Number of grid points for the tabulated potential

Example:

.. code-block::

   "helixBoundaries":{
     "type":["External","HelixBoundaries"],
     "parameters":{
       "helixPitch": 10.0,
       "helixRadius": 5.0,
       "eps": 1.0,
       "helixInnerRadius": 2.0,
       "K": 100.0,
       "nTurns": 5,
       "nPointsHelix": 1000,
       "nx": 64, "ny": 64, "nz": 64
     }
   }


.. note::
   This potential is particularly useful for simulating systems with helical symmetry, such as certain biological structures or chiral nanotubes.

.. warning::
   The potential calculation can be computationally intensive due to the discretization and tabulation process. Ensure that the number of grid points and helix discretization points are sufficient for accurate representation without being excessively large.

.. tip::
   The ``helixInnerRadius`` and ``K`` parameters can be adjusted to control the "softness" of the confinement. A larger ``K`` will result in a steeper, more rigid boundary.
