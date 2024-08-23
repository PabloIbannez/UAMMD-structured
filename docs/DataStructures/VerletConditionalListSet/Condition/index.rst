Condition
=========

Using the conditions listed below, we can create various **Verlet lists**. 
Depending on the selected condition, one, two, or more lists will be generated. 
Each created list ensures that all pairs of particles within it satisfy a specific condition.

For example, if we select the *intramolecular-intermolecular* condition, two Verlet lists are created. 
The first list ensures all pairs consist of particles with the same ``modelId`` value, 
while the second list ensures the ``modelId`` value differs for each particle in the pair. 
These lists are associated with labels: "intra" and "inter" respectively.

To employ one of the created lists in a potential, select the corresponding label in the potential's ``condition`` parameter. 
For instance, to use the "inter" condition:

.. code-block:: json

   {
     "condition": "inter"
   }

To implement conditional Verlet lists, first create the conditional list, 
then select the desired condition (from those available) for use in the potential.

All Verlet lists created by a condition share the same Verlet radius. 
Consequently, all lists have ``cutOff`` and ``cutOffVerletFactor`` as parameters. 
The Verlet radius is determined by:

.. math::

   r_{\text{verlet}} = \text{cutOff} \times \text{cutOffVerletFactor}

----

The available conditions are:

.. toctree::
   :maxdepth: 1

   all
   intra_inter
   nonExcluded
   nonExclIntra_nonExclInter
   nonExclIntra_nonExclInter_nonExclCharged
   nonExclIntra_nonExclInter_nonExclInterCharged
   interDifferentType

----

.. tip::

    The cutOff parameter is typically not specified explicitly, as it is externally determined by the potentials present in the simulation.
