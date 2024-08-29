Groups
======

In certain situations, we want a specific interaction or a particular simulation step to apply only to a subset of particles. For such scenarios, groups are defined. We define a group by specifying its name, the type of selection to be made, and the selection itself:

.. code-block:: yaml

   topology:
     forceField:
       groups_example:
         type: ["Groups", "GroupsList"]
         labels: ["name", "type", "selection"]
         data:
           - ["ids_group_1", "Ids", [0, 1, 2, 3, 4]]
           - ["ids_group_2", "Ids", [5, 6, 7, 8, 9]]
       interaction1:
         parameters:
           group: "ids_group_1"
       interaction2:
         parameters:
           group: "ids_group_2"

In this example, we define two groups, "ids_groups_1" and "ids_groups_2". The definition of these groups is made through IDs, followed by indicating the selection. Subsequently, we specify in each interaction the group it acts upon using the `group` parameter. Similarly, we can define selections for simulation steps:

.. code-block:: yaml

   simulationStep:
     groups_example:
       type: ["Groups", "GroupsList"]
       labels: ["name", "type", "selection"]
       data:
         - ["types_group_1", "Types", ["A"]]
         - ["types_group_2", "notTypes", ["A"]]
     simStep1:
       parameters:
         group: "types_group_1"
     simStep2:
       parameters:
         group: "types_group_2"

In this example, group definitions are made using "Types". In the first selection "types_group_1", we specify that it comprises particles of type "A". The second defined group is "types_group_2", which consists of particles that are not of type "A".

Finally, we can also declare groups similarly in integrators. This can be useful for constructing complex boundaries that cannot be expressed well through a potential. We can create a system formed by two groups, "particles" and "boundaries", where only the first is integrated, but potentials apply to both (no group defined in the force field). Thus, only "particles" will advance in time but will feel the interactions of the "boundaries" particles.
