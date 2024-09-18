Scheduling the components
=========================

Both interactions and simulationSteps in **UAMMD-structured** share a set of common parameters that allow users to specify when these components should be activated or deactivated during the simulation. These parameters provide fine-grained control over the simulation process. The two key parameters are `startStep` and `endStep`.

The `startStep` parameter determines when a component (interaction or simulationStep) becomes active during the simulation. Its default value is 0, meaning the components are active from the beginning of the simulation.

Conversely, the `endStep` parameter specifies when a component should be deactivated during the simulation. Its default value is INFINITY, ensuring that the component remains active throughout the entire simulation unless explicitly set. When the simulation reaches the specified `endStep`, the component is deactivated.

These parameters can be used individually or in combination to achieve precise control over when components are active during the simulation. This flexibility allows for complex simulation scenarios where different interactions or steps are enabled or disabled at specific points in time.

For example, you might configure a Lennard-Jones interaction to become active at step 1000 and deactivate at step 5000 like this:

.. code-block::

   "lj":{
     "type":["NonBonded","LennardJonesType1"],
     "parameters":{
       "cutOffFactor":2.5,
       "condition":"all",
       "startStep":1000,
       "endStep":5000
     },
     "labels":["name_i", "name_j", "epsilon", "sigma"],
     "data":[
       ["A", "A", 1.0, 1.0],
       ["A", "B", 1.2, 0.9],
       ["B", "B", 0.8, 1.1]
     ]
   }

In this case, the Lennard-Jones interaction would only be active for steps 1000 through 4999 of the simulation. This feature can be particularly useful for modeling systems with time-dependent interactions or for gradually introducing or removing forces during a simulation.
