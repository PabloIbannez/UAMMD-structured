State
=====

In the *State* Section, the simulation is initialized by adding particles, each equipped with unique attributes such as initial position and velocity. It is mandatory to assign each particle an **identifier** (Id) and a position. The Id serves as a unique tag for each particle, enabling their identification throughout the simulation for tracking purposes.

In defining the particles present in the simulation, specifying an *Id*, *position*, and possibly other properties is equivalent to creating a particle. In the following example of *State*, when defining a particle, we assign it a position and a velocity:

.. code-block:: yaml

   state:
     labels: ["id", "position", "velocity"]
     data:
       - [0, [1.6, 2.7, 3.3], [0.1, -1.8, 2.3]]
       - [1, [2.5, -3.7, -1.9], [-1.2, 3.7, 1.1]]
       - [2, [-3.3, 1.0, 2.1], [-3.2, 1.3, -2.6]]
       # ...

An important available property is *direction*, which defines the orientation of a particle in space. UAMMD-structured incorporates various potentials, such as *patchy particles*, which require the definition of this property. Direction is specified using quaternions, which can be simplistically understood as a set of four numbers equivalent to a rotation. When a particle is assigned a direction, represented as q=[q0,q1,q2,q3], it corresponds to the rotation matrix:

.. math::

   R = \begin{bmatrix}
   1 - 2q_2^2 - 2q_3^2 & 2q_1q_2 - 2q_0q_3 & 2q_1q_3 + 2q_0q_2 \\
   2q_1q_2 + 2q_0q_3 & 1 - 2q_1^2 - 2q_3^2 & 2q_2q_3 - 2q_0q_1 \\
   2q_1q_3 - 2q_0q_2 & 2q_2q_3 + 2q_0q_1 & 1 - 2q_1^2 - 2q_2^2
   \end{bmatrix}

To utilize this information in encapsulating the orientation of the particle, the particle is assumed to carry a co-moving reference frame. To express the orientation of this frame, the laboratory reference system is used. Thus, the co-moving reference system can be described as a rotation of the laboratory reference system, precisely given by the above equation.

.. figure:: /img/frame.png
   :alt: Quaternion representation for particle orientation

   Particle orientation is encapsulated as a quaternion

In some cases, properties like mass or charge may vary significantly between particles, making it impossible to express them uniformly through types. Hence, these properties can also be defined within the State section. If mass or charge are specified in both *State* and *Types*, the values in *State* take precedence. This flexibility allows for the specification of individual particle characteristics when necessary.
