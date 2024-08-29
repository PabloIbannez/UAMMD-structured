Topology
========

In the *Global* Section, we establish the general conditions of the simulation, setting aspects like the units, particle types, and the ensemble. In the *State* Section, we introduce particles and define some of their properties, such as position and orientation. The *Integrators* Section is then used to determine the evolution of these particles over time. At this stage, the information input is quite generic, and without further detail, all systems would appear similar. The key to **differentiating and detailing a unique physical system** lies in the *Topology* Section.

The *Topology* Section is where we define the properties of particles that stay constant throughout the simulation, including their type or position within the structure, and the interactions that take place in the system. *Topology* Section consists of a Data Entry, called *structure*, that specifies types and the structural information (this is where its name comes from), and a Section named *force field* for detailing the interactions:

.. code-block:: yaml

   topology:
     structure:
       # ...
     forceField:
       # ...

Structure
---------

In many simulations, particles are characterized by specific features, rather than being anonymous entities. As explained in `Types <Global.html#types>`_, each particle is part of a broader category defined by its type. Additionally, in numerous systems, **each particle holds a specific position within the structure to which it belongs**. Certain models of proteins require potentials for the interaction of two beads belonging to different proteins (inter) different from those used for the interaction between beads of the same protein (intra). In this type of system particles need to be classified according to the structure or substructure they belong to. This categorization process is managed in the *structure* Data Entry. Here, every particle must be assigned a type, which renders various properties such as mass, charge, and radius, dependent on the selected type style. Assigning a type is compulsory, however, one may provide additional structural details allowing for several degrees of structural information. The structure within UAMMD-structured is illustrated in the following example:

.. code-block:: yaml

   topology:
     structure:
       labels: ["id", "type", "modelId"]
       data:
         - [0, "A", 0]
         - [1, "A", 0]
         - [2, "A", 1]
         - [3, "A", 1]
         # ...

UAMMD-structured uses a hierarchical classification system inherited from protein notation. In this system, particles are organized into residues ("**resId**"), chains ("**chainId**") and models ("**modelId**"). These additional data do not modify the particles themselves, but can be used to define particular interactions at different structural levels, e.g., two particles may interact with one type of potential if they belong to the same model, but with a different one if they are located in distinct models. In this example, each particle is assigned not just a type, such as "A", but is also linked to a specific model. For instance, particles 0 and 1 are part of model 0, while particles 2 and 3 belong to model 1. This configuration enables the calculation of unique potentials between particles from different models. Assigning both type and structural information requires including the particle identifier, "id".

The ability of UAMMD-structured to define interactions on various hierarchical levels, based on the type and structural position of particles, makes it a highly effective tool for simulating coarse-grained models.

Force Field
-----------

To fully customize our system, it is essential to specify the types of interactions that take place. The *force field* Section is where the interactions present in the simulation are defined. These interactions vary based on the number of particles involved, their nature, among other factors. In this Section, different Data Entries are incorporated, each representing an interaction. We can also introduce certain Data Entries that activate not interactions per se, rather specific data structures used by the interactions:

.. code-block:: yaml

   topology:
     forceField:
       interaction1:
         # ...
       interaction2:
         # ...
       interaction3:
         # ...
       # ...

The wide variety of available interactions, along with the numerous different integration schemes, are the two cornerstones of UAMMD/UAMMD-structured. When we refer to "interaction," we are discussing a broad concept: an interaction is a specific way to group particles before evaluating a potential on them. For instance, in a bonding interaction, one might assess a harmonic potential or a Morse potential. For short-range non-bonded interactions, a potential like Lennard-Jones might be used. In the case of external interactions, a potential from an electric field or a surface could be applied. We will examine them in detail in the following sections, offering examples of the available potentials when describing an interaction.
