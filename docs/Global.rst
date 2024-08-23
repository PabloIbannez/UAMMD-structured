Global
======

The second section to be processed is *Global*. As its name suggests, the *Global* section sets certain overall aspects of the simulation. Unlike the system section, the properties here have a physical significance. Within the global section, we can add four types of entries: **fundamental**, **units**, **types**, and **ensemble**. It only makes sense to add one entry of each type (for example, it does not make sense to indicate two different units systems); doing otherwise will result in an error and prevent the simulation from starting. The following code shows a typical configuration example for the *Global* section:

.. code-block:: yaml

   global:
     fundamental:
       type: ["Fundamental", "Time"]
     ensemble:
       type: ["Ensemble", "NVT"]
       labels: ["box", "temperature"]
       data:
         - [[100.0, 100.0, 100.0], 300.0]
     types:
       type: ["Types", "Basic"]
       labels: ["name", "mass", "radius", "charge"]
       data:
         - ["A", 1.0, 0.5, -1.0]
         - ["B", 1.0, 0.5, 0.0]
         - ["C", 1.0, 0.5, 1.0]
     units:
       type: ["Units", "KcalMol_A"]

.. _global_fundamental:

Fundamental
-----------

This Data Entry is somewhat abstract. In fundamental, we can set the *primary property governing the simulation*. A common example would be to have "Time" as the fundamental type, meaning time is the main controlling quantity, as is typical in molecular dynamics simulations. There are other options for specific cases, but generally, "Time" is used, and if this entry is not specified, "Time" is considered the default.

.. _global_units:

Units
-----

The UAMMD structure is equipped with a unit system featuring various options, including the **Kcal/mol Å** system, commonly known as AKMA in the literature. When some unit system is selected, all input data is interpreted according to its standards. For example in the case of AKMA (widely used in molecular dynamics simulations of proteins, DNA, RNA ...), defines specific units for physical quantities like length, mass, temperature, energy, and others. The units within the AKMA system are as follows:

========================  ====================================  ====================
Quantity                  Unit                                  Symbol/Formula
========================  ====================================  ====================
Length                    Angstrom                              Å
Mass                      Dalton (Atomic Mass Unit)             Da (amu)
Temperature               Kelvin                                K
Charge                    Electron Charge                       e
Energy                    Kilocalorie per mole                  kcal/mol
Time                      AKMA Time Unit                        AKMA Time Unit
Force                     Kilocalorie per mole per Angstrom     kcal/(mol·Å)
Pressure                  Atmospheres                           atm
Velocity                  Angstroms per AKMA Time Unit          Å/(AKMA Time unit)
========================  ====================================  ====================

The AKMA Time Unit is a derived unit in this system, calculated from the primary units of mass, length, and energy. It is equivalent to 4.888821 × 10^-14 seconds. For practical purposes, however, time in this units system is usually represented in picoseconds, with 20 AKMA time units equating to approximately 0.978 picoseconds. But the input and the output are unchanged and when this unit system is used then input and output are assumed to be in AKMA time units.

"None" units system can also be selected, in this case the interpretation of units is left entirely to the user. Internally, when we set a units system, what happens is that the value of certain physical constants are fixed, such as the Boltzmann constant or the vacuum permittivity, when "None" units system is select all physical constants are considered to be equal to one. Currently, only these two units systems ("KcalMol_A" and "None") are available, but new units systems can be added.

.. _global_types:

Types
-----

In certain cases, particles can be **grouped** based on their mass, charge, etc. In this case, it is said that they belong to a certain type of particles. Therefore, instead of assigning certain properties to each particle, we can associate it with a type and have its properties dictated by it. The "Types" entry allows for the selection of particle types for use in the simulation. Choosing a specific type enables certain characteristics to be assigned to each particle. UAMMD-structured currently incorporates only one class of particle type, known as "Basic". When this type is selected, one can define particle types and specify their mass, radius, and charge. The example demonstrates how three particle types, "A", "B", and "C", are defined, all having the same mass and radius but differing in charge. It is important to note that these values are interpreted according to the units system, with mass in Dalton, radius in Angstrom, and charge in multiples of the electron charge, for this example. Importantly, this section defines available types, but the association between types and particles is performed in other input Sections, particularly the "structure" Section within *Topology*.

.. _global_ensemble:

Ensemble
--------

*Many of the physical processes of interest do not occur in a vacuum but rather under specific conditions*, such as a predetermined volume or temperature. These types of global variables are commonly referred to as ensembles.

The *Ensemble* Data Entry in UAMMD-structured is where the ensemble for a simulation is established. Defining the ensemble involves setting certain macroscopic variables of the system. For example, the NVT ensemble choice allows for specifying temperature and volume thereby defining the simulation box, the unit cell where PBC are applied.

However, sometimes the options available under the ensemble category are stretched to include global variables that typically would not be categorized as ensembles in scientific literature. An example within UAMMD-structured is the NVTλ ensemble. In addition to defining the volume and temperature, this ensemble specifies the parameter λ. This parameter is utilized in techniques like **Thermodynamic Integration** (TI), demonstrating a broader application of the ensemble concept in UAMMD-structured simulations.
