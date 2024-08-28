GeometricalMeasure
==================

The GeometricalMeasure simulation steps are designed to analyze and record various geometric properties of the simulated system over time. These measures provide valuable insights into the structural and dynamical behavior of particles or molecular systems.

These simulation steps typically calculate and output specific geometric quantities at regular intervals during the simulation. They can be used to track changes in system size, particle distributions, molecular conformations, and more.

Some key features of GeometricalMeasure steps include:

- Calculation of system-wide properties like radius of gyration or mean radius
- Tracking of particle displacements and rotations
- Analysis of specific particle groups or subsets
- Recording of time-dependent geometric changes

The following GeometricalMeasure steps are available:

.. toctree::
   :maxdepth: 1

   MeanRadius
   GyrationRadius
   MeanSquareDisplacement
   MeanAngularCorrelation
   CenterOfMassPosition
   EscapeTime
   Height
   DistanceBetweenCentersOfMass
