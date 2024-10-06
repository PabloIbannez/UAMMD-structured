#pragma once

#define COMPUTABLES \
((energy)(Energy))\
((force)(Force))\
((torque)(Torque))\
((virial)(Virial))\
((stress)(Stress))\
((magneticField)(MagneticField))\
((transitionProbability)(TransitionProbability))\
((hessian)(Hessian))\
((lambdaDerivative)(LambdaDerivative))\
((pairwiseForce)(PairwiseForce))\


#define COMPUTABLES_COMBINATIONS \
((energyForceTorque)(EnergyForceTorque)(energy)(force)(torque))\
((energyForce)(EnergyForce)(energy)(force))\
((forceTorqueMagneticField)(ForceTorqueMagneticField)(force)(torque)(magneticField))\
((forceTorque)(ForceTorque)(force)(torque))\

