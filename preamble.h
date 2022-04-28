/*Raul P.Pelaez 2021. Extensions preamble.
  This file will be the first uammd related file included if UAMMD_EXTENSIONS is defined at compilation.
  
  You may add preprocessor macros leveraging this assumption.

  Additionally, two special macros are defined here:

    EXTRA_PARTICLE_PROPERTIES: 
       -A list of additional particle properties needed for your extensions. See below for details.
       -Adding a property here will make it available along the rest in ParticleData.

    EXTRA_UPDATABLE_PARAMETERS:
       -A list of additional parameters that will be included in the ParameterUpdatable interface.
       -Refer to the wiki page for ParameterUpdatable for a list of parameters already available in UAMMD.

        

 */
#ifndef UAMMD_EXTENSIONS_PREAMBLE_H
#define UAMMD_EXTENSIONS_PREAMBLE_H

//Apend a list of particle properties you need to this list with the following syntax:
//(([Capitalized name],[lower case name], [type]))
//For example, if ((Custom, custom, real4)) is in this list then the property "Custom"
// will be available along the rest exposed by ParticleData (like pos, vel, force,...)
// meaning that if a ParticleData instance called "pd" is available, you may do:
// auto custom = pd->getCustom(access::gpu, access::read);
#define EXTRA_PARTICLE_PROPERTIES ((ResId  , resId  , int))  \
                                  ((ChainId, chainId, int))  \
                                  ((ModelId, modelId, int))  \
                                  ((SimulationId, simulationId, int))  \
                                  ((FrictionConstant, frictionConstant, real)) \
                                  ((TranslationalSelfDiffusion, TranslationalSelfDiffusion, real)) \
                                  ((RotationalSelfDiffusion, rotationalSelfDiffusion, real)) \
                                  ((SASA, SASA, real)) \
                                  ((InnerRadius, innerRadius, real)) \
                                  ((Epsilon, epsilon, real)) \
                                  ((Surface, surface, real)) \
                                  ((Torque, torque, real4))\
                                  ((Dir, dir, real4)) 

//Append to this list any extra parameter that you need included in the ParameterUpdatable interface
//Use the following syntax:
// (([Capitalized name], [type]))
//For example, if ((MyParameter, real)) is in this list,
// You may safely assume any ParameterUpdatable entity (such as Interactor) exposes a function with the signature
// void updateMyParameter(real newMyParameter);
//You may use this function equivalently to the ones already available in UAMMD (like updateTemperature of updateBox)
#define EXTRA_UPDATABLE_PARAMETERS ((WriteBackup, std::ofstream&)) \
                                   ((ReadBackup, std::ofstream&))


#endif
