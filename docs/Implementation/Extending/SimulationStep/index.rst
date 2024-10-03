Extending SimulationStep
========================

Adding a new simulationStep is necessary if your simulation requires an 
output type that is not implemented in UAMMD-structured. However, it is also useful 
if you need to execute any type of code at specific ``intervalStep``, 
a parameter it receives as input.

A ``SimulationStep`` class inherits from ``SimulationStepBase``, and 
have 3 common methods: 

- **init()** Executes code when the first time the SimulationStep is called. 
- **applyStep()** Executes code every ``intervalStep``. 
- **update()** Executes code every step. It is not mandatory that your SimulationStep 
  have this function defined. **Default do nothing**. 

When programming yout SimulationSteps, take into account that the constructor 
is called with ``ParticleGroup``, ``IntegratorManager``, ``ForceField`` and 
``DataEntry``. So you can access all the information of the simulation in the 
SimulationStep methods. 

.. warning::
   Take into account that you are accessing to pointers of the simulation information,
   so you can actually modify this information. This could be very useful, but dangerous. 
   Make sure you know what your are doing when modify this information or your simulation 
   might give incorrent results. 

One option to manage this potential problem is to create a memory buffer. Make 
a copy of the simulation data at the begining of the simulationStep so you don't 
modify it content.

.. code-block:: cpp

    #include "System/ExtendedSystem.cuh"
    #include "GlobalData/GlobalData.cuh"
    #include "ParticleData/ExtendedParticleData.cuh"
    #include "ParticleData/ParticleGroup.cuh"

    #include "SimulationStep/SimulationStep.cuh"
    #include "SimulationStep/SimulationStepFactory.cuh"

    namespace uammd{
    namespace structured{
    namespace SimulationStep{
    namespace Family{

    class myStep: public SimulationStepBase{

        public:

            myStep(std::shared_ptr<ParticleGroup>  pg,
                      std::shared_ptr<IntegratorManager> integrator,
                      std::shared_ptr<ForceField>    ff,
                      DataEntry& data,
                      std::string name):SimulationStepBase(pg,integrator,ff,data,name){

                      // Initialize your class
            }

            void init(cudaStream_t st) override{

                // Executed once, the first time you call the SimulationStep
            }

            void applyStep(ullint step, cudaStream_t st) override{

                // Executed every intervalStep
            }

            void update(cudaStream_t st) override{

                // Executed every step
            }
    };

    }}}}

    REGISTER_SIMULATION_STEP(
            myType,myStep,
            uammd::structured::SimulationStep::Family::myStep
    )

To register your own SimulationStep create the file
``src/SimulationStep/Family/myType/myStep.cu`` and add to
the ``Components.json``.

.. code-block:: json
   :emphasize-lines: 5

   {
   "SimulationStep":
        "Family":[
                ["..."],
                ["myType","myStep","myStep.cu"]
                ]
   }

