Extending Integrator
====================

UAMMD-structured uses UAMMD Integrators as start point.
So it is very recomendable to read
`UAMMD <https://uammd.readthedocs.io/en/latest/Integrator/index.html>`_
documentation before this section.

Every UAMMD-structured Integrator inherits from ``IntegratorBase``
(or ``IntegratorBaseNVT`` that add some common features for
NVT Integrators, such as ``KbT``). ``void forwardTime()`` is the only mandatory function
of an Integrator and is responsible for advancing the simulation's time.
You can think of forwardTime() as a single iteration of a loop within the simulation.
Usually, it computes all forces, energies, velocities, and updates the positions
of objects in the simulation. Interactors can be invoked from this function.

The simplest example of an Integrator is:

.. code-block:: cpp

    #include "System/ExtendedSystem.cuh"
    #include "GlobalData/GlobalData.cuh"
    #include "ParticleData/ExtendedParticleData.cuh"
    #include "ParticleData/ParticleGroup.cuh"

    #include "Integrator/IntegratorBase.cuh"
    #include "Integrator/IntegratorFactory.cuh"

    namespace uammd{
    namespace structured{
    namespace Integrator{
    namespace myType{
    namespace myIntegrator{

        class  myIntegrator: public IntegratorBase{

            public:

                myIntegrator(std::shared_ptr<GlobalData> gd,
                             std::shared_ptr<ParticleGroup>   pg,
                             DataEntry& data,
                             std::string name): IntegratorBase(gd,pg,data,name){
                }

                void forwardTime() override {

                    // Do a step of the simulation

                    // Advance CurrentStep and SimulationTime is not mandatory, but very recomendable.
                    this->gd->getFundamental()->setCurrentStep(this->gd->getFundamental()->getCurrentStep()+1);
                    this->gd->getFundamental()->setSimulationTime(this->gd->getFundamental()->getSimulationTime()+this->dt);
                }

        };

    }}}}}

    REGISTER_INTEGRATOR(
            myType,myIntegrator,
            uammd::structured::Integrator::myType::myIntegrator
            )

To register your own Integrator create the file
``src/Integrator/Family/myType/myIntegrator.cu`` and add to
the ``Components.json``.

.. code-block:: json
   :emphasize-lines: 5

   {
   "Integractor":
        "Family":[
                 ["..."],
                 ["myType","myIntegrator","myIntegrator.cu"]
                 ]
   }

