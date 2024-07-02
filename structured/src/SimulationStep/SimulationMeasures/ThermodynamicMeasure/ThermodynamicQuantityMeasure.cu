#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "SimulationStep/SimulationStep.cuh"
#include "SimulationStep/SimulationStepFactory.cuh"

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{

class ThermodynamicQuantityMeasure: public SimulationStepBase_EnergyForceTorque{

        std::string    outputFilePath;
        std::ofstream  outputFile;

    public:

        ThermodynamicQuantityMeasure(std::shared_ptr<ParticleGroup>             pg,
                                     std::shared_ptr<IntegratorManager> integrator,
                                     std::shared_ptr<ForceField>                ff,
                                     DataEntry& data,
                                     std::string name):SimulationStepBase_EnergyForceTorque(pg,
                                                                                            integrator,ff,
                                                                                            data,name){

            //Read parameters

            outputFilePath = data.getParameter<std::string>("outputFilePath");
        }


        void init(cudaStream_t st) override{

            bool isFileEmpty = Backup::openFile(this->sys,outputFilePath,outputFile);

            if(isFileEmpty){
                outputFile << "#Step ";
                outputFile << "ParticleNumber Volume ";

                for(auto& interactor : this->topology->getInteractors()){
                    outputFile << "Energy(" << interactor.first << ") ";
                }

                outputFile << "KineticEnergy TotalPotentialEnergy TotalEnergy ";
                outputFile << "Temperature ";
                outputFile << "Virial " << std::endl;
            }
        }

        void applyStep(ullint step, cudaStream_t st) override{

            //Reset energy, force
            this->setZero(st);

            //Write current step to output file
            outputFile << step << " ";

            //Write particle number to output file
            outputFile << this->pg->getNumberParticles() << " ";

            //Write volume to output file
            Box box = gd->getEnsemble()->getBox();
            real V = box.boxSize.x*box.boxSize.y*box.boxSize.z;
            outputFile << V << " ";

            real totalPotentialEnergy = 0.0;
            //Write energy to output file
            {
                //Iterate over all interactors
                for(auto& interactor : this->topology->getInteractors()){
                    //Fill energy with zeros
                    {
                        auto energy = this->pd->getEnergy(access::location::gpu, access::mode::write);
                        thrust::fill(thrust::cuda::par.on(st),
                                     energy.begin(),
                                     energy.end(),
                                     real(0));
                    }

                    //Create computable
                    uammd::Interactor::Computables compTmp;
                    compTmp.energy = true;

                    interactor.second->sum(compTmp,st);
                    real totalEnergy = Measures::totalPotentialEnergy(pg);
                    totalPotentialEnergy += totalEnergy;

                    outputFile << totalEnergy << " ";

                }
            }

            //Write kinetic energy to output file
            real kineticEnergy = Measures::totalKineticEnergy(pg);
            outputFile << kineticEnergy << " ";

            //Write total potential energy to output file
            outputFile << totalPotentialEnergy << " ";

            //Write total energy to output file
            outputFile << kineticEnergy + totalPotentialEnergy << " ";

            //Temperature
            real T = real(2.0/(3.0*this->pg->getNumberParticles()*gd->getUnits()->getBoltzmannConstant()))*kineticEnergy;
            outputFile << T << " ";

            //Compute virial
            {
                //Iterate over interactors and sum forces
                uammd::Interactor::Computables compTmp;
                compTmp.force = true;
                for(auto& interactor : this->topology->getInteractors()){
                    interactor.second->sum(compTmp,st);
                }

                {
                    Box box = gd->getEnsemble()->getBox();

                    auto force    = this->pd->getForce(access::location::cpu, access::mode::read);
                    auto position = this->pd->getPos(access::location::cpu, access::mode::read);

                    auto groupIndex = this->pg->getIndexIterator(access::location::cpu);

                    real virial = 0.0;
                    for(int i = 0; i < this->pg->getNumberParticles(); ++i){
                        int index = groupIndex[i];
                        virial += dot(make_real3(force[index]),box.apply_pbc(make_real3(position[index])));
                    }
                    virial *= (-0.5/pg->getNumberParticles());

                    outputFile << virial << " ";
                }
            }

            outputFile << std::endl;
        }
};

}}}}

REGISTER_SIMULATION_STEP(
    ThermodynamicMeasure,ThermodynamicQuantityMeasure,
    uammd::structured::SimulationStep::SimulationMeasures::ThermodynamicQuantityMeasure
)
