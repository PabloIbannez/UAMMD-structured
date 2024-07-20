#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "SimulationStep/SimulationStep.cuh"
#include "SimulationStep/SimulationStepFactory.cuh"

#include "Utils/Measures/MeasuresBasic.cuh"

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationMeasures{

class ThermodynamicIntegration: public SimulationStepBase{

        std::string    outputFilePath;
        std::ofstream  outputFile;

        ullint stepLambda;

        std::vector<real> lambdaValues;

    public:

        ThermodynamicIntegration(std::shared_ptr<ParticleGroup>             pg,
                                 std::shared_ptr<IntegratorManager> integrator,
                                 std::shared_ptr<ForceField>                ff,
                                 DataEntry& data,
                                 std::string name):SimulationStepBase(pg,
                                                                      integrator,ff,
                                                                      data,name){

            //Read parameters
            outputFilePath    = data.getParameter<std::string>("outputFilePath");

            //Lambda options
            stepLambda           = data.getParameter<ullint>("stepLambda");
            lambdaValues         = data.getParameter<std::vector<real>>("lambdaValues");

            //Write parameters

            System::log<System::MESSAGE>("[ThermodynamicIntegration] outputFilePath: %s",outputFilePath.c_str());
            System::log<System::MESSAGE>("[ThermodynamicIntegration] stepLambda: %i",stepLambda);

            std::string lambdaValuesString;
            for(auto& lambda:lambdaValues){
                lambdaValuesString += std::to_string(lambda) + " ";
            }
            System::log<System::MESSAGE>("[ThermodynamicIntegration] lambdaValues: %s",lambdaValuesString.c_str());
        }


        void init(cudaStream_t st) override {

            System::log<System::MESSAGE>("[ThermodynamicIntegration] Initializing...");

            //Check if the number of batches is larger than one
            GroupUtils::BatchGroupNumberCheck(this->pg,1);

            //Check initial lambda match the current system lambda
            real initLambda = this->lambdaValues[0];
            if(initLambda != this->gd->getEnsemble()->getLambda()){
                System::log<System::CRITICAL>("[ThermodynamicIntegration] Initial lambda (%f) does not match the current system lambda (%f).",
                                              initLambda,this->gd->getEnsemble()->getLambda());
            }

            // Write header

            bool isFileEmpty = Backup::openFile(this->sys,outputFilePath,outputFile);

        }

        void update(ullint step,cudaStream_t st) override { //This function is called every step

            ullint effectiveStep = step - startStep;

            if((effectiveStep % stepLambda == 0) and (effectiveStep/stepLambda < lambdaValues.size())){

                real lambda = lambdaValues[effectiveStep/stepLambda];
                System::log<System::MESSAGE>("[ThermodynamicIntegration] Lambda = %f",lambda);

                // Update lambda
                this->gd->getEnsemble()->setLambda(lambda);
                cudaDeviceSynchronize();
                // Write data
                outputFile << "# " << lambda << std::endl;
            }
        }

        void applyStep(ullint step, cudaStream_t st) override{

            if(step/stepLambda >= lambdaValues.size()){
                if(outputFile.is_open()){
                    System::log<System::MESSAGE>("[ThermodynamicIntegration] Finished!");
                    outputFile.close();
                }
            }

            //Reset lambda derivative
            {
                auto lambdaDerivative = this->pd->getLambdaDerivative(access::location::gpu,
                                                                      access::mode::write);

                thrust::fill(thrust::cuda::par.on(st),
                             lambdaDerivative.begin(),
                             lambdaDerivative.end(),
                             real(0.0));

                cudaStreamSynchronize(st);
            }

            //Compute lambdaDerivative
            {
				        uammd::Interactor::Computables comp;
				        comp.lambdaDerivative = true;

                //Sum energy
    				    for(auto& interactor: this->topology->getInteractors()){
                    interactor.second->sum(comp,st);
                }
            }

            ////DEBUG
            //{
            //    auto ld = this->pd->getLambdaDerivative(access::location::cpu,
            //                                            access::mode::read);

            //    std::ofstream ldFile("ldDebug.dat");

            //    for(int i = 0; i < ld.size(); ++i){
            //        ldFile << i << " " << ld[i] << std::endl;
            //    }

            //}

            {
                real totalLambdaDerivative = Measures::totalLambdaDerivative(this->pg,st);
                outputFile << totalLambdaDerivative << std::endl;
            }

        }
};

}}}}

REGISTER_SIMULATION_STEP(
    ThermodynamicMeasure,ThermodynamicIntegration,
    uammd::structured::SimulationStep::SimulationMeasures::ThermodynamicIntegration
)
