#ifndef __STEEPEST_DESCENT__
#define __STEEPEST_DESCENT__

#include "Integrator/Integrator.cuh"

namespace uammd{
namespace structured{
namespace Minimization{
namespace SteepestDescent{

    class SteepestDescent: public Integrator{

        private:

            cudaStream_t stream;

            std::shared_ptr<GlobalData> gd;

            double h;

            ullint nSteps = 0;
            ullint nStepsPrintProgress;

            real maxObjectiveForce;

            bool initialized = false;
            bool reached     = false;

            /////////////////////////////////////////

            thrust::device_vector<real4> posRef;

            /////////////////////////////////////////

            void copyToRef(){
                auto pos = pd->getPos(access::location::gpu, access::mode::readwrite);
                thrust::copy(thrust::cuda::par.on(stream),pos.begin(), pos.end(), posRef.begin());
                CudaSafeCall(cudaDeviceSynchronize());
                CudaCheckError();
            }

            void copyFromRef(){
                auto pos = pd->getPos(access::location::gpu, access::mode::readwrite);
                thrust::copy(thrust::cuda::par.on(stream),posRef.begin(), posRef.end(), pos.begin());
                CudaSafeCall(cudaDeviceSynchronize());
                CudaCheckError();
            }

            void resetForce(){
                auto force = pd->getForce(access::location::gpu, access::mode::readwrite);
                thrust::fill(thrust::cuda::par.on(stream), force.begin(), force.end(), make_real4(0));
                CudaSafeCall(cudaDeviceSynchronize());
                CudaCheckError();
            }

            void resetEnergy(){
                auto energy = pd->getEnergy(access::location::gpu, access::mode::readwrite);
                thrust::fill(thrust::cuda::par.on(stream), energy.begin(), energy.end(), real(0));
                CudaSafeCall(cudaDeviceSynchronize());
                CudaCheckError();
            }

            void updateForce(){
                for(auto forceComp: interactors) forceComp->sum({.force =true,
                                                                 .energy=false,
                                                                 .virial=false,},0);
                CudaSafeCall(cudaDeviceSynchronize());
                CudaCheckError();
            }

            void updateEnergy(){
                for(auto energyComp: interactors) energyComp->sum({.force =false,
                                                                   .energy=true,
                                                                   .virial=false,},0);
                CudaSafeCall(cudaDeviceSynchronize());
                CudaCheckError();
            }

        public:

            SteepestDescent(std::shared_ptr<GlobalData>    gd,
                            std::shared_ptr<ParticleGroup> pg,
                            DataEntry& data,
                            std::string name):Integrator(pg,name){

                stream = gd->getSystem()->getCudaStream();

                System::log<System::MESSAGE>("[SteepestDescent] Created integrator \"%s\"",name.c_str());

                h                 = data.getParameter<double>("h");
                maxObjectiveForce = data.getParameter<real>("maxObjectiveForce");

                nStepsPrintProgress = data.getParameter<ullint>("nStepsPrintProgress",0);

                System::log<System::MESSAGE>("[SteepestDescent] (%s) h : %f",name.c_str(),h);
                System::log<System::MESSAGE>("[SteepestDescent] (%s) Max Objective Force: %f",name.c_str(),maxObjectiveForce);

                System::log<System::MESSAGE>("[SteepestDescent] (%s) nStepsPrintProgress: %llu",name.c_str(),nStepsPrintProgress);

            }

            void init(){

                //Store the initial positions
                posRef.resize(pd->getNumParticles()); //Copy all

                initialized = true;

            }

            void forwardTime() override {

                CudaCheckError();

                if(not reached){

                    if(not initialized){
                        init();
                    }

                    int N = pg->getNumberParticles();

                    this->resetEnergy();
                    this->updateEnergy();

                    this->resetForce();
                    this->updateForce();

                    real energy   = Measures::totalPotentialEnergy(pg);
                    real maxForce = Measures::maxForce(pg);

                    this->copyToRef();

                    {

                        auto pos   = this->pd->getPos(access::location::gpu, access::mode::readwrite);
                        auto force = this->pd->getForce(access::location::gpu, access::mode::readwrite);

                        auto groupIterator = this->pg->getIndexIterator(access::location::gpu);

                        real4* pos_ptr   = pos.raw();
                        real4* force_ptr = force.raw();

                        double hOverMaxForce = h/maxForce;
                        thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                                [=] __host__ __device__ (int index){

                                    pos_ptr[index].x=double(pos_ptr[index].x+hOverMaxForce*force_ptr[index].x); //pos also stores type at .w
                                    pos_ptr[index].y=double(pos_ptr[index].y+hOverMaxForce*force_ptr[index].y);
                                    pos_ptr[index].z=double(pos_ptr[index].z+hOverMaxForce*force_ptr[index].z);

                                    force_ptr[index] = make_real4(0);
                                });

                    }

                    this->resetEnergy();
                    this->updateEnergy();

                    this->resetForce();
                    this->updateForce();

                    maxForce = Measures::maxForce(pg);

                    if(energy < Measures::totalPotentialEnergy(pg)){
                        h=real(1.001)*h;
                    } else {
                        this->copyFromRef();
                        h=real(0.999)*h;
                    }

                    nSteps++;

                    if(nStepsPrintProgress > 0){
                        if(nSteps%nStepsPrintProgress==0){
                            System::log<System::MESSAGE>("[SteepestDescent] (%s) Steepest descent (%i), maxForce: %f,"
                                                         " maxObjectiveForce: %f , h: %lf",
                                                          name.c_str(),
                                                          nSteps,
                                                          maxForce,
                                                          maxObjectiveForce,h);
                        }
                    }

                    if(maxForce < maxObjectiveForce){
                        reached = true;
                        System::log<System::MESSAGE>("[SteepestDescent] (%s) Steepest descent reached objective force, maxForce: %f,",
                                                      name.c_str(),
                                                      maxForce);
                    }

                }

                this->gd->getFundamental()->setCurrentStep(this->gd->getFundamental()->getCurrentStep()+1);
            }
    };

}}}}

#endif

