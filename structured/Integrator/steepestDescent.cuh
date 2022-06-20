#ifndef __STEEPEST_DESCENT__
#define __STEEPEST_DESCENT__

#include "Integrator/Integrator.cuh"

namespace uammd{
namespace structured{

    class SteepestDescent: public Integrator{

        private:
            
            int N;
            
            int steps=0;

            real maxForce;
            
            thrust::device_vector<real4> posRef;

            cudaStream_t stream;
            
        public:

            struct Parameters{
                double h=-1.0;
            
                int nStepsSteepestDescent=0;
                int nStepsSteepestDescentProgressInterval=0;
                
                real maxObjectiveForce=0;
            };

        private:

            double h;
                
            int nStepsSteepestDescent;
            int nStepsSteepestDescentProgressInterval;
            
            real maxObjectiveForce;
            
            void copyToRef(){
                
                auto pos = pd->getPos(access::location::gpu, access::mode::readwrite);     
                thrust::copy(thrust::cuda::par.on(stream),pos.begin(), pos.end(), posRef.begin());
            }

            void copyFromRef(){
                
                auto pos = pd->getPos(access::location::gpu, access::mode::readwrite);     
                thrust::copy(thrust::cuda::par.on(stream),posRef.begin(), posRef.end(), pos.begin());

            }

            void resetForce(){
                
                auto force = pd->getForce(access::location::gpu, access::mode::readwrite);     
                thrust::fill(thrust::cuda::par.on(stream), force.begin(), force.end(), make_real4(0));
                
                //CudaSafeCall(cudaStreamSynchronize(stream));
            }
            
            void resetEnergy(){
                
                auto energy = pd->getEnergy(access::location::gpu, access::mode::readwrite);     
                thrust::fill(thrust::cuda::par.on(stream), energy.begin(), energy.end(), real(0));
                
                //CudaSafeCall(cudaStreamSynchronize(stream));
            }


            Parameters inputFileToParam(InputFile& in){
                
                Parameters param;
                
                in.getOption("nStepsSteepestDescent",InputFile::Optional)
                              >>param.nStepsSteepestDescent;

                if(param.nStepsSteepestDescent>0){

                    in.getOption("h",InputFile::Required)>>param.h;
                              
                    in.getOption("nStepsSteepestDescentProgressInterval",InputFile::Required)
                                  >>param.nStepsSteepestDescentProgressInterval;

                    in.getOption("maxObjectiveForce",InputFile::Required)
                                  >>param.maxObjectiveForce;
                }

                return param;
            }


        public:
            
            SteepestDescent(shared_ptr<ParticleGroup> pg,
                            InputFile& in,
                            cudaStream_t stream):SteepestDescent(pg,inputFileToParam(in),stream){}

            SteepestDescent(shared_ptr<ParticleGroup> pg,
                            Parameters par,
                            cudaStream_t stream):Integrator(pg,"SteepestDescent"),
                                                 h(par.h),
                                                 nStepsSteepestDescent(par.nStepsSteepestDescent),
                                                 nStepsSteepestDescentProgressInterval(par.nStepsSteepestDescentProgressInterval),
                                                 maxObjectiveForce(par.maxObjectiveForce),
                                                 stream(stream){
                
                sys->log<System::MESSAGE>("[SteepestDescent] h: %lf", 
                                           h);
                sys->log<System::MESSAGE>("[SteepestDescent] nStepsSteepestDescent: %i", 
                                           nStepsSteepestDescent);
                sys->log<System::MESSAGE>("[SteepestDescent] nStepsSteepestDescentProgressInterval: %i",
                                           nStepsSteepestDescentProgressInterval);
                sys->log<System::MESSAGE>("[SteepestDescent] maxObjectiveForce: %f", 
                                           maxObjectiveForce);
                
                N = pg->getNumberParticles();
                
                posRef.resize(pd->getNumParticles()); //Copy all
            
            }

            cudaStream_t getIntegratorStream(){
                return stream;
            }

            void start(){

                if(nStepsSteepestDescent!=0){
                    this->resetForce();
                    this->updateForce();

                    maxForce = Measures::maxForce(pg,stream);
                    cudaStreamSynchronize(stream);

                    sys->log<System::MESSAGE>("[%s] Steepest descent (init), maxForce: %f,"
                                              " maxObjectiveForce: %f , h: %lf",
                                               name.c_str(),
                                               maxForce,
                                               maxObjectiveForce,h);
                    
                    if(maxForce < maxObjectiveForce){
                        sys->log<System::MESSAGE>("[%s] Initial force is lower than maxObjectiveForce, doing nothing ...",
                                                   name.c_str());
                        return;
                    }
                            

                    int t;
                    for(t=0;t<nStepsSteepestDescent;t++){
                        this->forwardTime();
                        if(maxForce < maxObjectiveForce){
                            break;
                        }
                        if(t%nStepsSteepestDescentProgressInterval==0){
                            sys->log<System::MESSAGE>("[%s] Steepest descent (%i/%i), maxForce: %f,"
                                                      " maxObjectiveForce: %f , h: %lf",
                                                       name.c_str(),
                                                       t,nStepsSteepestDescent,
                                                       maxForce,
                                                       maxObjectiveForce,h);
                        }

                    }
                    
                    sys->log<System::MESSAGE>("[%s] Steepest descent (%i/%i), maxForce: %f,"
                                              " maxObjectiveForce: %f , h: %lf",
                                              name.c_str(),
                                              t,nStepsSteepestDescent,
                                              maxForce,
                                              maxObjectiveForce,h);
                }

            }
            
            void forwardTime() {
                
                CudaCheckError();

                steps++;
                sys->log<System::DEBUG1>("[%s] Performing Steepest descent step %d", name.c_str(), steps);

                this->resetEnergy();
                this->updateEnergy();
                
                this->updateForce();
                cudaStreamSynchronize(stream);

                real Eprev    = Measures::totalPotentialEnergy(pg,stream);
                     maxForce = Measures::maxForce(pg,stream);
                cudaStreamSynchronize(stream);

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
                cudaStreamSynchronize(stream);
                
                real Epost = Measures::totalPotentialEnergy(pg,stream);
                cudaStreamSynchronize(stream);
                //sys->log<System::MESSAGE>("[SteepestDescent] Epost/N:%f Eprev/N:%f", 
                //                           Epost/N, Eprev/N);

                if(Epost < Eprev){
                    h=real(1.001)*h;
                } else {
                    this->copyFromRef();
                    h=real(0.999)*h;
                }
            
            }
            
            void updateForce(){
                for(auto forceComp: interactors) forceComp->sum({.force =true,
                                                                 .energy=false,
                                                                 .virial=false,},stream);
                //CudaSafeCall(cudaStreamSynchronize(stream));
            }

            void updateEnergy(){
                for(auto energyComp: interactors) energyComp->sum({.force =false,
                                                                   .energy=true,
                                                                   .virial=false,},stream);
                //CudaSafeCall(cudaStreamSynchronize(stream));
            }
    };

}}

#endif

