#ifndef __GJF_NVT__
#define __GJF_NVT__

namespace uammd{
namespace structured{
namespace NVT{

    namespace Langevin{

        class GJF: public IntegratorBasicNVT{

            private:

                real frictionConstant;

                ///////////////////////////

                bool initialized = false;
                bool resetVelocities = true;

            public:

                GJF(std::shared_ptr<GlobalData>    gd,
                    std::shared_ptr<ParticleGroup> pg,
                    DataEntry& data,
                    std::string name):IntegratorBasicNVT(gd,pg,data,name){

                    System::log<System::MESSAGE>("[GJF] Created GJF integrator \"%s\"",name.c_str());

                    frictionConstant  = data.getParameter<real>("frictionConstant");
                    resetVelocities   = data.getParameter<bool>("resetVelocities",true);

                    System::log<System::MESSAGE>("[GJF] (%s) frictionConstant : %f",name.c_str(),frictionConstant);
                    System::log<System::MESSAGE>("[GJF] (%s) resetVelocities : %i",name.c_str(),resetVelocities);
                }

                void init(){

                    System::log<System::MESSAGE>("[GJF] Performing initialization step");

                    //Inital velocities for some T
                    if(resetVelocities){
                        IntegratorUtils::generateVelocity(this->pg,
                                                          this->kBT,
                                                          this->gd->getSystem()->getSeed(),
                                                          this->stream);
                    }

                    initialized = true;

                }

                void forwardTime() override{

                    if(not initialized){
                        this->init();
                        this->updateForce();
                    }

                    System::log<System::DEBUG1>("[GJF] (%s) Performing integration step %llu",name.c_str(), this->gd->getFundamental()->getCurrentStep());

                    this->integrationStep1();
                    this->updateForce();
                    this->integrationStep2();

                    this->gd->getFundamental()->setCurrentStep(this->gd->getFundamental()->getCurrentStep()+1);
                    this->gd->getFundamental()->setSimulationTime(this->gd->getFundamental()->getSimulationTime()+this->dt);
                }

                void integrationStep1(){

                    int N = this->pg->getNumberParticles();

                    uint step_temp = this->gd->getFundamental()->getCurrentStep();
                    uint seed_temp = this->gd->getSystem()->getSeed();

                    real comboGJ = 0.5*this->dt*frictionConstant;
                    real sigmaGaussPrefactor = sqrt(2*frictionConstant*this->kBT*this->dt);
                    real bGJ = 1.0/(1.0+comboGJ);
                    real aGJ = (1.0-comboGJ)/(1.0+comboGJ);
                    real cr2 = bGJ*this->dt;
                    real cr3 = bGJ*this->dt*this->dt*0.5;
                    real cr4 = bGJ*this->dt*0.5;
                    real vr1 = aGJ;
                    real vr2 = this->dt*0.5*aGJ;
                    real vr4 = bGJ;

                    {
                        auto mass  = this->pd->getMass(access::location::gpu, access::mode::read);

                        auto pos   = this->pd->getPos(access::location::gpu, access::mode::readwrite);
                        auto vel   = this->pd->getVel(access::location::gpu, access::mode::readwrite);
                        auto force = this->pd->getForce(access::location::gpu, access::mode::readwrite);

                        auto groupIterator = this->pg->getIndexIterator(access::location::gpu);

                        real*  mass_ptr = mass.raw();

                        real4* pos_ptr   = pos.raw();
                        real3* vel_ptr   = vel.raw();
                        real4* force_ptr = force.raw();

                        thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                                [=] __host__ __device__ (int index){
                                const real m = mass_ptr[index];

                                real sigma   = sqrt(m)*sigmaGaussPrefactor;

                                Saru rng_s(index, step_temp, seed_temp);
                                real3 noise = make_real3(rng_s.gf(real(0.0), real(1.0)), rng_s.gf(real(0.0),real(1.0)).x);

                                noise = sigma*noise;

                                real3 f = make_real3(force_ptr[index]);

                                pos_ptr[index].x=pos_ptr[index].x+cr2*vel_ptr[index].x+(cr3/m)*f.x+(cr4/m)*noise.x; //pos also stores type at .w
                                pos_ptr[index].y=pos_ptr[index].y+cr2*vel_ptr[index].y+(cr3/m)*f.y+(cr4/m)*noise.y;
                                pos_ptr[index].z=pos_ptr[index].z+cr2*vel_ptr[index].z+(cr3/m)*f.z+(cr4/m)*noise.z;

                                vel_ptr[index].x=vr1*vel_ptr[index].x+(vr2/m)*f.x+(vr4/m)*noise.x;
                                vel_ptr[index].y=vr1*vel_ptr[index].y+(vr2/m)*f.y+(vr4/m)*noise.y;
                                vel_ptr[index].z=vr1*vel_ptr[index].z+(vr2/m)*f.z+(vr4/m)*noise.z;

                                force_ptr[index] = make_real4(0);});
                    }
                    CudaSafeCall(cudaStreamSynchronize(stream));
                }

                void integrationStep2(){

                    int N = this->pg->getNumberParticles();

                    real vr3 = this->dt*0.5;

                    {
                        auto mass  = this->pd->getMass(access::location::gpu, access::mode::read);

                        auto vel   = this->pd->getVel(access::location::gpu, access::mode::readwrite);
                        auto force = this->pd->getForce(access::location::gpu, access::mode::readwrite);

                        auto groupIterator = this->pg->getIndexIterator(access::location::gpu);

                        real*  mass_ptr = mass.raw();

                        real3* vel_ptr   = vel.raw();
                        real4* force_ptr = force.raw();

                        thrust::for_each(thrust::cuda::par.on(stream),groupIterator,groupIterator + N,
                                [=] __host__ __device__ (int index){
                                const real m = mass_ptr[index];

                                real3 f = make_real3(force_ptr[index]);

                                vel_ptr[index].x=vel_ptr[index].x+(vr3/m)*f.x;
                                vel_ptr[index].y=vel_ptr[index].y+(vr3/m)*f.y;
                                vel_ptr[index].z=vel_ptr[index].z+(vr3/m)*f.z;

                                });
                    }
                    CudaSafeCall(cudaStreamSynchronize(stream));

                }

        };

    }

}}}

#endif
