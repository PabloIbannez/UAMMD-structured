#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"
#include "ParticleGroup/ParticleGroupUtils.cuh"

#include "Integrator/IntegratorBase.cuh"
#include "Integrator/IntegratorFactory.cuh"
#include "Integrator/IntegratorUtils.cuh"

namespace uammd{
namespace structured{
namespace Integrator{
namespace NVT{

    namespace Brownian{

        class EulerMaruyama : public IntegratorBaseNVT{

            public:

                EulerMaruyama(std::shared_ptr<GlobalData>           gd,
                              std::shared_ptr<ParticleGroup>        pg,
                              DataEntry& data,
                              std::string name):IntegratorBaseNVT(gd,pg,data,name){

                    System::log<System::MESSAGE>("[EulerMaruyama] Created EulerMaruyama integrator \"%s\"",name.c_str());

                    if(!data.isParameterAdded("viscosity")){
                        if(pd->isTranslationalSelfDiffusionAllocated()){
                            System::log<System::WARNING>("[EulerMaruyama] (%s) viscosity not specified, using translational self diffusion from particle data",name.c_str());
                        } else {
                            System::log<System::CRITICAL>("[EulerMaruyama] (%s) viscosity not specified and no translational self diffusion in particle data",name.c_str());
                        }
                    } else {
                        real viscosity = data.getParameter<real>("viscosity");
                        System::log<System::MESSAGE>("[EulerMaruyama] (%s) viscosity : %f", name.c_str(),viscosity);
                        IntegratorUtils::loadMobility(pg,viscosity);
                    }


                }

                ~EulerMaruyama(){}

                void forwardTime() override {

                    System::log<System::DEBUG1>("[EulerMaruyama] (%s) Performing integration step %llu",name.c_str(), this->gd->getFundamental()->getCurrentStep());

                    this->updateForce();
                    this->integrationStep();

                    this->gd->getFundamental()->setCurrentStep(this->gd->getFundamental()->getCurrentStep()+1);
                    this->gd->getFundamental()->setSimulationTime(this->gd->getFundamental()->getSimulationTime()+this->dt);
                }

                void integrationStep(){

                    uint N = this->pg->getNumberParticles();

                    auto M = this->pd->getTranslationalSelfDiffusion(access::location::gpu, access::mode::read);

                    auto pos   = this->pd->getPos(access::location::gpu, access::mode::readwrite);
                    auto force = this->pd->getForce(access::location::gpu, access::mode::readwrite);

                    auto groupIterator = this->pg->getIndexIterator(access::location::gpu);

                    real* M_ptr = M.raw();

                    real4* pos_ptr   = pos.raw();
                    real4* force_ptr = force.raw();

                    uint step_temp = this->gd->getFundamental()->getCurrentStep();
                    uint seed_temp = this->gd->getSystem()->getSeed();
                    real dt_temp   = this->dt;
                    real kBT_temp  = this->kBT;
                    thrust::for_each(thrust::cuda::par.on(this->stream),groupIterator,groupIterator + N,
                            [=] __device__ (int index){real m     = M_ptr[index];
                                                       real sigma = sqrt(real(2.0)*kBT_temp*m*dt_temp);

                                                       Saru rng_s(index, step_temp, seed_temp);
                                                       const real3 noise = make_real3(rng_s.gf(0, real(1.0)), rng_s.gf(0,real(1.0)).x);

                                                       real3 f = (dt_temp*m)*make_real3(force_ptr[index])+sigma*noise;

                                                       pos_ptr[index].x+=f.x; //pos also stores type at .w
                                                       pos_ptr[index].y+=f.y;
                                                       pos_ptr[index].z+=f.z;

                                                       force_ptr[index] = make_real4(0);});
                }
        };

    }

}}}}

REGISTER_INTEGRATOR(
    Brownian,EulerMaruyama,
    uammd::structured::Integrator::NVT::Brownian::EulerMaruyama
)
