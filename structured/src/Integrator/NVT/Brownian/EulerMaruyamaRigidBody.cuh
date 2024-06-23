#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"
#include "ParticleGroup/ParticleGroupUtils.cuh"

#include "Integrator/IntegratorBase.cuh"
#include "Integrator/IntegratorFactory.cuh"
#include "Integrator/IntegratorUtils.cuh"

//#include "Interactor/PatchyParticles/PatchyParticlesInteractor.cuh"

namespace uammd{
namespace structured{
namespace Integrator{
namespace NVT{

    namespace Brownian{

        namespace EulerMaruyamaRigidBody_ns {

            VECATTR real4 rotateSecondOrderExpansion(const real4& quat,const real3& omega) {

                real a = quat.x;
                real b = quat.y;
                real c = quat.z;
                real d = quat.w;

                real4 rot;

                real w1=omega.x;
                real w2=omega.y;
                real w3=omega.z;

                real factor = -dot(omega,omega)/real(8.0);

                rot.x = a  +  real(0.5)*(-b*w1-c*w2-d*w3)  +  factor*a;
                rot.y = b  +  real(0.5)*(+a*w1+d*w2-c*w3)  +  factor*b;
                rot.z = c  +  real(0.5)*(-d*w1+a*w2+b*w3)  +  factor*c;
                rot.w = d  +  real(0.5)*(+c*w1-b*w2+a*w3)  +  factor*d;

                return rot;

            }

        }

        class EulerMaruyamaRigidBody : public IntegratorBaseNVT{

            public:

                EulerMaruyamaRigidBody(std::shared_ptr<GlobalData>           gd,
                                       std::shared_ptr<ParticleGroup>        pg,
                                       DataEntry& data,
                                       std::string name):IntegratorBaseNVT(gd,pg,data,name){

                    System::log<System::MESSAGE>("[EulerMaruyamaRigidBody] Created EulerMaruyamaRigidBody integrator \"%s\"",name.c_str());

                    if(!data.isParameterAdded("viscosity")){
                        if(pd->isTranslationalSelfDiffusionAllocated() and
                           pd->isRotationalSelfDiffusionAllocated()){
                            System::log<System::WARNING>("[EulerMaruyamaRigidBody] (%s) viscosity not specified, using translational and rotational self diffusion from particle data",name.c_str());
                        } else {
                            System::log<System::CRITICAL>("[EulerMaruyamaRigidBody] (%s) viscosity not specified and no translational and rotational self diffusion in particle data",name.c_str());
                        }
                    } else {
                        real viscosity = data.getParameter<real>("viscosity");
                        System::log<System::MESSAGE>("[EulerMaruyamaRigidBody] (%s) viscosity : %f", name.c_str(),viscosity);
                        IntegratorUtils::loadMobility(pg,viscosity);
                    }

                }

                ~EulerMaruyamaRigidBody(){}

                void forwardTime() override {

                    System::log<System::DEBUG1>("[EulerMaruyamaRigidBody] (%s) Performing integration step %llu",name.c_str(), this->gd->getFundamental()->getCurrentStep());

                    this->updateForce();
                    this->integrationStep();

                    this->gd->getFundamental()->setCurrentStep(this->gd->getFundamental()->getCurrentStep()+1);
                    this->gd->getFundamental()->setSimulationTime(this->gd->getFundamental()->getSimulationTime()+this->dt);
                }

                void integrationStep(){

                    uint N = this->pg->getNumberParticles();

                    auto mt = this->pd->getTranslationalSelfDiffusion(access::location::gpu, access::mode::read);
                    auto mr = this->pd->getRotationalSelfDiffusion(access::location::gpu, access::mode::read);

                    auto pos   = this->pd->getPos(access::location::gpu, access::mode::readwrite);
                    auto force = this->pd->getForce(access::location::gpu, access::mode::readwrite);

                    auto dir    = this->pd->getDir(access::location::gpu, access::mode::readwrite);
                    auto torque = this->pd->getTorque(access::location::gpu, access::mode::readwrite);

                    auto groupIterator = this->pg->getIndexIterator(access::location::gpu);

                    real* mt_ptr = mt.raw();
                    real* mr_ptr = mr.raw();

                    real4* pos_ptr   = pos.raw();
                    real4* force_ptr = force.raw();

                    real4* dir_ptr    = dir.raw();
                    real4* torque_ptr = torque.raw();

                    uint step_temp = this->gd->getFundamental()->getCurrentStep();
                    uint seed_temp = this->gd->getSystem()->getSeed();
                    real dt_temp   = this->dt;
                    real kBT_temp  = this->kBT;

                    thrust::for_each(thrust::cuda::par.on(this->stream),groupIterator,groupIterator + N,
                            [=] __device__ (int index){
                                                       Saru rng_s(index, step_temp, seed_temp);

                                                       real2 rndAux = make_real2(rng_s.gf(0, real(1.0)));
                                                       real3 noiseTrans = make_real3(rng_s.gf(0, real(1.0)),rndAux.x);
                                                       real3 noiseRot   = make_real3(rng_s.gf(0, real(1.0)),rndAux.y);

                                                       //Translation
                                                       real sigma;
                                                       sigma = sqrt(real(2.0)*kBT_temp*mt_ptr[index]/dt_temp);

                                                       real3 f = make_real3(force_ptr[index]);
                                                       real3 v = mt_ptr[index]*f+sigma*noiseTrans;

                                                       pos_ptr[index].x+=dt_temp*v.x;
                                                       pos_ptr[index].y+=dt_temp*v.y;
                                                       pos_ptr[index].z+=dt_temp*v.z;

                                                       //Rotation
                                                       sigma = sqrt(real(2.0)*kBT_temp*mr_ptr[index]/dt_temp);

                                                       real3 torque = make_real3(torque_ptr[index]);
                                                       real3 domega=mr_ptr[index]*torque+sigma*noiseRot;

                                                       dir[index] = EulerMaruyamaRigidBody_ns::rotateSecondOrderExpansion(dir[index],dt_temp*domega);
                                                       dir[index] = normalize(dir[index]);

                                                       force_ptr[index]  = make_real4(0);
                                                       torque_ptr[index] = make_real4(0);});
                }
        };

        //class EulerMaruyamaRigidBodyPatchesState : public EulerMaruyamaRigidBody {

        //    private:

        //        void updateTransitionProbability() {

        //            uammd::Interactor::Computables comp;
        //            comp.transitionProbability = true;

        //            for(auto forceComp: this->interactors) forceComp->sum(comp,stream);

        //            //CudaSafeCall(cudaStreamSynchronize(stream));
        //            CudaSafeCall(cudaDeviceSynchronize());
        //            CudaCheckError();
        //        }

        //        bool firstStep = true;

        //        std::vector<std::shared_ptr<Interactor::PatchyParticles::PatchyParticles>> patchyParticles;

        //    public:

        //        EulerMaruyamaRigidBodyPatchesState(std::shared_ptr<GlobalData>           gd,
        //                                           std::shared_ptr<ParticleGroup>        pg,
        //                                           DataEntry& data,
        //                                           std::string name):EulerMaruyamaRigidBody(gd,pg,data,name){}

        //        ~EulerMaruyamaRigidBodyPatchesState(){}

        //        void forwardTime() override {

        //            if(firstStep){

        //                for(auto& inte : this->interactors){
        //                    auto ff = std::dynamic_pointer_cast<ForceField>(inte);
        //                    if(ff){
        //                        auto top = ff->getTopology();
        //                        std::map<std::string,std::shared_ptr<typename uammd::Interactor>> topInteractors = top->getInteractors();

        //                        for(auto& topInte : topInteractors){
        //                            auto patchyParticlesInte = std::dynamic_pointer_cast<Interactor::PatchyParticles::PatchyParticles>(topInte.second);
        //                            if(patchyParticlesInte){
        //                                this->patchyParticles.push_back(patchyParticlesInte);
        //                            }
        //                        }
        //                    }

        //                    auto patchyParticlesInte = std::dynamic_pointer_cast<Interactor::PatchyParticles::PatchyParticles>(inte);
        //                    if(patchyParticlesInte){
        //                        patchyParticles.push_back(patchyParticlesInte);
        //                    }
        //                }

        //                for(auto& inte : this->patchyParticles){
        //                    System::log<System::MESSAGE>("[EulerMaruyamaRigidBodyPatchesState] PatchyParticles interactor found, name: %s",inte->getName().c_str());
        //                }

        //                if(patchyParticles.size() == 0){
        //                    System::log<System::CRITICAL>("[EulerMaruyamaRigidBodyPatchesState] No PatchyParticles interactor found, use EulerMaruyamaRigidBody instead");
        //                }

        //                //TODO, check patches group consistency

        //                firstStep = false;
        //            }

        //            System::log<System::DEBUG1>("[EulerMaruyamaRigidBodyPatchesState] (%s) Performing integration step %llu",
        //                                        name.c_str(), this->gd->getFundamental()->getCurrentStep());

        //            this->updateForce();
        //            EulerMaruyamaRigidBody::integrationStep();
        //            this->updateTransitionProbability();
        //            this->integrationStep();

        //            this->gd->getFundamental()->setCurrentStep(this->gd->getFundamental()->getCurrentStep()+1);
        //            this->gd->getFundamental()->setSimulationTime(this->gd->getFundamental()->getSimulationTime()+this->dt);
        //        }

        //        void integrationStep(){
        //            for(auto& patches : this->patchyParticles){

        //                auto patchesPd = patches->getPatchesParticleData();

        //                uint N = patchesPd->getNumParticles();

        //                //Current state
        //                int4* patchesState_ptr = patchesPd->getState(access::location::gpu, access::mode::readwrite).raw();

        //                //Next state
        //                int4* patchesTentativeState_ptr = patchesPd->getTentativeState(access::location::gpu, access::mode::readwrite).raw();
        //                real* transitionProbability_ptr = patchesPd->getTransitionProbability(access::location::gpu, access::mode::readwrite).raw();

        //                uint step_temp = this->gd->getFundamental()->getCurrentStep();
        //                uint seed_temp = this->gd->getSystem()->getSeed();
        //                real dt_temp = this->dt;

        //                thrust::for_each(thrust::cuda::par.on(this->stream),
        //                                 thrust::make_counting_iterator(uint(0)),thrust::make_counting_iterator(N),
        //                        [=] __device__ (int index){

        //                                                   Saru rng_s(index, step_temp, seed_temp);
        //                                                   real rnd = rng_s.d();

        //                                                   real p = transitionProbability_ptr[index];
        //                                                   //Fix to dt
        //                                                   p = real(1.0)-expf(-p*dt_temp);

        //                                                   if(rnd < p){
        //                                                        patchesState_ptr[index] = patchesTentativeState_ptr[index];
        //                                                   }

        //                                                   patchesTentativeState_ptr[index] = {-1,-1,0,0};
        //                                                   transitionProbability_ptr[index] = real(0);

        //                                                   });

        //                cudaStreamSynchronize(this->stream);
        //            }
        //        }

        //};



    }

}}}}

