#ifndef __AFM_SPHERICAL_TIP_POTENTIAL__
#define __AFM_SPHERICAL_TIP_POTENTIAL__

namespace uammd{
namespace structured{
namespace Potentials{
namespace AFM{

    struct SphericalTip_{

        struct ComputationalData {

            real*  radius;
            real4* pos;

            Box  box;

            real   dt;
            ullint step;
        };

        struct StorageData {};

        struct AFMParameters {

            real  epsilon;
            real  sigma;

            real  K;
            real  Kxy;

            real  tipVelocity;

            real3  startChipPosition;
            ullint indentationStartStep;
            ullint indentationBackwardStep;
        };

        //Computational data getter
        static ComputationalData getComputationalData(std::shared_ptr<GlobalData>           gd,
                                                      std::shared_ptr<ParticleGroup>        pg,
                                                      const StorageData&  storage,
                                                      const Computables& comp,
                                                      const cudaStream_t& st){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.radius = pd->getRadius(access::location::gpu, access::mode::read).raw();
            computational.pos    = pd->getPos(access::location::gpu, access::mode::read).raw();

            computational.box    = gd->getEnsemble()->getBox();

            computational.dt    = gd->getFundamental()->getTimeStep();
            computational.step  = gd->getFundamental()->getCurrentStep();

            return computational;
        }

        //Storage data reader
        static StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                          std::shared_ptr<ParticleGroup> pg,
                                          DataEntry& data){

            StorageData storage;
            return storage;
        }


        template<typename T>
        static AFMParameters processAFMParameters(std::map<std::string,T>& afmParametersMap){

            AFMParameters param;

            param.epsilon = afmParametersMap.at("epsilon");
            param.sigma   = afmParametersMap.at("sigma");

            param.K   = afmParametersMap.at("K");
            param.Kxy = afmParametersMap.at("Kxy");

            param.tipVelocity       = afmParametersMap.at("tipVelocity");

            param.startChipPosition    = afmParametersMap.at("startChipPosition");

            param.indentationStartStep    = afmParametersMap.at("indentationStartStep");
            param.indentationBackwardStep = afmParametersMap.at("indentationBackwardStep");

            if(param.indentationBackwardStep == 0){
                //Set to maximum value
                param.indentationBackwardStep = std::numeric_limits<ullint>::max();
            }

            System::log<System::MESSAGE>("[SphericalTip] AFM epsilon: %f", param.epsilon);
            System::log<System::MESSAGE>("[SphericalTip] AFM sigma: %f", param.sigma);

            System::log<System::MESSAGE>("[SphericalTip] AFM K: %f", param.K);
            System::log<System::MESSAGE>("[SphericalTip] AFM Kxy: %f", param.Kxy);

            System::log<System::MESSAGE>("[SphericalTip] AFM tipVelocity: %f", param.tipVelocity);

            System::log<System::MESSAGE>("[SphericalTip] AFM startChipPosition: %f %f %f",
                                         param.startChipPosition.x,
                                         param.startChipPosition.y,
                                         param.startChipPosition.z);

            System::log<System::MESSAGE>("[SphericalTip] AFM indentationStartStep: %llu"   , param.indentationStartStep);
            System::log<System::MESSAGE>("[SphericalTip] AFM indentationBackwardStep: %llu", param.indentationBackwardStep);

            //Check if indentationBackwardStep is larger than indentationStartStep
            if(param.indentationBackwardStep < param.indentationStartStep){
                System::log<System::CRITICAL>("[SphericalTip] indentationBackwardStep must be larger than indentationStartStep");
            }


            return param;
        }

        //Energy and force definition

        static inline __device__ real energy(int index,
                                             int tipIndex,
                                             const ComputationalData &computational,
                                             const AFMParameters &afmParam){

            real e = real(0.0);
            if(index == tipIndex){
                const real4 tipPos    = computational.pos[tipIndex];

                real3 chipPos = afmParam.startChipPosition;
                if(computational.step >= afmParam.indentationStartStep and computational.step < afmParam.indentationBackwardStep){
                    real time = computational.dt*(computational.step - afmParam.indentationStartStep);
                    chipPos.z    += afmParam.tipVelocity*time;
                }
                if(computational.step >= afmParam.indentationBackwardStep){
                    real time              = computational.dt*(computational.step - afmParam.indentationBackwardStep);
                    real timeUntilBackward = computational.dt*(afmParam.indentationBackwardStep - afmParam.indentationStartStep);
                    chipPos.z    += afmParam.tipVelocity*timeUntilBackward - afmParam.tipVelocity*time;
                }

                //Tip-chip interaction

                const real3 dr = chipPos - make_real3(tipPos);

                real3  K;
                K.x = afmParam.Kxy;
                K.y = afmParam.Kxy;
                K.z = afmParam.K;

                const real3 r0 = make_real3(0.0);
                e += BasicPotentials::HarmonicAnisotropic::energy(dr,K,r0);

            } else {
                //Tip-particle interaction
                const real3 pos    = make_real3(computational.pos[index]);
                const real3 tipPos = make_real3(computational.pos[tipIndex]);

                real eps = afmParam.epsilon;

                real A = real( 0.0);
                real B = real( 1.0);

                if( eps < real(0.0) ){
                    eps = -eps;
                    A = real( 1.0);
                    B = real(-2.0);
                }

                e += BasicPotentials::Tip::GenericSphericalTip::energy(pos,tipPos,A,B,
                                                                       computational.radius[tipIndex],
                                                                       eps,afmParam.sigma,
                                                                       computational.box);

            }

            return e;
        }

        static inline __device__ real3 force(int index,
                                             int tipIndex,
                                             const ComputationalData &computational,
                                             const AFMParameters &afmParam){

            real3 f = make_real3(0.0);
            if(index == tipIndex){
                const real4 tipPos    = computational.pos[tipIndex];

                real3 chipPos = afmParam.startChipPosition;
                if(computational.step >= afmParam.indentationStartStep and computational.step < afmParam.indentationBackwardStep){
                    real time = computational.dt*(computational.step - afmParam.indentationStartStep);
                    chipPos.z    += afmParam.tipVelocity*time;
                }
                if(computational.step >= afmParam.indentationBackwardStep){
                    real time              = computational.dt*(computational.step - afmParam.indentationBackwardStep);
                    real timeUntilBackward = computational.dt*(afmParam.indentationBackwardStep - afmParam.indentationStartStep);
                    chipPos.z    += afmParam.tipVelocity*timeUntilBackward - afmParam.tipVelocity*time;
                }

                //Tip-chip interaction

                const real3 dr = chipPos - make_real3(tipPos);

                real3  K;
                K.x = afmParam.Kxy;
                K.y = afmParam.Kxy;
                K.z = afmParam.K;

                const real3 r0 = make_real3(0.0);

                f += BasicPotentials::HarmonicAnisotropic::force(dr,K,r0);

            } else {

                //Tip-particle interaction
                const real3 pos    = make_real3(computational.pos[index]);
                const real3 tipPos = make_real3(computational.pos[tipIndex]);

                real eps = afmParam.epsilon;

                real A = real( 0.0);
                real B = real( 1.0);

                if( eps < real(0.0) ){
                    eps = -eps;
                    A = real( 1.0);
                    B = real(-2.0);
                }

                f += BasicPotentials::Tip::GenericSphericalTip::force(pos,tipPos,A,B,
                                                                      computational.radius[tipIndex],
                                                                      eps,afmParam.sigma,
                                                                      computational.box);

            }

            return f;

        }

    };

    using SphericalTip = AFM_<SphericalTip_>;

}}}}

#endif //__AFM_SPHERICAL_TIP_POTENTIAL__
