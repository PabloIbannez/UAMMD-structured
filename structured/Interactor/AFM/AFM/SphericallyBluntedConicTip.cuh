#ifndef __AFM_SPHERICALLY_BLUNTED_CONIC_TIP_POTENTIAL__
#define __AFM_SPHERICALLY_BLUNTED_CONIC_TIP_POTENTIAL__

namespace uammd{
namespace structured{
namespace Potentials{
namespace AFM{

    struct SphericallyBluntedConicTip_{

        struct ComputationalData {

            real*  radius;
            real4* pos;

            Box  box;
            real time;

        };

        struct StorageData {};

        struct AFMParameters {

            real sigma;
            real epsilon;

            real K;
            real Kxy;

            real tipAngle;

            real  tipVelocity;
            real3 startChipPosition;

        };

        //Computational data getter
        static ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                      std::shared_ptr<ParticleGroup> pg,
                                                      const StorageData&  storage,
                                                      const Computables& comp){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.radius = pd->getRadius(access::location::gpu, access::mode::read).raw();
            computational.pos    = pd->getPos(access::location::gpu, access::mode::read).raw();

            computational.box    = gd->getEnsemble()->getBox();
            computational.time   = gd->getFundamental()->getSimulationTime();

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

            param.sigma   = afmParametersMap.at("sigma");
            param.epsilon = afmParametersMap.at("epsilon");

            param.K   = afmParametersMap.at("K");
            param.Kxy = afmParametersMap.at("Kxy");

            param.tipAngle = afmParametersMap.at("tipAngle");

            param.tipVelocity       = afmParametersMap.at("tipVelocity");
            param.startChipPosition = afmParametersMap.at("startChipPosition");

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
                chipPos.z    += afmParam.tipVelocity*computational.time;

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
                const real4 pos    = computational.pos[index];
                const real4 tipPos = computational.pos[tipIndex];

                const real3 rij = computational.box.apply_pbc(make_real3(tipPos) - make_real3(pos));
                const real  r2  = dot(rij,rij);
                const real  r   = sqrtf(r2);

                const real Rtip = computational.radius[tipIndex];

                const real sin_theta = sinf(afmParam.tipAngle);
                const real tan_theta = tanf(afmParam.tipAngle);

                const real zSep = -Rtip*sin_theta + tipPos.z;

                if(pos.z <= zSep){

                    const real alpha = abs(r-Rtip);

                    const real sinvr   = afmParam.sigma/alpha;
                    const real sinvr2  = sinvr*sinvr;
                    const real sinvr6  = sinvr2*sinvr2*sinvr2;
                    const real sinvr12 = sinvr6*sinvr6;

                    e += afmParam.epsilon*sinvr12;

                } else {

                    const real A     = sqrtf(rij.x*rij.x + rij.y*rij.y);
                    const real B     = -rij.z + Rtip/sin_theta;
                    const real alpha = abs(A - abs(B)*tan_theta);

                    const real sinvr   = afmParam.sigma/alpha;
                    const real sinvr2  = sinvr*sinvr;
                    const real sinvr6  = sinvr2*sinvr2*sinvr2;
                    const real sinvr12 = sinvr6*sinvr6;

                    e += afmParam.epsilon*sinvr12;

                }

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
                chipPos.z    += afmParam.tipVelocity*computational.time;

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
                const real4 pos    = computational.pos[index];
                const real4 tipPos = computational.pos[tipIndex];

                const real3 rij = computational.box.apply_pbc(make_real3(tipPos) - make_real3(pos));
                const real  r2  = dot(rij,rij);
                const real  r   = sqrtf(r2);

                const real Rtip = computational.radius[tipIndex];

                const real sin_theta = sinf(afmParam.tipAngle);
                const real tan_theta = tanf(afmParam.tipAngle);

                const real zSep = -Rtip*sin_theta + tipPos.z;

                if(pos.z <= zSep){

                    const real alpha = abs(r-Rtip);

                    const real sinvr   = afmParam.sigma/alpha;
                    const real sinvr2  = sinvr*sinvr;
                    const real sinvr6  = sinvr2*sinvr2*sinvr2;
                    const real sinvr12 = sinvr6*sinvr6;

                    const real fmod = afmParam.epsilon*real(12.0)*sinvr12/alpha;
                               f   += -fmod*(rij/r);

                } else {

                    const real A     = sqrtf(rij.x*rij.x + rij.y*rij.y);
                    const real B     = -rij.z + Rtip/sin_theta;
                    const real alpha = abs(A - abs(B)*tan_theta);

                    const real sinvr   = afmParam.sigma/alpha;
                    const real sinvr2  = sinvr*sinvr;
                    const real sinvr6  = sinvr2*sinvr2*sinvr2;
                    const real sinvr12 = sinvr6*sinvr6;

                    const real fmod = afmParam.epsilon*real(12.0)*sinvr12/alpha;

                    f.x += -fmod*(rij.x/A);
                    f.y += -fmod*(rij.y/A);
                    f.z +=  fmod*(tan_theta*B/abs(B));

                }

            }

            return f;

        }

    };

    using SphericallyBluntedConicTip = AFM_<SphericallyBluntedConicTip_>;

}}}}

#endif
