#include "uammd.cuh"
#include "utils/quaternion.cuh"

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Pair/PairInteractor.cuh"
#include "Interactor/Pair/NonBonded/NonBonded.cuh"
#include "Interactor/InteractorFactory.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace NonBonded{

    namespace DipolarMagnetic_ns{
        __device__ real3 computeForce(real3 rij, real3 mi, real3 mj, real prefactor,
                real invr2, real invr5){
            real mi_dot_rij = dot(rij, mi);
            real mj_dot_rij = dot(rij, mj);
            real3 term1 = mi_dot_rij*mj;
            real3 term2 = mj_dot_rij*mi;
            real3 term3 = dot(mi,mj)*rij;
            real3 term4 = real(-5.0)*(mi_dot_rij*mj_dot_rij)*rij*invr2;
            return -prefactor*(term1+term2+term3+term4)*invr5;
        }
        __device__ real3 computeField(real3 rij, real3 mi, real prefactor, real invr3, real invr5){
            return prefactor*(real(3.0)*dot(rij, mi)*rij*invr5-mi*invr3);
        }
    }

    struct DipolarMagnetic_{

        //Computational data
        struct ComputationalData{

            real4* pos;
            real4* dir;
            real4* magnetization;

            Box box;

            real permeability;

            real cutOff;
        };

        //Potential parameters
        struct StorageData{

            real permeability;
            real cutOff;
        };

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>           gd,
                std::shared_ptr<ParticleGroup>        pg,
                const StorageData&  storage,
                const Computables& comp,
                const cudaStream_t& st){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.dir = pd->getDir(access::location::gpu, access::mode::read).raw();
            computational.magnetization = pd->getMagnetization(access::location::gpu, access::mode::read).raw();

            computational.box = gd->getEnsemble()->getBox();

            computational.permeability = storage.permeability;
            computational.cutOff = storage.cutOff;


            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                std::shared_ptr<ParticleGroup>        pg,
                DataEntry& data){
            StorageData storage;

            storage.permeability = data.getParameter<real>("permeability");
            storage.cutOff = data.getParameter<real>("cutOff");

            System::log<System::MESSAGE>("[Dipolar Magnetic] cutOff: %f" ,storage.cutOff);
            System::log<System::MESSAGE>("[Dipolar Magnetic] permeability: %f" ,storage.permeability);

            return storage;
        }


        static inline __device__ real energy(const int index_i,const int index_j,
                const ComputationalData& computational){
            return real(0.0);
        }

        static inline __device__ ForceTorqueMagneticField forceTorqueMagneticField(const int index_i,const int index_j,
                const ComputationalData& computational){
            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            ForceTorqueMagneticField ftf;
            ftf.force  = make_real4(0.0);
            ftf.torque = make_real4(0.0);
            ftf.magneticField = make_real4(0.0);

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));
            const real r2 = dot(rij, rij);

            if (r2 == real(0.0) or r2>computational.cutOff*computational.cutOff){
                return ftf;
            }

            const Quat diri = computational.dir[index_i];
            const Quat dirj = computational.dir[index_j];
            const real4 mi_and_Mi = computational.magnetization[index_i];
            const real4 mj_and_Mj = computational.magnetization[index_j];

            real3 mi = rotateVector(diri, make_real3(mi_and_Mi));
            real3 mj = rotateVector(dirj, make_real3(mj_and_Mj));
            const real invr2 = real(1.0)/r2;
            const real invr3 = sqrtf(invr2)*invr2;
            const real invr5 = invr3*invr2;
            real prefactorField = real(0.25)*computational.permeability*mj_and_Mj.w/real(M_PI);
            real prefactorForce = real(3.0)*prefactorField*mi_and_Mi.w;
            real3 f = DipolarMagnetic_ns::computeForce(rij, mi, mj, prefactorForce,
                    invr2, invr5);
            real3 field = DipolarMagnetic_ns::computeField(rij, mj, prefactorField,
                    invr3, invr5);
            real3 t = mi_and_Mi.w*cross(mi, field);
            ftf.force = make_real4(f,real(0));
            ftf.torque = make_real4(t, real(0));
            ftf.magneticField = make_real4(field, real(0));
            return ftf;

        }

        static inline __device__ real4 magneticField(const int index_i,const int index_j,
                const ComputationalData& computational){
            return forceTorqueMagneticField(index_i, index_j, computational).magneticField;
        }

        static inline __device__ ForceTorque forceTorque(const int index_i,const int index_j,
                const ComputationalData& computational){

            ForceTorqueMagneticField ftf = forceTorqueMagneticField(index_i, index_j, computational);

            ForceTorque ft;

            ft.force  = ftf.force;
            ft.torque = ftf.torque;

            return ft;
        }
    };

    using DipolarMagnetic = NonBonded_<DipolarMagnetic_>;

}}}}

REGISTER_NONBONDED_INTERACTOR(
    NonBonded,DipolarMagnetic,
    uammd::structured::Interactor::PairInteractor<uammd::structured::Potentials::NonBonded::DipolarMagnetic>
)
