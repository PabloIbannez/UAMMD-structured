#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Bonds/BondsInteractor.cuh"
#include "Interactor/Bonds/Bond2/Bond2.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials/Harmonic.cuh"
#include "Interactor/BasicPotentials/RotationalAlignmentPotential.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond2{

    struct HarmonicRAP_{

        struct ComputationalData{
            real4* pos;
            real4* dir;

            Box    box;
        };

        //Potential parameters
        struct StorageData{};

        struct BondParameters{
            real Khrm;
            real r0;
            real3 leftConnection;
            real3 rightConnection;

            real Krap;
            Quat R;// Rotation encoded as a quaternion
        };

        //Computational data getter

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& computables,
                                                               const cudaStream_t& st){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.dir = pd->getDir(access::location::gpu, access::mode::read).raw();

            computational.box = gd->getEnsemble()->getBox();

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>           gd,
                                                   std::shared_ptr<ParticleGroup>        pg,
                                                   DataEntry& data){

            StorageData storage;
            return storage;
        }

        //Bond parameters reader

        template<typename T>
        static __host__ BondParameters processBondParameters(std::shared_ptr<GlobalData> gd,
                                                             std::map<std::string,T>& bondParametersMap){

            BondParameters param;
            param.Khrm            = bondParametersMap.at("Khrm");
            param.r0              = bondParametersMap.at("r0");
            param.leftConnection  = bondParametersMap.at("leftConnection");
            param.rightConnection = bondParametersMap.at("rightConnection");

            param.Krap = bondParametersMap.at("Krap");
            real4 R    = bondParametersMap.at("R");
            param.R    = Quat(R);
            return param;
        }

        //Energy and force definition
        static inline __device__ real energy(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){

            ////////////////////////////////////

            tensor3 R = MatrixOperations::quat2mat(bondParam.R);

            int otherParticleIndex;
            if(index_j == currentParticleIndex){
                otherParticleIndex = index_i;
                R = R.transpose();
            } else {
                otherParticleIndex = index_j;
            }

            ////////////////////////////////////

            const real3 pos_i = make_real3(computational.pos[index_i]);
            const real3 pos_j = make_real3(computational.pos[index_j]);

            const Quat q_i = Quat(computational.dir[index_i]);
            const Quat q_j = Quat(computational.dir[index_j]);

            const tensor3 A = MatrixOperations::quat2mat(computational.dir[currentParticleIndex]);
            const tensor3 B = MatrixOperations::quat2mat(computational.dir[otherParticleIndex]);

            // Harmonic section
            const real3 leftConnection  = bondParam.leftConnection;
            const real3 rightConnection = bondParam.rightConnection;

            const real3 leftPos  = rotateVector(q_i, leftConnection)  + pos_i;
            const real3 rightPos = rotateVector(q_j, rightConnection) + pos_j;

            const real3 rij = computational.box.apply_pbc(rightPos - leftPos);
            const real  r2  = dot(rij, rij);

            const real Khrm   = bondParam.Khrm;
            const real r0     = bondParam.r0;
            const real e_harm = BasicPotentials::Harmonic::energy(rij,r2,Khrm,r0);

            // RAP section
            real  Krap = bondParam.Krap;
            real e_rap = BasicPotentials::RAP::energy(A, B, R, Krap);
            //

            return e_harm + e_rap;
        }

        static inline __device__ ForceTorque forceTorque(int index_i, int index_j,
                                                         int currentParticleIndex,
                                                         const ComputationalData &computational,
                                                         const BondParameters &bondParam){
            ForceTorque forceTorque;

            ////////////////////////////////////

            tensor3 R = MatrixOperations::quat2mat(bondParam.R);

            int otherParticleIndex;
            if(index_j == currentParticleIndex){
                otherParticleIndex = index_i;
                R = R.transpose();
            } else {
                otherParticleIndex = index_j;
            }

            ////////////////////////////////////

            const real3 pos_i = make_real3(computational.pos[index_i]);
            const real3 pos_j = make_real3(computational.pos[index_j]);

            const Quat q_i = Quat(computational.dir[index_i]);
            const Quat q_j = Quat(computational.dir[index_j]);

            const tensor3 A = MatrixOperations::quat2mat(computational.dir[currentParticleIndex]);
            const tensor3 B = MatrixOperations::quat2mat(computational.dir[otherParticleIndex]);

            // Harmonic section
            const real3 leftConnection  = bondParam.leftConnection;
            const real3 rightConnection = bondParam.rightConnection;

            const real3 leftLocal  = rotateVector(q_i, leftConnection);
            const real3 rightLocal = rotateVector(q_j, rightConnection);

            const real3 leftPos  = leftLocal  + pos_i;
            const real3 rightPos = rightLocal + pos_j;

            const real3 rij = computational.box.apply_pbc(rightPos - leftPos);
            const real r2   = dot(rij, rij);

            const real Kharm = bondParam.Khrm;
            const real r0    = bondParam.r0;

            const real3 f_harm = BasicPotentials::Harmonic::force(rij,r2,Kharm,r0);

            if        (currentParticleIndex == index_i){
                forceTorque.force  = make_real4( f_harm);
                forceTorque.torque = make_real4(cross( leftLocal, f_harm));
            } else if (currentParticleIndex == index_j){
                forceTorque.force  = make_real4(-f_harm);
                forceTorque.torque = make_real4(cross(rightLocal,-f_harm));
            }

            // RAP section
            real   Krap = bondParam.Krap;
            real3 t_rap = BasicPotentials::RAP::torque(A, B, R, Krap);

            forceTorque.torque += make_real4(t_rap);

            return forceTorque;
        }
    };

    using HarmonicRAP = Bond2Torque_<HarmonicRAP_>;

}}}}

REGISTER_BOND_INTERACTOR(
    Bond2,HarmonicRAP,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::HarmonicRAP>
)
