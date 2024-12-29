#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Bonds/BondsInteractor.cuh"
#include "Interactor/Bonds/Bond2/Bond2.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials/RotationalAlignmentPotential.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace Bond2{

    struct RAP_{

        struct ComputationalData{
            real4* pos;
            real4* dir;

            Box    box;
        };

        //Potential parameters
        struct StorageData{};

        struct BondParameters{
            real K;
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
            param.K   = bondParametersMap.at("K");
            real4 R   = bondParametersMap.at("R");
            param.R   = Quat(R);
            return param;
        }

        //Energy and force definition
        static inline __device__ real energy(int index_i, int index_j,
                                             int currentParticleIndex,
                                             const ComputationalData &computational,
                                             const BondParameters &bondParam){



            real    K = bondParam.K;
            tensor3 R = MatrixOperations::quat2mat(bondParam.R);

            int otherParticleIndex;
            if(index_j == currentParticleIndex){
                otherParticleIndex = index_i;
                R = R.transpose();
            } else {
                otherParticleIndex = index_j;
            }

            const tensor3 A = MatrixOperations::quat2mat(computational.dir[currentParticleIndex]);
            const tensor3 B = MatrixOperations::quat2mat(computational.dir[otherParticleIndex]);

            real e = K*BasicPotentials::RAP::energy(A, B, R);

            return e;
        }

        static inline __device__ ForceTorque forceTorque(int index_i, int index_j,
                                                         int currentParticleIndex,
                                                         const ComputationalData &computational,
                                                         const BondParameters &bondParam){

            real    K = bondParam.K;
            tensor3 R = MatrixOperations::quat2mat(bondParam.R);

            int otherParticleIndex;
            if(index_j == currentParticleIndex){
                otherParticleIndex = index_i;
                R = R.transpose();
            } else {
                otherParticleIndex = index_j;
            }

            const tensor3 A = MatrixOperations::quat2mat(computational.dir[currentParticleIndex]);
            const tensor3 B = MatrixOperations::quat2mat(computational.dir[otherParticleIndex]);

            ForceTorque forceTorque;
            forceTorque.force = make_real4(0.0);

            real3 t = K*BasicPotentials::RAP::torque(A, B, R);
            forceTorque.torque = make_real4(t);

            //real e = K*BasicPotentials::RAP::energy(A, B, R);
            //printf("curr %i, i:%i, j:%i, Energy: %f, Torque: %f %f %f\n", currentParticleIndex, index_i, index_j, e, t.x, t.y, t.z);
            //printf("Torque: %f %f %f\n", t.x, t.y, t.z);
            //printf("qR: n:%f x:%f y:%f z:%f, Rmatrix: %f %f %f\n %f %f %f\n %f %f %f\n\n\n",
            //       bondParam.R.n, bondParam.R.v.x, bondParam.R.v.y, bondParam.R.v.z,
            //       R.xx, R.xy, R.xz,
            //       R.yx, R.yy, R.yz,
            //       R.zx, R.zy, R.zz);

            return forceTorque;
        }


    };

    using RAP = Bond2Torque_<RAP_>;

}}}}

REGISTER_BOND_INTERACTOR(
    Bond2,RAP,
    uammd::structured::Interactor::BondsInteractor<uammd::structured::Potentials::Bond2::RAP>
)
