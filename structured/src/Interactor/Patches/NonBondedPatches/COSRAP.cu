#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Pair/PairInteractor.cuh"
#include "Interactor/Patches/NonBondedPatches/NonBondedPatches.cuh"
#include "Interactor/Patches/PatchesFactory.cuh"

#include "Interactor/BasicPotentials/DistanceSwitch.cuh"
#include "Interactor/BasicPotentials/RotationalAlignmentPotential.cuh"

#include "Interactor/BasicParameters/Pair/SWTRAP.cuh"

#include "Utils/ParameterHandler/PairParameterHandler.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace NonBondedPatches{

    struct COSRAP_{

        using ParametersType        = typename BasicParameters::Pairs::SWTRAP;
        using ParameterPairsHandler = typename structured::PairParameterHandler<ParametersType>;

        using ParametersPairsIterator = typename ParameterPairsHandler::PairIterator;

        ///////////////////////////

        //Computational data
        struct ComputationalData{

            real4* dir;
            int*   patchesParentIndex;

            real4* patchesPos;
            real4* patchesVector;

            Box box;

            ParametersPairsIterator paramPairIterator;
        };

        //Potential parameters
        struct StorageData{

            std::shared_ptr<ParameterPairsHandler> cosrapParam;

            real cutOff;
        };

        static StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                          std::shared_ptr<ParticleGroup> pg,
                                          std::shared_ptr<GlobalData>    patchesGd,
                                          std::shared_ptr<ParticleGroup> patchesPg,
                                          DataEntry& data){

            StorageData storage;

            storage.cosrapParam = std::make_shared<ParameterPairsHandler>(patchesGd,patchesPg,
                                                                          data);

            /////////////////////////////////////////////////////////

            auto pairsParam = storage.cosrapParam->getPairParameters();

            real cutOff = 0.0;
            for(auto p : pairsParam){
                cutOff = std::max(cutOff, p.second.rc);
            }

            System::log<System::MESSAGE>("[COSRAP] cutOff: %f" ,cutOff);
            storage.cutOff = cutOff;

            /////////////////////////////////////////////////////////
            // Check pair parameters
            for(auto p : pairsParam){
                std::string name_i = std::get<0>(p.first);
                std::string name_j = std::get<1>(p.first);
                int batchId        = std::get<2>(p.first);

                real4 Rij_q = p.second.R;
                real4 Rji_q = pairsParam[std::make_tuple(name_j,name_i,batchId)].R;

                // Check if Rij.T=Rji
                tensor3 Rij = MatrixOperations::quat2mat(Rij_q);
                tensor3 Rji = MatrixOperations::quat2mat(Rji_q);

                tensor3 RijT = Rij.transpose();

                real tol   = 1e-6;
                bool check = true;

                if(abs(RijT.xx-Rji.xx)>tol){check=false;}
                if(abs(RijT.xy-Rji.xy)>tol){check=false;}
                if(abs(RijT.xz-Rji.xz)>tol){check=false;}

                if(abs(RijT.yx-Rji.yx)>tol){check=false;}
                if(abs(RijT.yy-Rji.yy)>tol){check=false;}
                if(abs(RijT.yz-Rji.yz)>tol){check=false;}

                if(abs(RijT.zx-Rji.zx)>tol){check=false;}
                if(abs(RijT.zy-Rji.zy)>tol){check=false;}
                if(abs(RijT.zz-Rji.zz)>tol){check=false;}

                if(!check){
                    System::log<System::CRITICAL>("[COSRAP] The transpose of the rotational matrix Rij (%f %f %f %f) is not equal to Rji (%f %f %f %f) for i: %s and j: %s",
                                                  Rij_q.x,Rij_q.y,Rij_q.z,Rij_q.w,
                                                  Rji_q.x,Rji_q.y,Rji_q.z,Rji_q.w,
                                                  name_i.c_str(),name_j.c_str());
                }

            }
            /////////////////////////////////////////////////////////

            return storage;
        }

        static ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                      std::shared_ptr<ParticleGroup> pg,
                                                      std::shared_ptr<GlobalData>    patchesGd,
                                                      std::shared_ptr<ParticleGroup> patchesPg,
                                                      const StorageData&  storage,
                                                      const Computables& comp,
                                                      const cudaStream_t& st){

            ComputationalData computational;

            computational.box               = patchesGd->getEnsemble()->getBox();
            computational.paramPairIterator = storage.cosrapParam->getPairIterator();

            std::shared_ptr<ParticleData> pd = pg->getParticleData();
            computational.dir                = pd->getDir(access::location::gpu, access::mode::read).raw();


            std::shared_ptr<ParticleData> patchesPd = patchesPg->getParticleData();
            computational.patchesParentIndex        = patchesPd->getParentIndex(access::location::gpu, access::mode::read).raw();
            computational.patchesPos                = patchesPd->getPos(access::location::gpu, access::mode::read).raw();
            computational.patchesVector             = patchesPd->getPatchVector(access::location::gpu, access::mode::read).raw();

            return computational;
        }

        //index_i is the current particle
        static inline __device__ EnergyForceTorque energyForceTorque(const int& index_i,const int& index_j,
                                                                     const ComputationalData& computational){

            const real4 patch_i_pos = computational.patchesPos[index_i];
            const real4 patch_j_pos = computational.patchesPos[index_j];

            const tensor3 ori_i = MatrixOperations::quat2mat(computational.dir[computational.patchesParentIndex[index_i]]);
            const tensor3 ori_j = MatrixOperations::quat2mat(computational.dir[computational.patchesParentIndex[index_j]]);

            const real3 rij = computational.box.apply_pbc(make_real3(patch_j_pos)-make_real3(patch_i_pos));

            const auto param = computational.paramPairIterator(index_i,index_j);

            const real E    = param.E;
            const real rc   = param.rc;
            const real Kswt = param.Kswt;
            const real Krap = param.Krap;
                  tensor3 R = MatrixOperations::quat2mat(param.R);

            const real r2 = dot(rij, rij);

            real  e = real(0.0);
            real3 f = make_real3(0.0);
            real3 t = make_real3(0.0);

            if(r2<=rc*rc){
                real4 fe_cos = BasicPotentials::DistanceSwitchCosine::forceEnergy(rij,r2,rc,Kswt);

                real  e_cos = fe_cos.w;
                real3 f_cos = make_real3(fe_cos);
                real3 t_cos = cross(make_real3(computational.patchesVector[index_i]),f_cos);

                EnergyForceTorque eFrcTrq_rap = BasicPotentials::RAP::Stiffness::energyForceTorque(ori_i,ori_j,R,Krap);

                real  e_rap = eFrcTrq_rap.energy;
                //real3 f_rap = make_real3(eFrcTrq_rap.force); //Not used since it is zero
                real3 t_rap = make_real3(eFrcTrq_rap.torque);

                e = E*(e_cos*e_rap-real(1.0));
                f = E*f_cos*e_rap;
                t = E*(t_cos*e_rap+t_rap*e_cos);
            }

            EnergyForceTorque eFrcTrq;

            eFrcTrq.energy = e;
            eFrcTrq.force  = make_real4(make_real3(f),0.0);
            eFrcTrq.torque = make_real4(make_real3(t),0.0);


            return eFrcTrq;

        }

    };

    using COSRAP = NonBondedDynamicallyBondedPatchyParticles_<COSRAP_>;

}}}}

REGISTER_NONBONDED_PATCHES_INTERACTOR(
    NonBondedPatches,COSRAP,
    uammd::structured::Interactor::PairInteractor<uammd::structured::Potentials::NonBondedPatches::COSRAP>
)
