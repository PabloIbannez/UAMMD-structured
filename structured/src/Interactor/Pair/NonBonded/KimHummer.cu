#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"
#include "ParticleData/ExtendedParticleData.cuh"
#include "ParticleData/ParticleGroup.cuh"

#include "Interactor/Pair/PairInteractor.cuh"
#include "Interactor/Pair/NonBonded/NonBonded.cuh"
#include "Interactor/InteractorFactory.cuh"

#include "Interactor/BasicPotentials/NonPolar.cuh"
#include "Interactor/BasicPotentials/DebyeHuckel.cuh"

#include "Interactor/BasicParameters/Pair/LennardJones.cuh"
#include "Utils/ParameterHandler/PairParameterHandler.cuh"

namespace uammd{
namespace structured{
namespace Potentials{
namespace NonBonded{

    namespace KimHummer_ns{
    namespace SasaModel{

        struct A{
            __host__ __device__ inline static real SASAweight(real SASAratio){
                return real(1.0);
            }
        };

        struct B{
            __host__ __device__ inline static real SASAweight(real SASAratio){
                real w = tanhf(real(10)*tanf(SASAratio*real(M_PI_2)));
                return (w<real(0.0))?real(1.0):w;
            }
        };

        struct C{
            __host__ __device__ inline static real SASAweight(real SASAratio){
                real w = tanhf(real(5)*tanf(SASAratio*real(M_PI_2)));
                return (w<real(0.0))?real(1.0):w;
            }
        };

        struct D{
            __host__ __device__ inline static real SASAweight(real SASAratio){
                real w = tanhf(real(2)*tanf(SASAratio*real(M_PI_2)));
                return (w<real(0.0))?real(1.0):w;
            }
        };

        struct E{
            __host__ __device__ inline static real SASAweight(real SASAratio){
                real w = (real(1.0)+tanhf(real(2)*tanf(SASAratio*real(M_PI_2))))/real(2.0);
                return (w<real(0.0))?real(1.0):w;
            }
        };

        struct F{
            __host__ __device__ inline static real SASAweight(real SASAratio){
                real w = (real(1.0)+tanhf(tanf(SASAratio*real(M_PI_2))))/real(2.0);
                return w<real(0.0)?real(1.0):w;
            }
        };
    }}

    struct KimHummer_{

        using ParametersType        = typename BasicParameters::Pairs::LennardJones;
        using ParameterPairsHandler = typename structured::PairParameterHandler<ParametersType>;

        using ParametersPairsIterator = typename ParameterPairsHandler::PairIterator;

        ///////////////////////////

        //Computational data
        struct ComputationalData{

            real* charge;
            real4* pos;

            real* SASAweight;

            Box box;

            ParametersPairsIterator paramPairIterator;

            real ELECOEF;

            real dielectricConstant;
            real debyeLength;

            real zeroEnergy;

            real cutOffNPFactor;
            real cutOffDHFactor;
        };

        //Potential parameters
        struct StorageData{

            std::shared_ptr<ParameterPairsHandler> statParam;

            real ELECOEF;

            real dielectricConstant;
            real debyeLength;

            real zeroEnergy;

            real cutOffNPFactor;
            real cutOffDHFactor;

            real cutOff;
        };

        static __host__ ComputationalData getComputationalData(std::shared_ptr<GlobalData>    gd,
                                                               std::shared_ptr<ParticleGroup> pg,
                                                               const StorageData&  storage,
                                                               const Computables& comp,
                                                               const cudaStream_t& st){

            ComputationalData computational;

            std::shared_ptr<ParticleData> pd = pg->getParticleData();

            computational.pos        = pd->getPos(access::location::gpu, access::mode::read).raw();
            computational.charge     = pd->getCharge(access::location::gpu, access::mode::read).raw();
            computational.SASAweight = pd->getSASAweight(access::location::gpu, access::mode::read).raw();

            computational.box = gd->getEnsemble()->getBox();

            computational.paramPairIterator = storage.statParam->getPairIterator();

            computational.ELECOEF = storage.ELECOEF;

            computational.dielectricConstant = storage.dielectricConstant;
            computational.debyeLength        = storage.debyeLength;

            computational.zeroEnergy = storage.zeroEnergy;

            computational.cutOffNPFactor = storage.cutOffNPFactor;
            computational.cutOffDHFactor = storage.cutOffDHFactor;

            return computational;
        }

        //Storage data reader

        static __host__ StorageData getStorageData(std::shared_ptr<GlobalData>    gd,
                                                   std::shared_ptr<ParticleGroup> pg,
                                                   DataEntry& data){

            StorageData storage;

            storage.ELECOEF = gd->getUnits()->getElectricConversionFactor();

            ///////////////////////////////////////////////////////////

            storage.cutOffNPFactor  = data.getParameter<real>("cutOffNPFactor");
            storage.cutOffDHFactor  = data.getParameter<real>("cutOffDHFactor");

            storage.dielectricConstant = data.getParameter<real>("dielectricConstant");
            storage.debyeLength        = data.getParameter<real>("debyeLength");

            storage.zeroEnergy = data.getParameter<real>("zeroEnergy",real(0.01));

            std::string sasaModel = data.getParameter<std::string>("sasaModel","A");

            ///////////////////////////////////////////////////////////

            storage.statParam = std::make_shared<ParameterPairsHandler>(gd, pg,
                                                                      data);

            ///////////////////////////////////////////////////////////

            auto pairsParam = storage.statParam->getPairParameters();

            real maxSigma = 0.0;
            for(auto p : pairsParam){
                maxSigma=std::max(maxSigma,p.second.sigma);
            }

            real cutOffNP = maxSigma*storage.cutOffNPFactor;
            real cutOffDH = storage.debyeLength*storage.cutOffDHFactor;

            System::log<System::MESSAGE>("[KimHummer] cutOffNP: %f" ,cutOffNP);
            System::log<System::MESSAGE>("[KimHummer] cutOffDH: %f" ,cutOffDH);

            storage.cutOff = std::max(cutOffNP,cutOffDH);

            ///////////////////////////////////////////////////////////

            //Check sasa model
            if(sasaModel != "A"){

                std::shared_ptr<ParticleData> pd = pg->getParticleData();

                if(!pd->isSASAAllocated()){
                    System::log<System::CRITICAL>("[KimHummer] SASA data must be added"
                                                  " for the Kim-Hummer potential if the SASA model"
                                                  " selected is not A");
                }

                auto groupIndex = pg->getIndexIterator(access::location::cpu);

                auto SASA       = pd->getSASA(access::location::cpu, access::mode::read);
                auto SASAweight = pd->getSASAweight(access::location::cpu, access::mode::write);

                fori(0,pg->getNumberParticles()){
                    int  index = groupIndex[i];

                    if        (sasaModel == "B"){
                        SASAweight[index]=KimHummer_ns::SasaModel::B::SASAweight(SASA[index]);
                    } else if (sasaModel == "C"){
                        SASAweight[index]=KimHummer_ns::SasaModel::C::SASAweight(SASA[index]);
                    } else if (sasaModel == "D"){
                        SASAweight[index]=KimHummer_ns::SasaModel::D::SASAweight(SASA[index]);
                    } else if (sasaModel == "E"){
                        SASAweight[index]=KimHummer_ns::SasaModel::E::SASAweight(SASA[index]);
                    } else if (sasaModel == "F"){
                        SASAweight[index]=KimHummer_ns::SasaModel::F::SASAweight(SASA[index]);
                    } else {
                        System::log<System::CRITICAL>("[KimHummerPotential] SASA model: %s is no available ",sasaModel.c_str());
                    }

                    if(SASAweight[index] > real(1.0)){
                        SASAweight[index] = real(1.0);
                    }

                }

            } else {

                auto groupIndex = pg->getIndexIterator(access::location::cpu);

                std::shared_ptr<ParticleData> pd = pg->getParticleData();

                auto SASAweight = pd->getSASAweight(access::location::cpu, access::mode::write);

                fori(0,pg->getNumberParticles()){
                    int  index = groupIndex[i];
                    SASAweight[index]=real(1.0);
                }
            }

            return storage;
        }

        static inline __device__ real energy(const int index_i,const int index_j,
                                             const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));
            const real r2   = dot(rij, rij);

            real e = real(0.0);

            const real sigma = computational.paramPairIterator(index_i,index_j).sigma;

            real cutOffNP2=sigma*computational.cutOffNPFactor;
                 cutOffNP2=cutOffNP2*cutOffNP2;
            if(r2>0 and r2<=cutOffNP2){

                const real eps   = computational.paramPairIterator(index_i,index_j).epsilon;

                e+=BasicPotentials::NonPolar::energy(rij,r2,eps,sigma,computational.zeroEnergy);
            }

            real cutOffDH2 = computational.debyeLength*computational.cutOffDHFactor;
                 cutOffDH2 = cutOffDH2*cutOffDH2;
            const real chgProduct = computational.charge[index_i]*computational.charge[index_j];
            if(r2>0 and r2<=cutOffDH2 and chgProduct != real(0.0)){

                e+=BasicPotentials::DebyeHuckel::DebyeHuckel::energy(rij,r2,
                                                                     computational.ELECOEF,
                                                                     chgProduct,
                                                                     computational.dielectricConstant,
                                                                     computational.debyeLength);
            }

            return computational.SASAweight[index_i]*computational.SASAweight[index_j]*e;
        }

        static inline __device__ real3 force(const int index_i,const int index_j,
                                             const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));
            const real r2   = dot(rij, rij);

            real3 f = make_real3(real(0.0));

            const real sigma = computational.paramPairIterator(index_i,index_j).sigma;

            real cutOffNP2=sigma*computational.cutOffNPFactor;
                 cutOffNP2=cutOffNP2*cutOffNP2;
            if(r2>0 and r2<=cutOffNP2){

                const real eps   = computational.paramPairIterator(index_i,index_j).epsilon;

                f+=BasicPotentials::NonPolar::force(rij,r2,eps,sigma,computational.zeroEnergy);
            }

            real cutOffDH2 = computational.debyeLength*computational.cutOffDHFactor;
                 cutOffDH2 = cutOffDH2*cutOffDH2;
            const real chgProduct = computational.charge[index_i]*computational.charge[index_j];
            if(r2>0 and r2<=cutOffDH2 and chgProduct != real(0.0)){

                f+=BasicPotentials::DebyeHuckel::DebyeHuckel::force(rij,r2,
                                                                    computational.ELECOEF,
                                                                    chgProduct,
                                                                    computational.dielectricConstant,computational.debyeLength);
            }

            return computational.SASAweight[index_i]*computational.SASAweight[index_j]*f;
        }


        static inline __device__ tensor3 hessian(const int index_i,const int index_j,
						 const ComputationalData& computational){

            const real4 posi = computational.pos[index_i];
            const real4 posj = computational.pos[index_j];

            const real3 rij = computational.box.apply_pbc(make_real3(posj)-make_real3(posi));
            const real r2   = dot(rij, rij);

            tensor3 H = tensor3();

            const real sigma = computational.paramPairIterator(index_i,index_j).sigma;

            real cutOffNP2=sigma*computational.cutOffNPFactor;
                 cutOffNP2=cutOffNP2*cutOffNP2;
            if(r2>0 and r2<=cutOffNP2){

                const real eps   = computational.paramPairIterator(index_i,index_j).epsilon;

                H+=BasicPotentials::NonPolar::hessian(rij,r2,eps,sigma,computational.zeroEnergy);
            }

            real cutOffDH2 = computational.debyeLength*computational.cutOffDHFactor;
                 cutOffDH2 = cutOffDH2*cutOffDH2;
            const real chgProduct = computational.charge[index_i]*computational.charge[index_j];
            if(r2>0 and r2<=cutOffDH2 and chgProduct != real(0.0)){

                H+=BasicPotentials::DebyeHuckel::DebyeHuckel::hessian(rij,r2,
								      computational.ELECOEF,
								      chgProduct,
								      computational.dielectricConstant,computational.debyeLength);
            }

            return computational.SASAweight[index_i]*computational.SASAweight[index_j]*H;
        }

    };

    using KimHummer = NonBonded_<KimHummer_>;

}}}}

REGISTER_NONBONDED_INTERACTOR(
    NonBonded,KimHummer,
    uammd::structured::Interactor::PairInteractor<uammd::structured::Potentials::NonBonded::KimHummer>
)
