#pragma once

#include "uammd.cuh"

#include "System/ExtendedSystem.cuh"
#include "GlobalData/GlobalData.cuh"

#include "Interactor/Interactor.cuh"

#include "Definitions/SFINAE.cuh"
#include "Utils/Containers/SetUtils.cuh"

#include "Definitions/Computations.cuh"
#include "Definitions/Types.cuh"

namespace uammd{
namespace structured{
namespace Interactor{

    namespace BondsInteractor_ns{

        template<class ComputationalData, class BondParameters, class PotentialTransverser>
        __global__ void transverseBondListThreadPerParticle(const int      nParticlesWithBonds,
                                                            const ComputationalData computational,
                                                            const int*  __restrict__ id2index,
                                                            const int*  __restrict__ partLocalIndex2id,
                                                            const BondParameters* __restrict__ bondList,
                                                            const int*  __restrict__ partLocalIndex2bondsStart,
                                                            const int*  __restrict__ partLocalIndex2bonds,
                                                            const int*  __restrict__ partLocalIndex2nBondsInvolved,
                                                            PotentialTransverser potTransverser) {

            int partLocalIndex = blockIdx.x*blockDim.x + threadIdx.x;

            if(partLocalIndex>=nParticlesWithBonds){return;}

            const int start = partLocalIndex2bondsStart[partLocalIndex];
            const int nBondsInvolved = partLocalIndex2nBondsInvolved[partLocalIndex];

            typename PotentialTransverser::resultType quantityLocal = potTransverser.zero();

            const int currentIndex = id2index[partLocalIndex2id[partLocalIndex]];
            for(int b = 0; b < nBondsInvolved; b++) {

                const int bondIndex = partLocalIndex2bonds[start+b];

                const BondParameters bondParam = bondList[bondIndex];

                potTransverser.accumulate(quantityLocal,potTransverser.compute(currentIndex,computational,bondParam));
            }

            potTransverser.set(currentIndex,quantityLocal);
        }
    }

    template<class BondType, int THREADS_PER_BLOCK = 256>
    class BondsInteractor: public Interactor{

        private:

            int nPart = BondType::nPart;

            std::shared_ptr<GlobalData> gd;
            std::shared_ptr<BondType>   pot;

            ////////////////////////////////////

            int nParticlesWithBonds = 0;

            //Bond list is an array of size nBonds * nPart.
            //Each bond is represented by a list of particles.
            thrust::device_vector<typename BondType::BondParameters> bondList_d;

            thrust::device_vector<int> partLocalIndex2id_d;

            thrust::device_vector<int> partLocalIndex2bondsStart_d;
            thrust::device_vector<int> partLocalIndex2bonds_d;
            thrust::device_vector<int> partLocalIndex2nBondsInvolved_d;

            ////////////////////////////////////
            //Warnings

            bool warningEnergy           = false;
            bool warningForce            = false;
            bool warningLambdaDerivative = false;
            bool warningStress           = false;
            bool warningMagneticField    = false;
            bool warningHessian          = false;
            bool warningPairwiseForces   = false;

        public:

            BondsInteractor(std::shared_ptr<GlobalData>           gd,
                            std::shared_ptr<ParticleGroup>        pg,
                            DataEntry& data,
                            std::string name):Interactor(pg,"BondsInteractor: \"" +name+"\""),
                                              gd(gd){

                this->pot = std::make_shared<BondType>(gd,pg,data);

                // Get the labels for every particle in the bond
                std::vector<std::string> partLabels = pot->getParticleBondLabels();

                //Check partLabels has the correct size
                if(partLabels.size() != nPart){
                    System::log<System::CRITICAL>("[BondsInteractor] (%s) particle labels size is not correct. Expected %d, but got %d",name.c_str(),nPart,partLabels.size());
                }

                //Check if partLabels are all different
                for(int i = 0; i < nPart; i++){
                    for(int j = i+1; j < nPart; j++){
                        if(partLabels[i] == partLabels[j]){
                            System::log<System::CRITICAL>("[BondsInteractor] (%s) particle labels are not all different. %s is repeated",name.c_str(),partLabels[i].c_str());
                        }
                    }
                }

                //Check if partLabels is contained in data labels
                {
                  std::vector<std::string> dataLabels = data.getLabels();
                  for(int i = 0; i < nPart; i++){
                    if(std::find(dataLabels.begin(), dataLabels.end(), partLabels[i]) == dataLabels.end()){
                      System::log<System::CRITICAL>("[BondsInteractor] (%s) particle label %s is not contained in data labels",name.c_str(),partLabels[i].c_str());
                    }
                  }
                }

                //Print the given labels, order matters!!
                {
                    std::string labels = "";
                    for(int i = 0; i < nPart; i++){
                      if(i == nPart - 1){
                        labels += partLabels[i];
                      }else{
                        labels += partLabels[i] + "----";
                      }
                    }
                    System::log<System::MESSAGE>("[BondsInteractor] (%s) Particles in the bond: %s",name.c_str(),labels.c_str());
                }

                //Load bonds from data into bondList_h
                auto bondsData = data.getDataMap();

                if(bondsData.size() == 0){
                    System::log<System::WARNING>("[BondsInteractor] (%s) No bonds found in data. This interactor will do nothing.",name.c_str());
                } else {
                    //Declare host temporary vectors
                    thrust::host_vector<typename BondType::BondParameters> bondList_h;
                    thrust::host_vector<int> partLocalIndex2id_h;

                    thrust::host_vector<int> partLocalIndex2bondsStart_h;
                    thrust::host_vector<int> partLocalIndex2bonds_h;
                    thrust::host_vector<int> partLocalIndex2nBondsInvolved_h;


                    //Create a set of the particles which belong to one bond at least
                    std::set<int> bondedParticles;
                    for(auto bondData : bondsData){
                      for(auto partLabel : partLabels){
                        bondedParticles.insert(int(bondData.at(partLabel)));
                      }
                    }

                    nParticlesWithBonds = bondedParticles.size();

                    //Sort the set
                    std::vector<int> bondedParticlesSorted(bondedParticles.begin(), bondedParticles.end());
                    std::sort(bondedParticlesSorted.begin(), bondedParticlesSorted.end());

                    partLocalIndex2id_h.resize(bondedParticlesSorted.size());

                    std::map<int, int> id2partLocalIndex;
                    //Fill id2partLocalIndex and partLocalIndex2id_h
                    for(int i = 0; i < bondedParticlesSorted.size(); i++){
                      id2partLocalIndex[bondedParticlesSorted[i]] = i;
                      partLocalIndex2id_h[i] = bondedParticlesSorted[i];
                    }

                    //Resize bondList_h
                    bondList_h.resize(data.getDataSize()*nPart);

                    std::map<int, std::vector<int>> partLocalIndex2bonds_tmp;
                    for(auto bondData : bondsData){
                      try{
                        bondList_h.push_back(pot->processBondParameters(bondData));
                      } catch (std::exception& e){
                        System::log<System::CRITICAL>("[BondsInteractor] (%s) Error while processing bond parameters: %s",name.c_str(),e.what());
                      }
                      for(auto partLabel : partLabels){
                        partLocalIndex2bonds_tmp[id2partLocalIndex[bondData.at(partLabel)]].push_back(bondList_h.size() - 1);
                      }
                    }

                    //Resize partLocalIndex2bondsStart_h and partLocalIndex2bonds_h
                    partLocalIndex2bondsStart_h.resize(bondedParticlesSorted.size());
                    partLocalIndex2bonds_h.resize(bondList_h.size()*nPart);
                    //Resize partLocalIndex2nBondsInvolved_h and fill it with zeros
                    partLocalIndex2nBondsInvolved_h.resize(bondedParticlesSorted.size(), 0);

                    //Fill partLocalIndex2bondsStart_h, partLocalIndex2bonds_h and partLocalIndex2nBondsInvolved_h
                    int bondStart = 0;
                    for(int i = 0; i < bondedParticlesSorted.size(); i++){
                      partLocalIndex2bondsStart_h[i] = bondStart;
                      for(int j = 0; j < partLocalIndex2bonds_tmp[i].size(); j++){
                        partLocalIndex2bonds_h[bondStart + j] = partLocalIndex2bonds_tmp[i][j];
                      }
                      bondStart += partLocalIndex2bonds_tmp[i].size();
                      partLocalIndex2nBondsInvolved_h[i] = partLocalIndex2bonds_tmp[i].size();
                    }

                    //Copy data to device
                    bondList_d = bondList_h;
                    partLocalIndex2id_d = partLocalIndex2id_h;

                    partLocalIndex2bondsStart_d     = partLocalIndex2bondsStart_h;
                    partLocalIndex2bonds_d          = partLocalIndex2bonds_h;
                    partLocalIndex2nBondsInvolved_d = partLocalIndex2nBondsInvolved_h;

                    //The way all vectors are used can be checked in the cuda kernel
                }
            }

            std::shared_ptr<BondType> getPotential(){
                return pot;
            }

            void sum(Computables comp,cudaStream_t st) override {

                if(nParticlesWithBonds == 0){
                    return;
                }

                //Get the device pointers
                const int* partLocalIndex2id_ptr = thrust::raw_pointer_cast(partLocalIndex2id_d.data());
                const typename BondType::BondParameters* bondList_ptr = thrust::raw_pointer_cast(bondList_d.data());

                const int* partLocalIndex2bondsStart_ptr      = thrust::raw_pointer_cast(partLocalIndex2bondsStart_d.data());
                const int* partLocalIndex2bonds_ptr           = thrust::raw_pointer_cast(partLocalIndex2bonds_d.data());
                const int* partLocalIndex2nBondsInvolved_ptr  = thrust::raw_pointer_cast(partLocalIndex2nBondsInvolved_d.data());

                if(comp.energy == true){

                    if constexpr (has_getEnergyTransverser<BondType>::value){

                        auto id2index = pd->getIdOrderedIndices(access::location::gpu);

                        int Nthreads=THREADS_PER_BLOCK;
                        int Nblocks=nParticlesWithBonds/Nthreads + ((nParticlesWithBonds%Nthreads)?1:0);

                        BondsInteractor_ns::transverseBondListThreadPerParticle<
                                            typename BondType::ComputationalData,
                                            typename BondType::BondParameters,
                                            typename BondType::EnergyTransverser>
                        <<<Nblocks,Nthreads,0,st>>>(nParticlesWithBonds,
                                                    pot->getComputationalData(comp,st),
                                                    id2index,
                                                    partLocalIndex2id_ptr,
                                                    bondList_ptr,
                                                    partLocalIndex2bondsStart_ptr,
                                                    partLocalIndex2bonds_ptr,
                                                    partLocalIndex2nBondsInvolved_ptr,
                                                    pot->getEnergyTransverser());
                        CudaCheckError();


                    } else {
                        if(!warningEnergy){
                            System::log<System::WARNING>("[BondsInteractor] (%s) Requested non-implemented transverser (energy)",
                                                         name.c_str());
                            warningEnergy = true;
                        }
                    }

                }

                if(comp.force == true){

                    if constexpr (has_getForceTransverser<BondType>::value){

                        auto id2index = pd->getIdOrderedIndices(access::location::gpu);

                        int Nthreads=THREADS_PER_BLOCK;
                        int Nblocks=nParticlesWithBonds/Nthreads + ((nParticlesWithBonds%Nthreads)?1:0);

                        BondsInteractor_ns::transverseBondListThreadPerParticle<
                                            typename BondType::ComputationalData,
                                            typename BondType::BondParameters,
                                            typename BondType::ForceTransverser>
                        <<<Nblocks,Nthreads,0,st>>>(nParticlesWithBonds,
                                                    pot->getComputationalData(comp,st),
                                                    id2index,
                                                    partLocalIndex2id_ptr,
                                                    bondList_ptr,
                                                    partLocalIndex2bondsStart_ptr,
                                                    partLocalIndex2bonds_ptr,
                                                    partLocalIndex2nBondsInvolved_ptr,
                                                    pot->getForceTransverser());
                        CudaCheckError();


                    } else {
                        if(!warningForce){
                            System::log<System::WARNING>("[BondsInteractor] (%s) Requested non-implemented transverser (force)",
                                                         name.c_str());
                            warningForce = true;
                        }
                    }

                }

                if(comp.lambdaDerivative == true){

                    if constexpr (has_getLambdaTransverser<BondType>::value){

                        auto id2index = pd->getIdOrderedIndices(access::location::gpu);

                        int Nthreads=THREADS_PER_BLOCK;
                        int Nblocks=nParticlesWithBonds/Nthreads + ((nParticlesWithBonds%Nthreads)?1:0);

                        BondsInteractor_ns::transverseBondListThreadPerParticle<
                                            typename BondType::ComputationalData,
                                            typename BondType::BondParameters,
                                            typename BondType::LambdaTransverser>
                        <<<Nblocks,Nthreads,0,st>>>(nParticlesWithBonds,
                                                    pot->getComputationalData(comp,st),
                                                    id2index,
                                                    partLocalIndex2id_ptr,
                                                    bondList_ptr,
                                                    partLocalIndex2bondsStart_ptr,
                                                    partLocalIndex2bonds_ptr,
                                                    partLocalIndex2nBondsInvolved_ptr,
                                                    pot->getLambdaTransverser());
                        CudaCheckError();

                    } else {
                        if(!warningLambdaDerivative){
                            System::log<System::WARNING>("[BondsInteractor] (%s) Requested non-implemented transverser (lambdaDerivative)",
                                                         name.c_str());
                            warningLambdaDerivative = true;
                        }
                    }
                }

                if(comp.stress == true){

                    if constexpr (has_getStressTransverser<BondType>::value){

                        auto id2index = pd->getIdOrderedIndices(access::location::gpu);

                        int Nthreads=THREADS_PER_BLOCK;
                        int Nblocks=nParticlesWithBonds/Nthreads + ((nParticlesWithBonds%Nthreads)?1:0);

                        BondsInteractor_ns::transverseBondListThreadPerParticle<
                                            typename BondType::ComputationalData,
                                            typename BondType::BondParameters,
                                            typename BondType::StressTransverser>
                        <<<Nblocks,Nthreads,0,st>>>(nParticlesWithBonds,
                                                    pot->getComputationalData(comp,st),
                                                    id2index,
                                                    partLocalIndex2id_ptr,
                                                    bondList_ptr,
                                                    partLocalIndex2bondsStart_ptr,
                                                    partLocalIndex2bonds_ptr,
                                                    partLocalIndex2nBondsInvolved_ptr,
                                                    pot->getStressTransverser());
                        CudaCheckError();


                    } else {
                        if(!warningStress){
                            System::log<System::WARNING>("[BondsInteractor] (%s) Requested non-implemented transverser (stress)",
                                                         name.c_str());
                        }
                    }

                }

                if(comp.magneticField == true){

                    if constexpr (has_getMagneticFieldTransverser<BondType>::value){

                        auto id2index = pd->getIdOrderedIndices(access::location::gpu);

                        int Nthreads=THREADS_PER_BLOCK;
                        int Nblocks=nParticlesWithBonds/Nthreads + ((nParticlesWithBonds%Nthreads)?1:0);

                        BondsInteractor_ns::transverseBondListThreadPerParticle<
                                            typename BondType::ComputationalData,
                                            typename BondType::BondParameters,
                                            typename BondType::MagneticFieldTransverser>
                        <<<Nblocks,Nthreads,0,st>>>(nParticlesWithBonds,
                                                    pot->getComputationalData(comp,st),
                                                    id2index,
                                                    partLocalIndex2id_ptr,
                                                    bondList_ptr,
                                                    partLocalIndex2bondsStart_ptr,
                                                    partLocalIndex2bonds_ptr,
                                                    partLocalIndex2nBondsInvolved_ptr,
                                                    pot->getMagneticFieldTransverser());
                        CudaCheckError();


                    } else {
                        if(!warningMagneticField){
                            System::log<System::WARNING>("[BondsInteractor] (%s) Requested non-implemented transverser (magneticField)",
                                                         name.c_str());
                        }
                    }

                }

                if(comp.hessian == true){

                    if constexpr (has_getHessianTransverser<BondType>::value){

                        auto id2index = pd->getIdOrderedIndices(access::location::gpu);

                        int Nthreads=THREADS_PER_BLOCK;
                        int Nblocks=nParticlesWithBonds/Nthreads + ((nParticlesWithBonds%Nthreads)?1:0);

                        BondsInteractor_ns::transverseBondListThreadPerParticle<
                                            typename BondType::ComputationalData,
                                            typename BondType::BondParameters,
                                            typename BondType::HessianTransverser>
                        <<<Nblocks,Nthreads,0,st>>>(nParticlesWithBonds,
                                                    pot->getComputationalData(comp,st),
                                                    id2index,
                                                    partLocalIndex2id_ptr,
                                                    bondList_ptr,
                                                    partLocalIndex2bondsStart_ptr,
                                                    partLocalIndex2bonds_ptr,
                                                    partLocalIndex2nBondsInvolved_ptr,
                                                    pot->getHessianTransverser());
                        CudaCheckError();


                    } else {
                        if(!warningHessian){
                            System::log<System::WARNING>("[BondsInteractor] (%s) Requested non-implemented transverser (hessian)",
                                                         name.c_str());
                        }
		    }
		}
		if(comp.pairwiseForce == true){
		  if constexpr (has_getPairwiseForceTransverser<BondType>::value){

		    auto id2index = pd->getIdOrderedIndices(access::location::gpu);

		    int Nthreads=THREADS_PER_BLOCK;
		    int Nblocks=nParticlesWithBonds/Nthreads + ((nParticlesWithBonds%Nthreads)?1:0);

		    BondsInteractor_ns::transverseBondListThreadPerParticle<
		      typename BondType::ComputationalData,
		      typename BondType::BondParameters,
		      typename BondType::PairwiseForceTransverser>
		      <<<Nblocks,Nthreads,0,st>>>(nParticlesWithBonds,
						  pot->getComputationalData(comp,st),
						  id2index,
						  partLocalIndex2id_ptr,
						  bondList_ptr,
						  partLocalIndex2bondsStart_ptr,
						  partLocalIndex2bonds_ptr,
						  partLocalIndex2nBondsInvolved_ptr,
						  pot->getPairwiseForceTransverser());
		    CudaCheckError();


		  } else {
		    if(!warningPairwiseForces){
		      System::log<System::WARNING>("[BondsInteractor] (%s) Requested non-implemented transverser (pairwiseForces)",
						   name.c_str());
		    }
		  }
		}
	    }
    };
}}}
