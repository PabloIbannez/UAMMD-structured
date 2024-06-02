#ifndef __PAIR_INTERACTION__
#define __PAIR_INTERACTION__

#include"Interactor/NeighbourList/CellList.cuh"

namespace uammd{
namespace structured{
namespace Interactor{

    namespace PairInteractor_ns{

    template<class ComputationalData, class NeighListData, class PotentialTransverser>
    __global__ void transverseNeighboursThreadPerParticle(int numberParticles,
                                                          ParticleGroup::MaskIterator groupMask,
                                                          ComputationalData computational,
                                                          NeighListData nlData,
                                                          PotentialTransverser potTransverser){
        int id = blockIdx.x*blockDim.x + threadIdx.x;

        if(id>=numberParticles){return;}

        const int nn   = nlData.numberNeighbours[id];

        if(nn>int(1)){

            int const * nl = nlData.neighbourList + nlData.neighbourStart[id];

            const int i_global = nl[0];

            typename PotentialTransverser::resultType quantityLocal = potTransverser.zero();

            for(int n = int(1); n<nn and groupMask[i_global]; n += int(1)) {

                const int j_global = nl[n*numberParticles];
                if(!groupMask[j_global]){
                    continue;
                }

                potTransverser.accumulate(quantityLocal,potTransverser.compute(i_global,
                                                                               j_global,
                                                                               computational));
            }

            potTransverser.set(i_global,quantityLocal);
        }
    }

    }

    template<class PotentialType, class NeighbourList, int THREADS_PER_BLOCK = 256>
    class PairInteractor: public Interactor{

        protected:

            std::shared_ptr<GlobalData>     gd;

            std::shared_ptr<PotentialType> pot;
            std::shared_ptr<NeighbourList> nl;

            std::string conditionInteractionName;

            ////////////////////////////////////
            //Warnings

            bool warningEnergy           = false;
            bool warningForce            = false;
            bool warningForceMagnetic    = false;
            bool warningLambdaDerivative = false;
            bool warningStress           = false;
            bool warningHessian          = false;
            bool warningMagneticField    = false;
            bool warningPairwiseForces   = false;

        public:

            PairInteractor(std::shared_ptr<GlobalData>           gd,
                           std::shared_ptr<ParticleGroup>        pg,
                           DataEntry& data,
                           std::shared_ptr<PotentialType> pot,
                           std::shared_ptr<NeighbourList> nl,
                           std::string name):Interactor(pg,"PairInteractor: \"" +name+"\""),
                                             gd(gd),
                                             pot(pot),
                                             nl(nl){

                conditionInteractionName = data.getParameter<std::string>("condition");
            }

            std::shared_ptr<PotentialType> getPotential(){
                return pot;
            }

            void sum(Computables comp,cudaStream_t st) override {

                nl->update(st);

                if(comp.energy == true){

                    if constexpr (has_getEnergyTransverser<PotentialType>::value){

                        auto nlData   = nl->getNeighbourList(conditionInteractionName);

                        int  numberParticles    = nlData.N; //numberParticles in the neig list
                        auto groupMask = pg->getGroupIndexMask(access::location::gpu);

                        int Nthreads=THREADS_PER_BLOCK;
                        int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);

                        PairInteractor_ns::transverseNeighboursThreadPerParticle<
                                             typename PotentialType::ComputationalData,
                                             typename NeighbourList::NeighbourListData,
                                             typename PotentialType::EnergyTransverser>
                        <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                     groupMask,
                                                     pot->getComputationalData(comp),
                                                     nlData,
                                                     pot->getEnergyTransverser());

                        CudaCheckError();


                    } else {
                        if(!warningEnergy){
                            System::log<System::WARNING>("[PairInteractor] (%s) Requested non-implemented transverser (energy)",
                                                         name.c_str());
                            warningEnergy = true;
                        }
                    }

                }

                if(comp.force == true){

		  if(comp.magneticField == false){

		    if constexpr (has_getForceTransverser<PotentialType>::value){

		      auto nlData   = nl->getNeighbourList(conditionInteractionName);

		      int  numberParticles    = nlData.N; //numberParticles in the neig list
                            auto groupMask = pg->getGroupIndexMask(access::location::gpu);

                            int Nthreads=THREADS_PER_BLOCK;
                            int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);

                            PairInteractor_ns::transverseNeighboursThreadPerParticle<
                                                 typename PotentialType::ComputationalData,
                                                 typename NeighbourList::NeighbourListData,
                                                 typename PotentialType::ForceTransverser>
                            <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                         groupMask,
                                                         pot->getComputationalData(comp),
                                                         nlData,
                                                         pot->getForceTransverser());

                            CudaCheckError();


                        } else {
                            if(!warningForce){
                                System::log<System::WARNING>("[PairInteractor] (%s) Requested non-implemented transverser (force)",
                                                             name.c_str());
                                warningForce = true;
                            }
                        }

                    } else {

		      if constexpr (has_getForceTorqueMagneticFieldTransverser<PotentialType>::value){

			auto nlData   = nl->getNeighbourList(conditionInteractionName);

			int  numberParticles    = nlData.N; //numberParticles in the neig list
                            auto groupMask = pg->getGroupIndexMask(access::location::gpu);

                            int Nthreads=THREADS_PER_BLOCK;
                            int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);

                            PairInteractor_ns::transverseNeighboursThreadPerParticle<
                                                 typename PotentialType::ComputationalData,
                                                 typename NeighbourList::NeighbourListData,
                                                 typename PotentialType::ForceTorqueMagneticFieldTransverser>
                            <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                         groupMask,
                                                         pot->getComputationalData(comp),
                                                         nlData,
                                                         pot->getForceTorqueMagneticFieldTransverser());

                            CudaCheckError();


		      } else {
			if(!warningForceMagnetic){
			  System::log<System::WARNING>("[PairInteractor] (%s) Requested non-implemented transverser (forceTorqueMagneticField)",
						       name.c_str());
			  warningForceMagnetic = true;
			}
		      }
                    }
                }

		if(comp.magneticField == true and comp.force == false){
		  if constexpr (has_getMagneticFieldTransverser<PotentialType>::value){

			auto nlData   = nl->getNeighbourList(conditionInteractionName);

			int  numberParticles    = nlData.N; //numberParticles in the neig list
			auto groupMask = pg->getGroupIndexMask(access::location::gpu);

			int Nthreads=THREADS_PER_BLOCK;
			int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);

			PairInteractor_ns::transverseNeighboursThreadPerParticle<
			  typename PotentialType::ComputationalData,
			  typename NeighbourList::NeighbourListData,
			  typename PotentialType::MagneticFieldTransverser>
			  <<<Nblocks,Nthreads,0, st>>>(numberParticles,
						       groupMask,
						       pot->getComputationalData(comp),
						       nlData,
						       pot->getMagneticFieldTransverser());

			CudaCheckError();


		  } else {
		    if(!warningMagneticField){
		      System::log<System::WARNING>("[PairInteractor] (%s) Requested non-implemented transverser (MagneticField)",
						   name.c_str());
		      warningForceMagnetic = true;
		    }
		  }
		}

                if(comp.lambdaDerivative == true){

                    if constexpr (has_getLambdaTransverser<PotentialType>::value){

                        auto nlData   = nl->getNeighbourList(conditionInteractionName);

                        int  numberParticles    = nlData.N; //numberParticles in the neig list
                        auto groupMask = pg->getGroupIndexMask(access::location::gpu);

                        int Nthreads=THREADS_PER_BLOCK;
                        int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);

                        PairInteractor_ns::transverseNeighboursThreadPerParticle<
                                             typename PotentialType::ComputationalData,
                                             typename NeighbourList::NeighbourListData,
                                             typename PotentialType::LambdaTransverser>
                        <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                     groupMask,
                                                     pot->getComputationalData(comp),
                                                     nlData,
                                                     pot->getLambdaTransverser());

                        CudaCheckError();

                    } else {
                        if(!warningLambdaDerivative){
                            System::log<System::WARNING>("[PairInteractor] (%s) Requested non-implemented transverser (lambdaDerivative)",
                                                         name.c_str());
                            warningLambdaDerivative = true;
                        }
                    }

                }

                if(comp.stress == true){

                    if constexpr (has_getStressTransverser<PotentialType>::value){

                        auto nlData   = nl->getNeighbourList(conditionInteractionName);

                        int  numberParticles    = nlData.N; //numberParticles in the neig list
                        auto groupMask = pg->getGroupIndexMask(access::location::gpu);

                        int Nthreads=THREADS_PER_BLOCK;
                        int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);

                        PairInteractor_ns::transverseNeighboursThreadPerParticle<
                                             typename PotentialType::ComputationalData,
                                             typename NeighbourList::NeighbourListData,
                                             typename PotentialType::StressTransverser>
                        <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                     groupMask,
                                                     pot->getComputationalData(comp),
                                                     nlData,
                                                     pot->getStressTransverser());

                        CudaCheckError();


                    } else {
                        if(!warningStress){
                            System::log<System::WARNING>("[PairInteractor] (%s) Requested non-implemented transverser (stress)",
                                                         name.c_str());
                            warningStress = true;
                        }
                    }
                }

		if(comp.hessian == true){

                    if constexpr (has_getHessianTransverser<PotentialType>::value){

                        auto nlData   = nl->getNeighbourList(conditionInteractionName);

                        int  numberParticles    = nlData.N; //numberParticles in the neig list
                        auto groupMask = pg->getGroupIndexMask(access::location::gpu);

                        int Nthreads=THREADS_PER_BLOCK;
                        int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);

                        PairInteractor_ns::transverseNeighboursThreadPerParticle<
                                             typename PotentialType::ComputationalData,
                                             typename NeighbourList::NeighbourListData,
                                             typename PotentialType::HessianTransverser>
                        <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                     groupMask,
                                                     pot->getComputationalData(comp),
                                                     nlData,
                                                     pot->getHessianTransverser());

                        CudaCheckError();


                    } else {
                        if(!warningHessian){
                            System::log<System::WARNING>("[PairInteractor] (%s) Requested non-implemented transverser (Hessian)",
                                                         name.c_str());
                            warningHessian = true;
                        }
                    }
                }

		if(comp.pairwiseForce == true){

                    if constexpr (has_getPairwiseForceTransverser<PotentialType>::value){

                        auto nlData   = nl->getNeighbourList(conditionInteractionName);

                        int  numberParticles    = nlData.N; //numberParticles in the neig list
                        auto groupMask = pg->getGroupIndexMask(access::location::gpu);

                        int Nthreads=THREADS_PER_BLOCK;
                        int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);

                        PairInteractor_ns::transverseNeighboursThreadPerParticle<
                                             typename PotentialType::ComputationalData,
                                             typename NeighbourList::NeighbourListData,
                                             typename PotentialType::PairwiseForceTransverser>
                        <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                     groupMask,
                                                     pot->getComputationalData(comp),
                                                     nlData,
                                                     pot->getPairwiseForceTransverser());

                        CudaCheckError();


                    } else {
                        if(!warningPairwiseForces){
                            System::log<System::WARNING>("[PairInteractor] (%s) Requested non-implemented transverser (PairwiseForce)",
                                                         name.c_str());
                            warningPairwiseForces = true;
                        }
                    }
                }

                if(comp.magneticField == true){

                    if constexpr (has_getMagneticFieldTransverser<PotentialType>::value){

                        auto nlData   = nl->getNeighbourList(conditionInteractionName);

                        int  numberParticles    = nlData.N; //numberParticles in the neig list
                        auto groupMask = pg->getGroupIndexMask(access::location::gpu);

                        int Nthreads=THREADS_PER_BLOCK;
                        int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);

                        PairInteractor_ns::transverseNeighboursThreadPerParticle<
                                             typename PotentialType::ComputationalData,
                                             typename NeighbourList::NeighbourListData,
                                             typename PotentialType::MagneticFieldTransverser>
                        <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                     groupMask,
                                                     pot->getComputationalData(comp),
                                                     nlData,
                                                     pot->getMagneticFieldTransverser());

                        CudaCheckError();


                    } else {
                        if(!warningMagneticField){
                            System::log<System::WARNING>("[PairInteractor] (%s) Requested non-implemented transverser (magneticField)",
                                                         name.c_str());
                            warningMagneticField = true;
                        }
                    }
                }

            }
    };

}}}

#endif
