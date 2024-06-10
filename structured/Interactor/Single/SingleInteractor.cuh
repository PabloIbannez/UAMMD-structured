#ifndef __SINGLE_INTERACTION__
#define __SINGLE_INTERACTION__

namespace uammd{
namespace structured{
namespace Interactor{

    namespace SingleInteractor_ns{

    template<class ComputationalData, class PotentialTransverser>
    __global__ void transverseThreadPerParticle(int numberParticles,
                                                ParticleGroup::IndexIterator groupIndex,
                                                ComputationalData computational,
                                                PotentialTransverser potTransverser){

        int index = blockIdx.x*blockDim.x + threadIdx.x;

        if(index>=numberParticles){return;}

        int indexGlobal = groupIndex[index];

        typename PotentialTransverser::resultType quantityLocal = potTransverser.zero();

        potTransverser.accumulate(quantityLocal,potTransverser.compute(indexGlobal,computational));
        potTransverser.set(indexGlobal,quantityLocal);

    }

    }

    template<class PotentialType, int THREADS_PER_BLOCK = 256>
    class SingleInteractor: public Interactor{

        protected:

            std::shared_ptr<GlobalData>     gd;

            std::shared_ptr<PotentialType>  pot;

            ////////////////////////////////////
            //Warnings

            bool warningEnergy           = false;
            bool warningForce            = false;
            bool warningForceMagnetic    = false;
            bool warningMagneticField    = false;
            bool warningHessian          = false;
            bool warningLambdaDerivative = false;

        public:

            SingleInteractor(std::shared_ptr<GlobalData>           gd,
                             std::shared_ptr<ParticleGroup>        pg,
                             DataEntry& data,
                             std::shared_ptr<PotentialType> pot,
                             std::string name):Interactor(pg,"SingleInteractor: \"" +name+"\""),
                                               gd(gd),
                                               pot(pot){
            }

            std::shared_ptr<PotentialType> getPotential(){
                return pot;
            }

            void sum(Computables comp,cudaStream_t st) override {

                if(comp.energy == true){

                    if constexpr (has_getEnergyTransverser<PotentialType>::value){

                        int numberParticles     = pg->getNumberParticles();
                        auto groupIndexIterator = pg->getIndexIterator(access::location::gpu);

                        int Nthreads=THREADS_PER_BLOCK;
                        int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);

                        SingleInteractor_ns::transverseThreadPerParticle
                        <typename PotentialType::ComputationalData,
                         typename PotentialType::EnergyTransverser>
                        <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                     groupIndexIterator,
                                                     pot->getComputationalData(comp,st),
                                                     pot->getEnergyTransverser());
                        CudaCheckError();
                    } else {
                        if(!warningEnergy){
                            System::log<System::WARNING>("[SingleInteractor] (%s) Requested non-implemented transverser (energy)",
                                                         name.c_str());
                            warningEnergy = true;
                        }
                    }
                }

                if(comp.force == true){

                    if(comp.magneticField == false){

                        if constexpr (has_getForceTransverser<PotentialType>::value){

                            int numberParticles     = pg->getNumberParticles();
                            auto groupIndexIterator = pg->getIndexIterator(access::location::gpu);

                            int Nthreads=THREADS_PER_BLOCK;
                            int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);

                            SingleInteractor_ns::transverseThreadPerParticle
                            <typename PotentialType::ComputationalData,
                             typename PotentialType::ForceTransverser>
                            <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                         groupIndexIterator,
                                                         pot->getComputationalData(comp,st),
                                                         pot->getForceTransverser());
                            CudaCheckError();

                        } else {
                            if(!warningEnergy){
                                System::log<System::WARNING>("[SingleInteractor] (%s) Requested non-implemented transverser (force)",
                                                             name.c_str());
                                warningForce = true;
                            }
                        }

                    } else {

		      if constexpr (has_getForceTorqueMagneticFieldTransverser<PotentialType>::value){

                            int numberParticles     = pg->getNumberParticles();
                            auto groupIndexIterator = pg->getIndexIterator(access::location::gpu);

                            int Nthreads=THREADS_PER_BLOCK;
                            int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);

                            SingleInteractor_ns::transverseThreadPerParticle
                            <typename PotentialType::ComputationalData,
                             typename PotentialType::ForceTorqueMagneticFieldTransverser>
                            <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                         groupIndexIterator,
                                                         pot->getComputationalData(comp,st),
                                                         pot->getForceTorqueMagneticFieldTransverser());
                            CudaCheckError();

                        } else {
                            if(!warningForceMagnetic){
                                System::log<System::WARNING>("[SingleInteractor] (%s) Requested non-implemented transverser (forceMagneticField)",
                                                             name.c_str());

                                warningForceMagnetic = true;
                            }
                        }
                    }
                }

                if(comp.magneticField == true){

                    if constexpr (has_getMagneticFieldTransverser<PotentialType>::value){

                        int numberParticles     = pg->getNumberParticles();
                        auto groupIndexIterator = pg->getIndexIterator(access::location::gpu);

                        int Nthreads=THREADS_PER_BLOCK;
                        int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);

                        SingleInteractor_ns::transverseThreadPerParticle
                        <typename PotentialType::ComputationalData,
                         typename PotentialType::MagneticFieldTransverser>
                        <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                     groupIndexIterator,
                                                     pot->getComputationalData(comp,st),
                                                     pot->getMagneticFieldTransverser());
                        CudaCheckError();

                    } else {
                        if(!warningMagneticField){
                            System::log<System::WARNING>("[SingleInteractor] (%s) Requested non-implemented transverser (magneticField)",
                                                         name.c_str());
                            warningMagneticField = true;
                        }
                    }
                }

		if(comp.hessian == true){

		  if constexpr (has_getHessianTransverser<PotentialType>::value){

                        int numberParticles     = pg->getNumberParticles();
                        auto groupIndexIterator = pg->getIndexIterator(access::location::gpu);

                        int Nthreads=THREADS_PER_BLOCK;
                        int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);

                        SingleInteractor_ns::transverseThreadPerParticle
                        <typename PotentialType::ComputationalData,
                         typename PotentialType::HessianTransverser>
                        <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                     groupIndexIterator,
                                                     pot->getComputationalData(comp,st),
                                                     pot->getHessianTransverser());
                        CudaCheckError();

                    } else {
                        if(!warningHessian){
                            System::log<System::WARNING>("[SingleInteractor] (%s) Requested non-implemented transverser (hessian)",
                                                         name.c_str());
                            warningHessian = true;
                        }
                    }
                }

                if(comp.lambdaDerivative == true){

                    if constexpr (has_getLambdaTransverser<PotentialType>::value){

                        int numberParticles     = pg->getNumberParticles();
                        auto groupIndexIterator = pg->getIndexIterator(access::location::gpu);

                        int Nthreads=THREADS_PER_BLOCK;
                        int Nblocks=numberParticles/Nthreads + ((numberParticles%Nthreads)?1:0);

                        SingleInteractor_ns::transverseThreadPerParticle
                        <typename PotentialType::ComputationalData,
                         typename PotentialType::LambdaTransverser>
                        <<<Nblocks,Nthreads,0, st>>>(numberParticles,
                                                     groupIndexIterator,
                                                     pot->getComputationalData(comp,st),
                                                     pot->getLambdaTransverser());
                        CudaCheckError();

                    } else {
                        if(!warningLambdaDerivative){
                            System::log<System::WARNING>("[SingleInteractor] (%s) Requested non-implemented transverser (lambdaDerivative)",
                                                         name.c_str());
                            warningLambdaDerivative = true;
                        }
                    }

                }
            }
    };

}}}

#endif
