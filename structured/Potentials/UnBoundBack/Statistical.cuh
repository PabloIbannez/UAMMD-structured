#ifndef __STATISTICAL_UNBOUND__
#define __STATISTICAL_UNBOUND__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace UnBound{
    
    struct StatisticalPotential: public ParameterUpdatable{
        
        using InteractionParameters = typename CommonParameters::StatisticalPotential::InteractionParameters;
        using ParameterPairsHandler = typename structured::PairParameterHandler<InteractionParameters>;
        
        std::shared_ptr<ParticleData> pd;
        
        std::shared_ptr<ParameterPairsHandler> paramPairsHandler;
        
        Box box;

        real cutOffDst;
            
        real n;
        real r0;

        real zeroEnergy;
        
        struct Parameters{
            
            std::shared_ptr<ParameterPairsHandler> paramPairsHandler;
            
            real cutOffDst;

            real n;
            real r0;
        
            real zeroEnergy = real(0.01);
        };

        StatisticalPotential(std::shared_ptr<ParticleData> pd,
                             Parameters par):pd(pd),
                                             cutOffDst(par.cutOffDst),
                                             n(par.n),r0(par.r0),
                                             zeroEnergy(par.zeroEnergy),
                                             paramPairsHandler(par.paramPairsHandler){}

        ~StatisticalPotential(){}

        void setZeroEnergy(real newZeroEnergy){
            zeroEnergy = newZeroEnergy;
        }

        struct forceTransverser{

            real4* force;
            
            real* radius;

            typename ParameterPairsHandler::PairIterator paramPairIterator;

            Box box;
            
            real cutOffDst2;

            real n;
            real r0;

            real zeroEnergy;
            
            forceTransverser(real4* force,
                             real* radius,
                             typename ParameterPairsHandler::PairIterator paramPairIterator,
                             Box  box,
                             real cutOffDst2,
                             real n,real r0,
                             real zeroEnergy):force(force),
                                              radius(radius),
                                              paramPairIterator(paramPairIterator),
                                              box(box),
                                              cutOffDst2(cutOffDst2),
                                              n(n),r0(r0),
                                              zeroEnergy(zeroEnergy){}

            using resultType=real4;

            inline __device__ resultType zero(){return make_real4(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                real3 f = make_real3(real(0.0));

                if(r2>0 and r2<=cutOffDst2){

                    const real sigma = radius[index_i]+radius[index_j];
                    const real eps   = paramPairIterator(int(posi.w),int(posj.w)).epsilon;
                    
                    f+=CommonPotentials::Contact::force(rij,r2,eps,sigma,n,r0,zeroEnergy);                
                } 

                return make_real4(f,real(0.0));
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){force[index_i]+=quantity;}

        };

        forceTransverser getForceTransverser(){
            
            real4* force = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();     
            
            real* radius = this->pd->getRadius(access::location::gpu, access::mode::readwrite).raw();     
                                                         
            return forceTransverser(force,
                                    radius,
                                    paramPairsHandler->getPairIterator(),
                                    box,
                                    cutOffDst*cutOffDst,
                                    n,r0,
                                    zeroEnergy);
        }
        
        struct virialTransverser{

            tensor3* virial;
            
            real* radius;

            typename ParameterPairsHandler::PairIterator paramPairIterator;

            Box box;
            
            real cutOffDst2;

            real n;
            real r0;

            real zeroEnergy;
            
            virialTransverser(tensor3* virial,
                              real* radius,
                              typename ParameterPairsHandler::PairIterator paramPairIterator,
                              Box  box,
                              real cutOffDst2,
                              real n,real r0,
                              real zeroEnergy):virial(virial),
                                               radius(radius),
                                               paramPairIterator(paramPairIterator),
                                               box(box),
                                               cutOffDst2(cutOffDst2),
                                               n(n),r0(r0),
                                               zeroEnergy(zeroEnergy){}

            using resultType=tensor3;

            inline __device__ resultType zero(){return tensor3(real(0.0));}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                tensor3 v = tensor3(real(0.0));

                if(r2>0 and r2<=cutOffDst2){

                    const real sigma = radius[index_i]+radius[index_j];
                    const real eps   = paramPairIterator(int(posi.w),int(posj.w)).epsilon;
                    
                    v+=CommonPotentials::Contact::virial(rij,r2,eps,sigma,n,r0,zeroEnergy);                
                } 

                return v;
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){virial[index_i]+=quantity;}

        };

        virialTransverser getVirialTransverser(){
            
            tensor3* virial = this->pd->getVirial(access::location::gpu, access::mode::readwrite).raw();     
            real*    radius = this->pd->getRadius(access::location::gpu, access::mode::readwrite).raw();     
                                                         
            return virialTransverser(virial,
                                     radius,
                                     paramPairsHandler->getPairIterator(),
                                     box,
                                     cutOffDst*cutOffDst,
                                     n,r0,
                                     zeroEnergy);
        }
        
        struct energyTransverser{

            real* energy;
            
            real* radius;

            typename ParameterPairsHandler::PairIterator paramPairIterator;

            Box box;
            
            real cutOffDst2;

            real n;
            real r0;

            real zeroEnergy;
            
            energyTransverser(real* energy,
                              real* radius,
                              typename ParameterPairsHandler::PairIterator paramPairIterator,
                              Box  box,
                              real cutOffDst2,
                              real n,real r0,
                              real zeroEnergy):energy(energy),
                                               radius(radius),
                                               paramPairIterator(paramPairIterator),
                                               box(box),
                                               cutOffDst2(cutOffDst2),
                                               n(n),r0(r0),
                                               zeroEnergy(zeroEnergy){}

            using resultType=real;

            inline __device__ resultType zero(){return real(0.0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                real e = real(0.0);

                if(r2>0 and r2<=cutOffDst2){

                    const real sigma = radius[index_i]+radius[index_j];
                    const real eps   = paramPairIterator(int(posi.w),int(posj.w)).epsilon;
                    
                    e+=CommonPotentials::Contact::energy(rij,r2,eps,sigma,n,r0,zeroEnergy);                
                } 

                return e;
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){energy[index_i]+=quantity;}

        };

        energyTransverser getEnergyTransverser(){
            
            real* energy = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();     
            
            real* radius = this->pd->getRadius(access::location::gpu, access::mode::readwrite).raw();     
                                                         
            return energyTransverser(energy,
                                     radius,
                                     paramPairsHandler->getPairIterator(),
                                     box,
                                     cutOffDst*cutOffDst,
                                     n,r0,
                                     zeroEnergy);
        }
        
        void updateBox(Box newBox) override {
            box=newBox;
        }
    };

}}}}

#endif
