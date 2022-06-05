#ifndef __RAVIKUMAR_UNBOUND__
#define __RAVIKUMAR_UNBOUND__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace UnBound{
    
    namespace RavikumarPotential_ns{
    namespace ModelType{

        struct A{

            static constexpr real lambda = 0.4;
            static constexpr real epsilon_0 = -1.3;
            
            static constexpr real gamma = 0.625;
        };
        
    }}

    namespace RavikumarPotential_ns{
        namespace NonPolar{
            
            namespace epsilonPositive{

                inline __device__ real3 force(const real3& rij, const real& r2,
                                              const real& d,
                                              const real& epsilon,const real& sigma){

                    const real r       = sqrt(r2);
                    const real invr2   = real(1.0)/r2;
                    const real sinvr2  = sigma*sigma*invr2;
                    const real sinvr6  = sinvr2*sinvr2*sinvr2;
                    const real sinvr12 = sinvr6*sinvr6;
                    
                    const real exponent = (r-sigma)/d;

                    const real d2 = d*d;

                    const real fmod = real(10.0)*epsilon*sinvr12*(exp(-exponent*exponent)*((real(6.0)*d2+r*(r-sigma))/d2)-real(6.0))*invr2;
                    
                    return fmod*rij;
                }
                
                inline __device__ tensor3 virial(const real3& rij, const real& r2,
                                                 const real& d,
                                                 const real& epsilon,const real& sigma){
                    return computeVirial(rij,force(rij,r2,d,epsilon,sigma));
                }
                
                inline __device__ real energy(const real3& rij, const real& r2,
                                              const real& d,
                                              const real& epsilon,const real& sigma){
                    
                    const real invr2   = real(1.0)/r2;
                    const real sinvr2  = sigma*sigma*invr2;
                    const real sinvr6  = sinvr2*sinvr2*sinvr2;
                    const real sinvr12 = sinvr6*sinvr6;

                    real exponent = (sqrt(r2)-sigma)/d;
                         exponent = -exponent*exponent;

                    real e = epsilon*(real(5.0)*sinvr12*(real(1.0)-exp(exponent)));
                    
                    return e/real(2.0);
                }
            }

            //Force
            inline __device__ real3 force(const real3& rij, const real& r2,const real& d, 
                                          const real& epsilon,const real& sigma,const real& zeroEnergy){
                
                if(epsilon >= real(0.0) ){
                    if(epsilon == real(0.0) ){
                        return epsilonPositive::force(rij,r2,d,zeroEnergy,sigma);
                    } else {
                        return epsilonPositive::force(rij,r2,d,epsilon,sigma);
                    }
                } else {
                        return CommonPotentials::LennardJones::Type3::force(rij,r2,abs(epsilon),sigma);
                }
            }
            
            //Virial
            inline __device__ tensor3 virial(const real3& rij, const real& r2,const real& d,
                                             const real& epsilon,const real& sigma,const real& zeroEnergy){
                return tensor3(0);
            }
            
            //Energy
            inline __device__ real energy(const real3& rij, const real& r2,const real& d, 
                                          const real& epsilon,const real& sigma,const real& zeroEnergy){
                
                if(epsilon >= real(0.0)){
                    if(epsilon == real(0.0) ){
                        return epsilonPositive::energy(rij,r2,d,zeroEnergy,sigma);
                    } else {
                        return epsilonPositive::energy(rij,r2,d,epsilon,sigma);
                    }
                } else {
                        return CommonPotentials::LennardJones::Type3::energy(rij,r2,abs(epsilon),sigma);
                }
            }

        }
    }

    template<typename Units_>
    struct RavikumarPotential: public ParameterUpdatable{
        
        using InteractionParameters = typename CommonParameters::StatisticalPotential::InteractionParameters;
        using ParameterPairsHandler = typename structured::PairParameterHandler<InteractionParameters>;
        
        std::shared_ptr<ParticleData> pd;
        
        std::shared_ptr<ParameterPairsHandler> paramPairsHandler;
        
        Box box;

        real dielectricConstant;
        real debyeLenght;
        
        real cutOffDstNP;
        real cutOffDstDH;

        real d;
        real gamma;

        real zeroEnergy;
        
        struct Parameters{
            
            std::shared_ptr<ParameterPairsHandler> paramPairsHandler;
            
            real dielectricConstant;
            real debyeLenght;
            
            real cutOffDstNP;
            real cutOffDstDH;
        
            real d;
            real gamma;
        
            real zeroEnergy = real(0.01);
        
        };

        RavikumarPotential(std::shared_ptr<ParticleData> pd,
                           Parameters par):pd(pd),
                                           dielectricConstant(par.dielectricConstant),
                                           debyeLenght(par.debyeLenght),
                                           cutOffDstNP(par.cutOffDstNP),
                                           cutOffDstDH(par.cutOffDstDH),
                                           d(par.d),
                                           gamma(par.gamma),
                                           zeroEnergy(par.zeroEnergy),
                                           paramPairsHandler(par.paramPairsHandler){}

        ~RavikumarPotential(){}

        void setZeroEnergy(real newZeroEnergy){
            zeroEnergy = newZeroEnergy;
        }

        struct forceTransverser{

            real4* force;
            
            real* radius;
            real* charge;

            typename ParameterPairsHandler::PairIterator paramPairIterator;

            Box box;
            
            real dielectricConstant;
            real debyeLenght;

            real cutOffDstNP2;
            real cutOffDstDH2;

            real d;
            real gamma;

            real zeroEnergy;
            
            forceTransverser(real4* force,
                             real* radius,
                             real* charge,
                             typename ParameterPairsHandler::PairIterator paramPairIterator,
                             Box  box,
                             real dielectricConstant,
                             real debyeLenght,
                             real cutOffDstNP2,
                             real cutOffDstDH2,
                             real d,
                             real gamma,
                             real zeroEnergy):force(force),
                                              radius(radius),
                                              charge(charge),
                                              paramPairIterator(paramPairIterator),
                                              box(box),
                                              dielectricConstant(dielectricConstant),
                                              debyeLenght(debyeLenght),
                                              cutOffDstNP2(cutOffDstNP2),
                                              cutOffDstDH2(cutOffDstDH2),
                                              d(d),
                                              gamma(gamma),
                                              zeroEnergy(zeroEnergy){}

            using resultType=real4;

            inline __device__ resultType zero(){return make_real4(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                real3 f = make_real3(real(0.0));

                if(r2>0 and r2<=cutOffDstNP2){

                    const real sigma = gamma*(radius[index_i]+radius[index_j]);
                    const real eps   = paramPairIterator(int(posi.w),int(posj.w)).epsilon;

                    f+=RavikumarPotential_ns::NonPolar::force(rij,r2,d,eps,sigma,zeroEnergy);
                } 
                
                const real chgProduct = charge[index_i]*charge[index_j];
                if(r2>0 and r2<=cutOffDstDH2 and chgProduct != real(0.0)){
                    
                    f+=CommonPotentials::DebyeHuckel::DebyeHuckel::force<Units_>(rij,r2,
                                                                                 chgProduct,
                                                                                 dielectricConstant,debyeLenght);
                }

                return make_real4(f,real(0.0));
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){force[index_i]+=quantity;}

        };

        forceTransverser getForceTransverser(){
            
            real4* force = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();     
            
            real* radius = this->pd->getRadius(access::location::gpu, access::mode::readwrite).raw();     
            real* charge = this->pd->getCharge(access::location::gpu, access::mode::readwrite).raw();     
                                                         
            return forceTransverser(force,
                                    radius,charge,
                                    paramPairsHandler->getPairIterator(),
                                    box,
                                    dielectricConstant,
                                    debyeLenght,
                                    cutOffDstNP*cutOffDstNP,
                                    cutOffDstDH*cutOffDstDH,
                                    d,
                                    gamma,
                                    zeroEnergy);
        }
        
        struct virialTransverser{

            tensor3* virial;
            
            real* radius;
            real* charge;

            typename ParameterPairsHandler::PairIterator paramPairIterator;

            Box box;
            
            real dielectricConstant;
            real debyeLenght;

            real cutOffDstNP2;
            real cutOffDstDH2;
            
            real d;
            real gamma;

            real zeroEnergy;
            
            virialTransverser(tensor3* virial,
                              real* radius,
                              real* charge,
                              typename ParameterPairsHandler::PairIterator paramPairIterator,
                              Box  box,
                              real dielectricConstant,
                              real debyeLenght,
                              real cutOffDstNP2,
                              real cutOffDstDH2,
                              real d,
                              real gamma,
                              real zeroEnergy):virial(virial),
                                               radius(radius),
                                               charge(charge),
                                               paramPairIterator(paramPairIterator),
                                               box(box),
                                               dielectricConstant(dielectricConstant),
                                               debyeLenght(debyeLenght),
                                               cutOffDstNP2(cutOffDstNP2),
                                               cutOffDstDH2(cutOffDstDH2),
                                               d(d),
                                               gamma(gamma),
                                               zeroEnergy(zeroEnergy){}

            using resultType=tensor3;

            inline __device__ resultType zero(){return tensor3(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                tensor3 v = tensor3(real(0.0));

                return v;
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){virial[index_i]+=quantity;}

        };

        virialTransverser getVirialTransverser(){
            
            tensor3* virial = this->pd->getVirial(access::location::gpu, access::mode::readwrite).raw();     
            
            real* radius = this->pd->getRadius(access::location::gpu, access::mode::readwrite).raw();     
            real* charge = this->pd->getCharge(access::location::gpu, access::mode::readwrite).raw();     
                                                         
            return virialTransverser(virial,
                                     radius,charge,
                                     paramPairsHandler->getPairIterator(),
                                     box,
                                     dielectricConstant,
                                     debyeLenght,
                                     cutOffDstNP*cutOffDstNP,
                                     cutOffDstDH*cutOffDstDH,
                                     d,gamma,
                                     zeroEnergy);
        }
        
        struct energyTransverser{

            real* energy;
            
            real* radius;
            real* charge;

            typename ParameterPairsHandler::PairIterator paramPairIterator;

            Box box;
            
            real dielectricConstant;
            real debyeLenght;

            real cutOffDstNP2;
            real cutOffDstDH2;

            real d;
            real gamma;

            real zeroEnergy;
            
            energyTransverser(real* energy,
                              real* radius,
                              real* charge,
                              typename ParameterPairsHandler::PairIterator paramPairIterator,
                              Box  box,
                              real dielectricConstant,
                              real debyeLenght,
                              real cutOffDstNP2,
                              real cutOffDstDH2,
                              real d,
                              real gamma,
                              real zeroEnergy):energy(energy),
                                               radius(radius),
                                               charge(charge),
                                               paramPairIterator(paramPairIterator),
                                               box(box),
                                               dielectricConstant(dielectricConstant),
                                               debyeLenght(debyeLenght),
                                               cutOffDstNP2(cutOffDstNP2),
                                               cutOffDstDH2(cutOffDstDH2),
                                               d(d),gamma(gamma),
                                               zeroEnergy(zeroEnergy){}

            using resultType=real;

            inline __device__ resultType zero(){return real(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                real e = real(0.0);

                if(r2>0 and r2<=cutOffDstNP2){

                    const real sigma = gamma*(radius[index_i]+radius[index_j]);
                    const real eps   = paramPairIterator(int(posi.w),int(posj.w)).epsilon;

                    e+=RavikumarPotential_ns::NonPolar::energy(rij,r2,d,eps,sigma,zeroEnergy);
                } 
                
                const real chgProduct = charge[index_i]*charge[index_j];
                if(r2>0 and r2<=cutOffDstDH2 and chgProduct != real(0.0)){
                    
                    e+=CommonPotentials::DebyeHuckel::DebyeHuckel::energy<Units_>(rij,r2,
                                                                                  chgProduct,
                                                                                  dielectricConstant,debyeLenght);
                }

                return e;
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){energy[index_i]+=quantity;}

        };

        energyTransverser getEnergyTransverser(){
            
            real* energy = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();     
            
            real* radius = this->pd->getRadius(access::location::gpu, access::mode::readwrite).raw();     
            real* charge = this->pd->getCharge(access::location::gpu, access::mode::readwrite).raw();     
                                                         
            return energyTransverser(energy,
                                     radius,charge,
                                     paramPairsHandler->getPairIterator(),
                                     box,
                                     dielectricConstant,
                                     debyeLenght,
                                     cutOffDstNP*cutOffDstNP,
                                     cutOffDstDH*cutOffDstDH,
                                     d,gamma,
                                     zeroEnergy);
        }
        
        void updateBox(Box newBox) override {
            box=newBox;
        }
    };

}}}}

#endif
