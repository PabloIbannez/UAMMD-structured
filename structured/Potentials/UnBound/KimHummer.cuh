#ifndef __KIM_HUMMER_UNBOUND__
#define __KIM_HUMMER_UNBOUND__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace UnBound{
    
    namespace KimHummerPotential_ns{
    namespace ModelType{

        struct A{

            static constexpr bool SASArequired = false;
            
            static constexpr real lambda = 0.159;
            static constexpr real epsilon_0 = -2.27;
            
            __host__ __device__ inline static real SASAweight(real SASAratio){
                return real(1.0);
            }
        };
        
        struct B{
            
            static constexpr bool SASArequired = true;
            
            static constexpr real lambda = 0.186;
            static constexpr real epsilon_0 = -1.95;
            
            __host__ __device__ inline static real SASAweight(real SASAratio){
                real w = tanhf(real(10)*tanf(SASAratio*real(M_PI_2)));
                return (w<real(0.0))?real(1.0):w;
            }
        };
        
        struct C{
            
            static constexpr bool SASArequired = true;
            
            static constexpr real lambda = 0.192;
            static constexpr real epsilon_0 = -1.85;
            
            __host__ __device__ inline static real SASAweight(real SASAratio){
                real w = tanhf(real(5)*tanf(SASAratio*real(M_PI_2)));
                return (w<real(0.0))?real(1.0):w;
            }
        };
        
        struct D{
            
            static constexpr bool SASArequired = true;
            
            static constexpr real lambda = 0.228;
            static constexpr real epsilon_0 = -1.67;
            
            __host__ __device__ inline static real SASAweight(real SASAratio){
                real w = tanhf(real(2)*tanf(SASAratio*real(M_PI_2)));
                return (w<real(0.0))?real(1.0):w;
            }
        };
        
        struct E{
            
            static constexpr bool SASArequired = true;
            
            static constexpr real lambda = 0.194;
            static constexpr real epsilon_0 = -2.00;
            
            __host__ __device__ inline static real SASAweight(real SASAratio){
                real w = (real(1.0)+tanhf(real(2)*tanf(SASAratio*real(M_PI_2))))/real(2.0);
                return (w<real(0.0))?real(1.0):w;
            }
        };
        
        struct F{
            
            static constexpr bool SASArequired = true;
            
            static constexpr real lambda = 0.223;
            static constexpr real epsilon_0 = -1.96;
            
            __host__ __device__ inline static real SASAweight(real SASAratio){
                real w = (real(1.0)+tanhf(tanf(SASAratio*real(M_PI_2))))/real(2.0);
                return w<real(0.0)?real(1.0):w;
            }
        };
        
        struct LQLHK{
            
            static constexpr bool SASArequired = false;
            
            static constexpr real lambda = 0.1243;
            static constexpr real epsilon_0 = -1.875;
            
            __host__ __device__ inline static real SASAweight(real SASAratio){
                return real(1.0);
            }
        };
    }}

    namespace KimHummerPotential_ns{
        namespace NonPolar{

            //Force
            inline __device__ real3 force(const real3& rij, const real& r2, const real& epsilon,const real& sigma,const real& zeroEnergy){
                
                if(epsilon == real(0.0) ){
                    return CommonPotentials::Steric::Steric::force<12>(rij,r2,zeroEnergy,sigma);
                }
            
                real3 f = CommonPotentials::LennardJones::Type1::force(rij,r2,abs(epsilon),sigma);
                //(2^(1/6)*sigma)^2
                const real r02 = (real(1.259921)*sigma*sigma);

                if(epsilon > real(0.0) and r2>=r02){
                    f=-f;
                }

                return f;
            }
            
            //Virial
            inline __device__ tensor3 virial(const real3& rij, const real& r2, const real& epsilon,const real& sigma,const real& zeroEnergy){
                return tensor3(0);
            }
            
            //Energy
            inline __device__ real energy(const real3& rij, const real& r2, const real& epsilon,const real& sigma,const real& zeroEnergy){
                
                if(epsilon == real(0.0) ){
                    return CommonPotentials::Steric::Steric::energy<12>(rij,r2,zeroEnergy,sigma);
                }
                
                real e  = CommonPotentials::LennardJones::Type1::energy(rij,r2,abs(epsilon),sigma);
            
                //(2^(1/6)*sigma)^2
                const real r02 = (real(1.259921)*sigma*sigma);

                if(epsilon > real(0.0)){
                       
                    if(r2<r02){
                        e+=epsilon; //Implicit 1/2
                    } else {
                        e=-e;
                    }
                }
            
                return e;
            }

        }
    }

    template<typename Units_>
    struct KimHummerPotential: public ParameterUpdatable{
        
        using InteractionParameters = typename CommonParameters::StatisticalPotential::InteractionParameters;
        using ParameterPairsHandler = typename structured::PairParameterHandler<InteractionParameters>;
        
        std::shared_ptr<ParticleData> pd;
        
        std::shared_ptr<ParameterPairsHandler> paramPairsHandler;
        
        Box box;

        real dielectricConstant;
        real debyeLenght;
        
        real cutOffDstNP;
        real cutOffDstDH;

        real zeroEnergy;
        
        struct Parameters{
            
            std::shared_ptr<ParameterPairsHandler> paramPairsHandler;
            
            real dielectricConstant;
            real debyeLenght;
            
            real cutOffDstNP;
            real cutOffDstDH;
        
            real zeroEnergy = real(0.01);
        
        };

        KimHummerPotential(std::shared_ptr<ParticleData> pd,
                  Parameters par):pd(pd),
                                  dielectricConstant(par.dielectricConstant),
                                  debyeLenght(par.debyeLenght),
                                  cutOffDstNP(par.cutOffDstNP),
                                  cutOffDstDH(par.cutOffDstDH),
                                  zeroEnergy(par.zeroEnergy),
                                  paramPairsHandler(par.paramPairsHandler){}

        ~KimHummerPotential(){}

        void setZeroEnergy(real newZeroEnergy){
            zeroEnergy = newZeroEnergy;
        }

        struct forceTransverser{

            real4* force;
            
            real* radius;
            real* charge;

            real*   sasa;

            typename ParameterPairsHandler::PairIterator paramPairIterator;

            Box box;
            
            real dielectricConstant;
            real debyeLenght;

            real cutOffDstNP2;
            real cutOffDstDH2;

            real zeroEnergy;
            
            forceTransverser(real4* force,
                             real* radius,
                             real* charge,
                             real* sasa,
                             typename ParameterPairsHandler::PairIterator paramPairIterator,
                             Box  box,
                             real dielectricConstant,
                             real debyeLenght,
                             real cutOffDstNP2,
                             real cutOffDstDH2,
                             real zeroEnergy):force(force),
                                              radius(radius),
                                              charge(charge),
                                              sasa(sasa),
                                              paramPairIterator(paramPairIterator),
                                              box(box),
                                              dielectricConstant(dielectricConstant),
                                              debyeLenght(debyeLenght),
                                              cutOffDstNP2(cutOffDstNP2),
                                              cutOffDstDH2(cutOffDstDH2),
                                              zeroEnergy(zeroEnergy){}

            using resultType=real4;

            inline __device__ resultType zero(){return make_real4(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                real3 f = make_real3(real(0.0));

                if(r2>0 and r2<=cutOffDstNP2){

                    const real sigma = radius[index_i]+radius[index_j];
                    const real eps   = paramPairIterator(int(posi.w),int(posj.w)).epsilon;

                    f+=KimHummerPotential_ns::NonPolar::force(rij,r2,eps,sigma,zeroEnergy);
                } 
                
                const real chgProduct = charge[index_i]*charge[index_j];
                if(r2>0 and r2<=cutOffDstDH2 and chgProduct != real(0.0)){
                    
                    f+=CommonPotentials::DebyeHuckel::DebyeHuckel::force<Units_>(rij,r2,
                                                                                 chgProduct,
                                                                                 dielectricConstant,debyeLenght);
                }

                return make_real4(sasa[index_j]*f,real(0.0));
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){force[index_i]+=sasa[index_i]*quantity;}

        };

        forceTransverser getForceTransverser(){
            
            real4* force = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();     
            
            real* radius = this->pd->getRadius(access::location::gpu, access::mode::readwrite).raw();     
            real* charge = this->pd->getCharge(access::location::gpu, access::mode::readwrite).raw();     
            real* sasa   = this->pd->getSASA(access::location::gpu, access::mode::readwrite).raw();     
                                                         
            return forceTransverser(force,
                                    radius,charge,sasa,
                                    paramPairsHandler->getPairIterator(),
                                    box,
                                    dielectricConstant,
                                    debyeLenght,
                                    cutOffDstNP*cutOffDstNP,
                                    cutOffDstDH*cutOffDstDH,
                                    zeroEnergy);
        }
        
        struct virialTransverser{

            tensor3* virial;
            
            real* radius;
            real* charge;
            real*   sasa;

            typename ParameterPairsHandler::PairIterator paramPairIterator;

            Box box;
            
            real dielectricConstant;
            real debyeLenght;

            real cutOffDstNP2;
            real cutOffDstDH2;

            real zeroEnergy;
            
            virialTransverser(tensor3* virial,
                              real* radius,
                              real* charge,
                              real* sasa,
                              typename ParameterPairsHandler::PairIterator paramPairIterator,
                              Box  box,
                              real dielectricConstant,
                              real debyeLenght,
                              real cutOffDstNP2,
                              real cutOffDstDH2,
                              real zeroEnergy):virial(virial),
                                               radius(radius),
                                               charge(charge),
                                               sasa(sasa),
                                               paramPairIterator(paramPairIterator),
                                               box(box),
                                               dielectricConstant(dielectricConstant),
                                               debyeLenght(debyeLenght),
                                               cutOffDstNP2(cutOffDstNP2),
                                               cutOffDstDH2(cutOffDstDH2),
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
            real* sasa   = this->pd->getSASA(access::location::gpu, access::mode::readwrite).raw();     
                                                         
            return virialTransverser(virial,
                                     radius,charge,sasa,
                                     paramPairsHandler->getPairIterator(),
                                     box,
                                     dielectricConstant,
                                     debyeLenght,
                                     cutOffDstNP*cutOffDstNP,
                                     cutOffDstDH*cutOffDstDH,
                                     zeroEnergy);
        }
        
        struct energyTransverser{

            real* energy;
            
            real* radius;
            real* charge;
            real*   sasa;

            typename ParameterPairsHandler::PairIterator paramPairIterator;

            Box box;
            
            real dielectricConstant;
            real debyeLenght;

            real cutOffDstNP2;
            real cutOffDstDH2;

            real zeroEnergy;
            
            energyTransverser(real* energy,
                              real* radius,
                              real* charge,
                              real* sasa,
                              typename ParameterPairsHandler::PairIterator paramPairIterator,
                              Box  box,
                              real dielectricConstant,
                              real debyeLenght,
                              real cutOffDstNP2,
                              real cutOffDstDH2,
                              real zeroEnergy):energy(energy),
                                               radius(radius),
                                               charge(charge),
                                               sasa(sasa),
                                               paramPairIterator(paramPairIterator),
                                               box(box),
                                               dielectricConstant(dielectricConstant),
                                               debyeLenght(debyeLenght),
                                               cutOffDstNP2(cutOffDstNP2),
                                               cutOffDstDH2(cutOffDstDH2),
                                               zeroEnergy(zeroEnergy){}

            using resultType=real;

            inline __device__ resultType zero(){return real(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                real e = real(0.0);

                if(r2>0 and r2<=cutOffDstNP2){

                    const real sigma = radius[index_i]+radius[index_j];
                    const real eps   = paramPairIterator(int(posi.w),int(posj.w)).epsilon;

                    e+=KimHummerPotential_ns::NonPolar::energy(rij,r2,eps,sigma,zeroEnergy);
                } 
                
                const real chgProduct = charge[index_i]*charge[index_j];
                if(r2>0 and r2<=cutOffDstDH2 and chgProduct != real(0.0)){
                    
                    e+=CommonPotentials::DebyeHuckel::DebyeHuckel::energy<Units_>(rij,r2,
                                                                                  chgProduct,
                                                                                  dielectricConstant,debyeLenght);
                }

                return sasa[index_j]*e;
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){energy[index_i]+=sasa[index_i]*quantity;}

        };

        energyTransverser getEnergyTransverser(){
            
            real* energy = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();     
            
            real* radius = this->pd->getRadius(access::location::gpu, access::mode::readwrite).raw();     
            real* charge = this->pd->getCharge(access::location::gpu, access::mode::readwrite).raw();     
            real* sasa   = this->pd->getSASA(access::location::gpu, access::mode::readwrite).raw();     
                                                         
            return energyTransverser(energy,
                                     radius,charge,sasa,
                                     paramPairsHandler->getPairIterator(),
                                     box,
                                     dielectricConstant,
                                     debyeLenght,
                                     cutOffDstNP*cutOffDstNP,
                                     cutOffDstDH*cutOffDstDH,
                                     zeroEnergy);
        }
        
        void updateBox(Box newBox) override {
            box=newBox;
        }
    };

}}}}

#endif
