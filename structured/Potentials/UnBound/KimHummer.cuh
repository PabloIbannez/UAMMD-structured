#ifndef __KIM_HUMMER_UNBOUND__
#define __KIM_HUMMER_UNBOUND__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace UnBound{
    
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

    template<class Topology>
    struct KimHummer_: public ParameterUpdatable{
        
        using ParametersType        = typename CommonParameters::StatisticalPotential::StatisticalPotential<Topology>;
        using ParameterPairsHandler = typename structured::PairParameterHandler<typename ParametersType::InteractionParameters>;
        
        std::shared_ptr<System>        sys;
        std::shared_ptr<ParticleData>  pd;
        std::shared_ptr<ParticleGroup> pg;
        std::shared_ptr<Topology>      top;

        std::shared_ptr<ParametersType> statParam;
        
        Box box;

        real dielectricConstant;
        real debyeLenght;
        
        real refTemperature;
        real epsilon_0;
        real lambda;
        std::string sasaModel;
        
        std::string label;
        std::string sasaLabel;
        
        real cutOffDstNP;
        real cutOffDstDH;
        
        real zeroEnergy;
        
        struct Parameters{
            
            real dielectricConstant;
            real debyeLenght;
            
            real refTemperature;
            real epsilon_0;
            real lambda;
            std::string sasaModel;
            
            std::string label;
            std::string sasaLabel;
            
            real cutOffDstNP;
            real cutOffDstDH;
        
            real zeroEnergy;
        
        };

        KimHummer_(std::shared_ptr<System>       sys,
                   std::shared_ptr<ParticleData>  pd,
                   std::shared_ptr<ParticleGroup> pg,
                   std::shared_ptr<Topology>     top,
                   Parameters par):sys(sys),
                                   pd(pd),
                                   pg(pg),
                                   top(top),
                                   dielectricConstant(par.dielectricConstant),
                                   debyeLenght(par.debyeLenght),
                                   refTemperature(par.refTemperature),
                                   epsilon_0(par.epsilon_0),
                                   lambda(par.lambda),
                                   sasaModel(par.sasaModel),
                                   label(par.label),
                                   sasaLabel(par.sasaLabel),
                                   cutOffDstNP(par.cutOffDstNP),
                                   cutOffDstDH(par.cutOffDstDH),
                                   zeroEnergy(par.zeroEnergy){
            
            typename ParametersType::Parameters param;

            param.label          = label;
            param.refTemperature = refTemperature;
            param.epsilon_0      = epsilon_0;
            param.lambda         = lambda;

            statParam = std::make_shared<ParametersType>(sys,pd,pg,top,param);

            //Check sasa model
            if(sasaModel != "A"){
                this->top->template loadProperty(pd,sasaLabel,pd->getSASA(access::location::cpu, access::mode::write));
                {
                    auto groupIndex = this->pg->getIndexIterator(access::location::cpu);
                    auto SASA = this->pd->getSASA(access::location::cpu, access::mode::readwrite);
                    
                    fori(0,this->pg->getNumberParticles()){
                        int  index = groupIndex[i];
                        
                        if(SASA[index] > real(1.0)){
                            SASA[index] = real(1.0);
                        }

                        if        (sasaModel == "B"){
                            SASA[index]=KimHummer_ns::SasaModel::B::SASAweight(SASA[index]);
                        } else if (sasaModel == "C"){
                            SASA[index]=KimHummer_ns::SasaModel::C::SASAweight(SASA[index]);
                        } else if (sasaModel == "D"){
                            SASA[index]=KimHummer_ns::SasaModel::D::SASAweight(SASA[index]);
                        } else if (sasaModel == "E"){
                            SASA[index]=KimHummer_ns::SasaModel::E::SASAweight(SASA[index]);
                        } else if (sasaModel == "F"){
                            SASA[index]=KimHummer_ns::SasaModel::F::SASAweight(SASA[index]);
                        } else {
                            this->sys->template log<System::CRITICAL>("[KimHummerPotential] Sasa model: %s is no available ",sasaModel.c_str());
                        }
                    }
                }
            } else {
                    
                auto groupIndex = this->pg->getIndexIterator(access::location::cpu);
                auto SASA = this->pd->getSASA(access::location::cpu, access::mode::readwrite);
                    
                fori(0,this->pg->getNumberParticles()){
                    int  index = groupIndex[i];
                    SASA[index]=real(1.0);
                }
            }
        }
        
        real getCutOffDst(){
            return std::max(cutOffDstNP,cutOffDstDH);
        }

        ~KimHummer_(){}

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

                    f+=CommonPotentials::NonPolar::force(rij,r2,eps,sigma,zeroEnergy);
                } 
                
                const real chgProduct = charge[index_i]*charge[index_j];
                if(r2>0 and r2<=cutOffDstDH2 and chgProduct != real(0.0)){
                    
                    f+=CommonPotentials::DebyeHuckel::DebyeHuckel::force<Topology::Units>(rij,r2,
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
                                    statParam->getParameters()->getPairIterator(),
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
                                     statParam->getParameters()->getPairIterator(),
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

                    e+=CommonPotentials::NonPolar::energy(rij,r2,eps,sigma,zeroEnergy);
                } 
                
                const real chgProduct = charge[index_i]*charge[index_j];
                if(r2>0 and r2<=cutOffDstDH2 and chgProduct != real(0.0)){
                    
                    e+=CommonPotentials::DebyeHuckel::DebyeHuckel::energy<Topology::Units>(rij,r2,
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
                                     statParam->getParameters()->getPairIterator(),
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
   
    template<class Topology>
    using KimHummer = KimHummer_<Topology>;

}}}}

#endif
