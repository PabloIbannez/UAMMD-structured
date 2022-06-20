#ifndef __DLVO_UNBOUND__
#define __DLVO_UNBOUND__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace UnBound{

    template<class Topology,class LennardJonesType>
    struct DLVO_: public ParameterUpdatable{
        
        using ParametersType        = typename CommonParameters::LennardJones::LennardJones<Topology>;
        using ParameterPairsHandler = typename structured::PairParameterHandler<typename ParametersType::InteractionParameters>;
        
        std::shared_ptr<System>       sys;
        std::shared_ptr<ParticleData>  pd;
        std::shared_ptr<ParticleGroup> pg;
        std::shared_ptr<Topology>     top;
        
        std::shared_ptr<ParametersType> ljParam;
        
        Box box;
            
        real dielectricConstant; 
        real debyeLength;

        real cutOffDstPolar;
        
        std::string label;
        real cutOffDstNonPolar;
            
        struct Parameters{

            real dielectricConstant; 
            real debyeLength;

            real cutOffDstPolar;
            
            std::string label;
            real cutOffDstNonPolar;
        };

        DLVO_(std::shared_ptr<ParticleGroup> pg,
              std::shared_ptr<Topology>     top,
              Parameters par):pg(pg),
                              pd(pg->getParticleData()),
                              sys(pg->getParticleData()->getSystem()),
                              dielectricConstant(par.dielectricConstant),
                              debyeLength(par.debyeLength),
                              cutOffDstPolar(par.cutOffDstPolar),
                              label(par.label),
                              cutOffDstNonPolar(par.cutOffDstNonPolar){
            
            typename ParametersType::Parameters param;

            param.label = label;

            ljParam = std::make_shared<ParametersType>(sys,pd,pg,top,param);
        }
        
        real getCutOffDst(){
            return std::max(cutOffDstPolar,cutOffDstNonPolar);
        }

        ~DLVO_(){}

        struct forceTransverser{

            real4* force;
            
            real* radius;
            real* charge;
            
            typename ParameterPairsHandler::PairIterator paramPairIterator;

            real dielectricConstant; 
            real debyeLength;

            Box box;

            real cutOffDstPolar2;
            real cutOffDstNonPolar2;

            forceTransverser(real4* force,
                             real* radius,
                             real* charge,
                             typename ParameterPairsHandler::PairIterator paramPairIterator,
                             real dielectricConstant, 
                             real debyeLength,
                             Box box,
                             real cutOffDstPolar2,
                             real cutOffDstNonPolar2):force(force),
                                                      radius(radius),
                                                      charge(charge),
                                                      paramPairIterator(paramPairIterator),
                                                      dielectricConstant(dielectricConstant),
                                                      debyeLength(debyeLength),
                                                      box(box),
                                                      cutOffDstPolar2(cutOffDstPolar2),
                                                      cutOffDstNonPolar2(cutOffDstNonPolar2){}

            using resultType=real4;

            inline __device__ resultType zero(){return make_real4(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType& current){total+=current;}
            
            inline __device__ resultType compute(const int& index_i,const int& index_j,const real4& posi,const real4& posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                real3 f = make_real3(0.0);

                if(r2<=cutOffDstPolar2){

                    const real chgProduct = charge[index_i]*charge[index_j];

                    f += CommonPotentials::DebyeHuckel::DebyeHuckelSpheres::force<Topology::Units>(rij,r2,chgProduct,
                                                                                                   radius[index_i],radius[index_j],
                                                                                                   dielectricConstant,debyeLength);
                } 
                
                if(r2<=cutOffDstNonPolar2){
                    
                    const real epsilon = paramPairIterator(int(posi.w),int(posj.w)).epsilon;
                    const real sigma   = paramPairIterator(int(posi.w),int(posj.w)).sigma;

                    f += LennardJonesType::force(rij,r2,epsilon,sigma);
                }

                return make_real4(f,0.0);
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){force[index_i]+=quantity;}

        };

        forceTransverser getForceTransverser(){
            
            real4* force  = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();     
            real*  radius = this->pd->getRadius(access::location::gpu, access::mode::readwrite).raw();     
            real*  charge = this->pd->getCharge(access::location::gpu, access::mode::readwrite).raw();     

            return forceTransverser(force,
                                    radius,
                                    charge,
                                    ljParam->getParameters()->getPairIterator(),
                                    dielectricConstant, 
                                    debyeLength,
                                    box,
                                    cutOffDstPolar*cutOffDstPolar,
                                    cutOffDstNonPolar*cutOffDstNonPolar);
        }
        
        struct energyTransverser{

            real* energy;
            
            real* radius;
            real* charge;
            
            typename ParameterPairsHandler::PairIterator paramPairIterator;

            real dielectricConstant; 
            real debyeLength;

            Box box;

            real cutOffDstPolar2;
            real cutOffDstNonPolar2;

            energyTransverser(real* energy,
                              real* radius,
                              real* charge,
                              typename ParameterPairsHandler::PairIterator paramPairIterator,
                              real dielectricConstant, 
                              real debyeLength,
                              Box box,
                              real cutOffDstPolar2,
                              real cutOffDstNonPolar2):energy(energy),
                                                       radius(radius),
                                                       charge(charge),
                                                       paramPairIterator(paramPairIterator),
                                                       dielectricConstant(dielectricConstant),
                                                       debyeLength(debyeLength),
                                                       box(box),
                                                       cutOffDstPolar2(cutOffDstPolar2),
                                                       cutOffDstNonPolar2(cutOffDstNonPolar2){}

            using resultType=real;

            inline __device__ resultType zero(){return real(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType& current){total+=current;}
            
            inline __device__ resultType compute(const int& index_i,const int& index_j,const real4& posi,const real4& posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                const real r2   = dot(rij, rij);
                
                real e = real(0.0);

                if(r2<=cutOffDstPolar2){

                    const real chgProduct = charge[index_i]*charge[index_j];

                    e += CommonPotentials::DebyeHuckel::DebyeHuckelSpheres::energy<Topology::Units>(rij,r2,chgProduct,
                                                                                                    radius[index_i],radius[index_j],
                                                                                                    dielectricConstant,debyeLength);
                } 
                
                if(r2<=cutOffDstNonPolar2){
                    
                    const real epsilon = paramPairIterator(int(posi.w),int(posj.w)).epsilon;
                    const real sigma   = paramPairIterator(int(posi.w),int(posj.w)).sigma;

                    e += LennardJonesType::energy(rij,r2,epsilon,sigma);
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
                                     radius,
                                     charge,
                                     ljParam->getParameters()->getPairIterator(),
                                     dielectricConstant, 
                                     debyeLength,
                                     box,
                                     cutOffDstPolar*cutOffDstPolar,
                                     cutOffDstNonPolar*cutOffDstNonPolar);
        }

        void updateBox(Box newBox) override {
            box=newBox;
        }

    };
    
    template<class Topology>
    using DLVOType1 = DLVO_<Topology,CommonPotentials::LennardJones::Type1>;
    template<class Topology>
    using DLVOType2 = DLVO_<Topology,CommonPotentials::LennardJones::Type2>;
    template<class Topology>
    using DLVOType3 = DLVO_<Topology,CommonPotentials::LennardJones::Type3>;

}}}}

#endif
