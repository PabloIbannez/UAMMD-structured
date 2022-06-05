#ifndef __LENNARD_JONES_POTENTIAL__
#define __LENNARD_JONES_POTENTIAL__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace UnBound{
    
    template<class Topology,class LennardJonesType>
    struct LennardJones_: public ParameterUpdatable{
        
        using ParametersType        = typename CommonParameters::LennardJones::LennardJones<Topology>;
        using ParameterPairsHandler = typename structured::PairParameterHandler<typename ParametersType::InteractionParameters>;
        
        std::shared_ptr<System>        sys;
        std::shared_ptr<ParticleData>  pd;
        std::shared_ptr<ParticleGroup> pg;
        std::shared_ptr<Topology>      top;

        std::shared_ptr<ParametersType> ljParam;
        
        Box box;

        std::string label;
        real cutOffDst;
        
        struct Parameters{
            
            std::string label;
            real cutOffDst;
        };

        LennardJones_(std::shared_ptr<System>       sys,
                      std::shared_ptr<ParticleData>  pd,
                      std::shared_ptr<ParticleGroup> pg,
                      std::shared_ptr<Topology>     top,
                      Parameters par):sys(sys),
                                      pd(pd),
                                      pg(pg),
                                      top(top),
                                      label(par.label),
                                      cutOffDst(par.cutOffDst){
            
            typename ParametersType::Parameters param;

            param.label = label;

            ljParam = std::make_shared<ParametersType>(sys,pd,pg,top,param);
        }

        real getCutOffDst(){
            return cutOffDst;
        }

        ~LennardJones_(){}

        struct forceTransverser{
            
            real4* force;
            
            typename ParameterPairsHandler::PairIterator paramPairIterator;

            Box box;
            
            real cutOffDst2;
            
            forceTransverser(real4* force,
                             typename ParameterPairsHandler::PairIterator paramPairIterator,
                             Box  box,
                             real cutOffDst2):force(force),
                                             paramPairIterator(paramPairIterator),
                                             box(box),
                                             cutOffDst2(cutOffDst2){}

            using resultType=real4;

            inline __device__ resultType zero(){return make_real4(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int& index_i,const int& index_j,const real4& posi,const real4& posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));

                const real epsilon = paramPairIterator(int(posi.w),int(posj.w)).epsilon;
                const real sigma   = paramPairIterator(int(posi.w),int(posj.w)).sigma;
            
                const real r2 = dot(rij, rij);
                
                real3 f = make_real3(0.0);

                if(r2<=cutOffDst2){
                    f = LennardJonesType::force(rij,r2,epsilon,sigma);
                }

                return make_real4(f,0.0);

            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){force[index_i]+=quantity;}

        };

        forceTransverser getForceTransverser(){

            real4* force = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();     
            
            return forceTransverser(force,
                                    ljParam->getParameters()->getPairIterator(),
                                    box,
                                    cutOffDst*cutOffDst);
        }
        
        struct virialTransverser{

            tensor3* virial;
            
            typename ParameterPairsHandler::PairIterator paramPairIterator;

            Box box;
            
            real cutOffDst2;
            
            virialTransverser(tensor3* virial,
                              typename ParameterPairsHandler::PairIterator paramPairIterator,
                              Box  box,
                              real cutOffDst2):virial(virial),
                                                       paramPairIterator(paramPairIterator),
                                                       box(box),
                                                       cutOffDst2(cutOffDst2){}
            using resultType=tensor3;

            inline __device__ resultType zero(){return tensor3(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                tensor3 v(0.0);

                return v;
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){virial[index_i]+=quantity;}

        };

        virialTransverser getVirialTransverser(){
            
            tensor3* virial = this->pd->getVirial(access::location::gpu, access::mode::readwrite).raw();     
            
            return virialTransverser(virial,
                                     ljParam->getParameters()->getPairIterator(),
                                     box,
                                     cutOffDst*cutOffDst);
        }
        
        struct energyTransverser{

            real* energy;
            
            typename ParameterPairsHandler::PairIterator paramPairIterator;

            Box box;
            
            real cutOffDst2;
            
            energyTransverser(real* energy,
                              typename ParameterPairsHandler::PairIterator paramPairIterator,
                              Box  box,
                              real cutOffDst2):energy(energy),
                                                       paramPairIterator(paramPairIterator),
                                                       box(box),
                                                       cutOffDst2(cutOffDst2){}

            using resultType=real;

            inline __device__ resultType zero(){return real(0.0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                    
                const real epsilon = paramPairIterator(int(posi.w),int(posj.w)).epsilon;
                const real sigma   = paramPairIterator(int(posi.w),int(posj.w)).sigma;

                const real r2 = dot(rij, rij);
                
                real e = real(0.0);

                if(r2<=cutOffDst2){
                    e = LennardJonesType::energy(rij,r2,epsilon,sigma);
                }

                return e;
                
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){energy[index_i]+=quantity;}

        };

        energyTransverser getEnergyTransverser(){
            
            real* energy = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();     
                                                         
            return energyTransverser(energy,
                                     ljParam->getParameters()->getPairIterator(),
                                     box,
                                     cutOffDst*cutOffDst);
        }
        
        void updateBox(Box newBox) override {
            box=newBox;
        }

    };
    
    template<class Topology>
    using LennardJonesType1 = LennardJones_<Topology,CommonPotentials::LennardJones::Type1>;
    template<class Topology>
    using LennardJonesType2 = LennardJones_<Topology,CommonPotentials::LennardJones::Type2>;
    template<class Topology>
    using LennardJonesType3 = LennardJones_<Topology,CommonPotentials::LennardJones::Type3>;
    
    template<class Topology>
    using WCAType1 = LennardJones_<Topology,CommonPotentials::WCA::Type1>;
    template<class Topology>
    using WCAType2 = LennardJones_<Topology,CommonPotentials::WCA::Type2>;
    template<class Topology>
    using WCAType3 = LennardJones_<Topology,CommonPotentials::WCA::Type3>;
    
    template<class Topology>
    using GeneralLennardJonesType1 = LennardJones_<Topology,CommonPotentials::GeneralLennardJones::Type1>;
    template<class Topology>
    using GeneralLennardJonesType2 = LennardJones_<Topology,CommonPotentials::GeneralLennardJones::Type2>;
    template<class Topology>
    using GeneralLennardJonesType3 = LennardJones_<Topology,CommonPotentials::GeneralLennardJones::Type3>;
    
    template<class Topology>
    using Steric6  = LennardJones_<Topology,CommonPotentials::Steric::Steric6>;
    template<class Topology>
    using Steric12 = LennardJones_<Topology,CommonPotentials::Steric::Steric12>;

}}}}

#endif
