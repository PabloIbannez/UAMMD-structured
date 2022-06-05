#ifndef __LJ_WCA_POTENTIAL__
#define __LJ_WCA_POTENTIAL__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace UnBound{
    
    namespace LJ_WCA_ns{
    
        struct InteractionParameters{
        
            struct InputPairParameters{
                real epsilon;
                real sigma;
            };
        
            struct PairParameters{
                real epsilon;
                real sigma;
            };
        
            static inline __host__ InputPairParameters readPairParameters(std::string& line){
            
                std::stringstream ss;
                
                InputPairParameters param;

                ss.str(line);

                ss  >> param.epsilon >> param.sigma;
                
                return param;
            
            }
        
            static inline __host__ PairParameters processPairParameters(InputPairParameters in_par){
        
                PairParameters params;
        
                params.epsilon = in_par.epsilon;
                params.sigma   = in_par.sigma;
        
                return params;
            }
        };
    }
    
    struct LJ_WCA: public ParameterUpdatable{
        
        using InteractionParameters = typename LJ_WCA_ns::InteractionParameters;
        using ParameterPairsHandler = typename structured::PairParameterHandler<InteractionParameters>;
        
        std::shared_ptr<ParticleData> pd;
        
        std::shared_ptr<ParameterPairsHandler> paramPairsHandler;
        
        Box box;

        real cutOffDst;
            
        struct Parameters{
            
            std::shared_ptr<ParameterPairsHandler> paramPairsHandler;
            
            real cutOffDst;
        
        };

        LJ_WCA(std::shared_ptr<ParticleData> pd,
               Parameters par):pd(pd),
                               cutOffDst(par.cutOffDst),
                               paramPairsHandler(par.paramPairsHandler){}

        ~LJ_WCA(){}

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
                    if(epsilon < real(0.0)){
                        f = CommonPotentials::LennardJones::Type3::force(rij,r2,abs(epsilon),sigma);
                    } else {
                        f = CommonPotentials::WCA::Type2::force(rij,r2,epsilon,sigma);
                    }
                }

                return make_real4(f,0.0);

            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){force[index_i]+=quantity;}

        };

        forceTransverser getForceTransverser(){

            real4* force = this->pd->getForce(access::location::gpu, access::mode::readwrite).raw();     
            
            return forceTransverser(force,
                                    paramPairsHandler->getPairIterator(),
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
                                     paramPairsHandler->getPairIterator(),
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
                    if(epsilon < real(0.0)){
                        e = CommonPotentials::LennardJones::Type3::energy(rij,r2,abs(epsilon),sigma);
                    } else {
                        e = CommonPotentials::WCA::Type2::energy(rij,r2,epsilon,sigma);
                    }
                }

                return e;
                
            }
            
            inline __device__ void set(const int& index_i,const resultType& quantity){energy[index_i]+=quantity;}

        };

        energyTransverser getEnergyTransverser(){
            
            real* energy = this->pd->getEnergy(access::location::gpu, access::mode::readwrite).raw();     
                                                         
            return energyTransverser(energy,
                                     paramPairsHandler->getPairIterator(),
                                     box,
                                     cutOffDst*cutOffDst);
        }
        
        void updateBox(Box newBox) override {
            box=newBox;
        }

    };

}}}}

#endif
