#ifndef __LENNARD_JONES_POTENTIAL__
#define __LENNARD_JONES_POTENTIAL__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace UnBound{
    
    namespace LennardJones_ns{
    
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
    
    struct LennardJones: public ParameterUpdatable{
        
        using InteractionParameters = typename LennardJones_ns::InteractionParameters;
        using ParameterPairsHandler = typename structured::PairParameterHandler<InteractionParameters>;
        
        std::shared_ptr<ParticleData> pd;
        
        std::shared_ptr<ParameterPairsHandler> paramPairsHandler;
        
        Box box;

        real cutOffDstLJ;
            
        struct Parameters{
            
            std::shared_ptr<ParameterPairsHandler> paramPairsHandler;
            
            real cutOffDstLJ;
        
        };

        LennardJones(std::shared_ptr<ParticleData> pd,
                     Parameters par):pd(pd),
                                     cutOffDstLJ(par.cutOffDstLJ),
                                     paramPairsHandler(par.paramPairsHandler){}

        ~LennardJones(){}

        struct forceTransverser{
            
            real4* force;
            
            typename ParameterPairsHandler::PairIterator paramPairIterator;

            Box box;
            
            real cutOffDstLJ2;
            
            forceTransverser(real4* force,
                             typename ParameterPairsHandler::PairIterator paramPairIterator,
                             Box  box,
                             real cutOffDstLJ2):force(force),
                                             paramPairIterator(paramPairIterator),
                                             box(box),
                                             cutOffDstLJ2(cutOffDstLJ2){}

            using resultType=real4;

            inline __device__ resultType zero(){return make_real4(0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int& index_i,const int& index_j,const real4& posi,const real4& posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));

                const real epsilon = paramPairIterator(int(posi.w),int(posj.w)).epsilon;
                const real sigma   = paramPairIterator(int(posi.w),int(posj.w)).sigma;
            
                const real r2 = dot(rij, rij);
                
                real3 f = make_real3(0.0);

                if(r2<=cutOffDstLJ2){
                    f = CommonPotentials::LennardJones::Type1::force(rij,r2,epsilon,sigma);
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
                                    cutOffDstLJ*cutOffDstLJ);
        }
        
        struct virialTransverser{

            tensor3* virial;
            
            typename ParameterPairsHandler::PairIterator paramPairIterator;

            Box box;
            
            real cutOffDstLJ2;
            
            virialTransverser(tensor3* virial,
                              typename ParameterPairsHandler::PairIterator paramPairIterator,
                              Box  box,
                              real cutOffDstLJ2):virial(virial),
                                                       paramPairIterator(paramPairIterator),
                                                       box(box),
                                                       cutOffDstLJ2(cutOffDstLJ2){}
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
                                     cutOffDstLJ*cutOffDstLJ);
        }
        
        struct energyTransverser{

            real* energy;
            
            typename ParameterPairsHandler::PairIterator paramPairIterator;

            Box box;
            
            real cutOffDstLJ2;
            
            energyTransverser(real* energy,
                              typename ParameterPairsHandler::PairIterator paramPairIterator,
                              Box  box,
                              real cutOffDstLJ2):energy(energy),
                                                       paramPairIterator(paramPairIterator),
                                                       box(box),
                                                       cutOffDstLJ2(cutOffDstLJ2){}

            using resultType=real;

            inline __device__ resultType zero(){return real(0.0);}
            
            inline __device__ void accumulate(resultType& total,const resultType current){total+=current;}
            
            inline __device__ resultType compute(const int index_i,const int index_j,const real4 posi,const real4 posj){
                
                const real3 rij = box.apply_pbc(make_real3(posj)-make_real3(posi));
                    
                const real epsilon = paramPairIterator(int(posi.w),int(posj.w)).epsilon;
                const real sigma   = paramPairIterator(int(posi.w),int(posj.w)).sigma;

                const real r2 = dot(rij, rij);
                
                real e = real(0.0);

                if(r2<=cutOffDstLJ2){
                    e = CommonPotentials::LennardJones::Type1::energy(rij,r2,epsilon,sigma);
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
                                     cutOffDstLJ*cutOffDstLJ);
        }
        
        void updateBox(Box newBox) override {
            box=newBox;
        }

    };

}}}}

#endif
