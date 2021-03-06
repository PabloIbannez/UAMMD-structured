#ifndef __GENERIC_SURFACE_POT__
#define __GENERIC_SURFACE_POT__
    
namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Surface{

    namespace GenericSurfacePotential_ns{
    namespace GenericSurfacePotentialModel{

        struct Steric{
            static constexpr real A = 1.0;
            static constexpr real B = 0.0;
            static constexpr real energyOffSet = 0.0;
            static constexpr real cutOff = INFINITY;
        };
        
        struct LennardJones{
            static constexpr real A =  1.0;
            static constexpr real B = -2.0;
            static constexpr real energyOffSet = 0.0;
            static constexpr real cutOff = INFINITY;
        };

    }}
    
    template<class SurfaceType> 
    class GenericSurfacePotential: public ParameterUpdatable{
    	
        private:
            
            real epsilon;
            real sigma;
            
            real surfacePosition;
    	        
        public:

            struct Parameters{
            
                real epsilon;
                real sigma;

                real surfacePosition;
            };
    	    
    	    GenericSurfacePotential(Parameters par):epsilon(par.epsilon),
                                                    sigma(par.sigma),
                                                    surfacePosition(par.surfacePosition){}     
    	    
    	    __device__ __forceinline__ real3 force(const real4 &pos){
                
                const real dz = abs(surfacePosition-pos.z);

                if(dz<=SurfaceType::cutOff){

                    const real invdz2  = real(1.0)/(dz*dz);
                    
                    const real sinvdz2  = sigma*sigma*invdz2;
                    const real sinvdz6  = sinvdz2*sinvdz2*sinvdz2;
                    const real sinvdz12 = sinvdz6*sinvdz6;

                    real fmod = epsilon*real(6.0)*(real(2.0)*SurfaceType::A*sinvdz12+SurfaceType::B*sinvdz6)*invdz2; 

                    return {real(0),
                            real(0),
                            -fmod*(surfacePosition-pos.z)};
                }

                return make_real3(0);
            
            }
    	    
            __device__ __forceinline__ real  energy(const real4 &pos){

                const real dz = abs(surfacePosition-pos.z);

                if(dz<=SurfaceType::cutOff){

                    const real sinvdz2  = sigma*sigma/(dz*dz);
                    const real sinvdz6  = sinvdz2*sinvdz2*sinvdz2;
                    const real sinvdz12 = sinvdz6*sinvdz6;

                    real e = epsilon*(SurfaceType::A*sinvdz12+SurfaceType::B*sinvdz6)+SurfaceType::energyOffSet; 

                    return e;
                }

                return 0;
            }
            __device__ __forceinline__ real virial(const real4 &pos){return real(0);}
            __device__ __forceinline__ tensor3 stress(const real4 &pos){return tensor3(0);}

            __device__ __forceinline__ ForceEnergyVirial sum(Interactor::Computables comp,const real4& pos){

                real3 f = comp.force?force(pos):real3();
                real  e = comp.energy?energy(pos):real(0.0);
                real  v = comp.virial?virial(pos):real(0.0);

                return {f,e,v};
            }
    	    
    	    std::tuple<const real4 *> getArrays(ParticleData *pd){
    	    	
                auto pos = pd->getPos(access::location::gpu, access::mode::read);
    	    	
                return std::make_tuple(pos.raw());
    	    }

            real getSurfacePosition(){return surfacePosition;}
            void setSurfacePosition(real newSurfacePosition){surfacePosition=newSurfacePosition;}

    };
    
}}}}

#endif
