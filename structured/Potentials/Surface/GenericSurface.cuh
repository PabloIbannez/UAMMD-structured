#ifndef __GENERIC_SURFACE_POT__
#define __GENERIC_SURFACE_POT__
    
namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Surface{
    
    class GenericSurfacePotential: public ParameterUpdatable{
    	
        private:
            
            real epsilonSurf;
            real sigmaSurf;
            
            real Asurf;
            real Bsurf;

            real offSetSurf;
            real cutOffSurf;
    	    
            real surfacePosition;
    	        
        public:

            struct Parameters{
            
                real epsilonSurf;
                real sigmaSurf;

                real Asurf;
                real Bsurf;
            
                real offSetSurf;
                real cutOffSurf;
    	        
                real surfacePosition;
            };
    	    
    	    GenericSurfacePotential(Parameters par):epsilonSurf(par.epsilonSurf),
                                                    sigmaSurf(par.sigmaSurf),
                                                    Asurf(par.Asurf),
                                                    Bsurf(par.Bsurf),      
                                                    offSetSurf(par.offSetSurf),
                                                    cutOffSurf(par.cutOffSurf),
                                                    surfacePosition(par.surfacePosition){}     
    	    
    	    __device__ __forceinline__ real3 force(const real4 &pos){
                
                const int type = int(pos.w);
    	    	
                const real dz = abs(surfacePosition-pos.z);

                if(dz<=cutOffSurf){

                    const real invdz2  = real(1.0)/(dz*dz);
                    
                    const real sinvdz2  = sigmaSurf*sigmaSurf*invdz2;
                    const real sinvdz6  = sinvdz2*sinvdz2*sinvdz2;
                    const real sinvdz12 = sinvdz6*sinvdz6;

                    real fmod = epsilonSurf*real(6.0)*(real(2.0)*Asurf*sinvdz12+Bsurf*sinvdz6)*invdz2; 

                    return {real(0),
                            real(0),
                            -fmod*(surfacePosition-pos.z)};
                }

                return make_real3(0);
            
            }
    	    
            __device__ __forceinline__ real  energy(const real4 &pos){

                const int type = int(pos.w);
    	    	
                const real dz = abs(surfacePosition-pos.z);

                if(dz<=cutOffSurf){

                    const real sinvdz2  = sigmaSurf*sigmaSurf/(dz*dz);
                    const real sinvdz6  = sinvdz2*sinvdz2*sinvdz2;
                    const real sinvdz12 = sinvdz6*sinvdz6;

                    real e = epsilonSurf*(Asurf*sinvdz12+Bsurf*sinvdz6)+offSetSurf; 

                    return e;
                }

                return 0;
            }
            __device__ __forceinline__ real3 virial(const real4 &pos){return make_real3(0);}
    	    
    	    std::tuple<const real4 *> getArrays(ParticleData *pd){
    	    	
                auto pos = pd->getPos(access::location::gpu, access::mode::read);
    	    	
                return std::make_tuple(pos.raw());
    	    }
    };
    
}}}}

#endif
