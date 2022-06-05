#ifndef __KARANICOLAS_BROOKS_SURFACE_POTENTIAL__
#define __KARANICOLAS_BROOKS_SURFACE_POTENTIAL__
    
namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Surface{
    
    class KaranicolasBrooksSurfacePotential: public ParameterUpdatable{
    	
        private:

            static constexpr real theta1 = 0.2340;
            static constexpr real theta2 = 0.4936;
            static constexpr real theta3 = 0.1333;
            static constexpr real thetaS = 0.0067;
            static constexpr real thetaP = 0.0333;
            
            real epsilonSurface;
            real sigmaSurface;
                
            real rho;
                
            real chiSurface;
    	    
            real surfacePosition;

            bool useInnerRadius;
    	        
        public:

            struct Parameters{
                
                real epsilonSurface;
                real sigmaSurface;
    	        
                real rho;
                
                real chiSurface;
                
                real surfacePosition;

                bool useInnerRadius=false;
            
            };
    	    
    	    KaranicolasBrooksSurfacePotential(Parameters par):epsilonSurface(par.epsilonSurface),     
                                                              sigmaSurface(par.sigmaSurface),
                                                              rho(par.rho),
                                                              chiSurface(par.chiSurface),
                                                              surfacePosition(par.surfacePosition),
        useInnerRadius(par.useInnerRadius){}     
    	    
    	    __device__ __forceinline__ real3 force(const real4& pos,const real& radius,const real& epsilon, real& surface){
                
                const real dz = abs(surfacePosition-pos.z);
                
                const real sigma     = (sigmaSurface+radius)/real(2.0);
                const real epsilonLB = sqrt(epsilonSurface*epsilon);

                const real sigma3  = sigma*sigma*sigma;
                const real sinvdz3 = sigma3/(dz*dz*dz);
                const real sinvdz7 = sinvdz3*sinvdz3*(sigma/dz);
                const real sinvdz9 = sinvdz3*sinvdz3*sinvdz3;

                const real ccommon = real(M_PI)*rho*sigma3*epsilonLB;
                const real c9      = ccommon*theta1;
                const real c7      = ccommon*theta2;
                const real c3      = ccommon*theta3;
                const real c3hydro = ccommon*(thetaS*chiSurface+thetaP*surface);

                const real invdz2  = real(1.0)/(dz*dz);

                real fmod = (real(9.0)*c9*sinvdz9-
                             real(7.0)*c7*sinvdz7+
                             real(3.0)*c3*sinvdz3-
                             real(3.0)*c3hydro*sinvdz3)*invdz2; 
                
                /*
                printf("rad: %f, sigmaS: %f," 
                       "eps: %f, epsS: %f,"
                       "theta1: %f, theta2: %f, theta3: %f,"
                       "thetaS: %f, thetaP: %f,"
                       "chiSurface: %f, chiParticle: %f," 
                       "surfacePos: %f\n",
                        radius,sigmaSurface,
                        epsilon,epsilonSurface,
                        theta1,theta2,theta3,
                        thetaS,thetaP,
                        chiSurface,surface,
                        surfacePosition);
    	        */
                
                return {real(0),
                        real(0),
                        -fmod*(surfacePosition-pos.z)};
            }
    	    
            __device__ __forceinline__ real energy(const real4& pos,const real& radius,const real& epsilon, real& surface){
    	        
                const real dz = abs(surfacePosition-pos.z);
                
                const real sigma     = (sigmaSurface+radius)/real(2.0);
                const real epsilonLB = sqrt(epsilonSurface*epsilon);

                const real sigma3  = sigma*sigma*sigma;
                const real sinvdz3 = sigma3/(dz*dz*dz);
                const real sinvdz7 = sinvdz3*sinvdz3*(sigma/dz);
                const real sinvdz9 = sinvdz3*sinvdz3*sinvdz3;

                const real ccommon = real(M_PI)*rho*sigma3*epsilonLB;
                const real c9      = ccommon*theta1;
                const real c7      = ccommon*theta2;
                const real c3      = ccommon*theta3;
                const real c3hydro = ccommon*(thetaS*chiSurface+thetaP*surface);
                
                real e = (c9*sinvdz9-
                          c7*sinvdz7+
                          c3*sinvdz3-
                          c3hydro*sinvdz3); 

                return e;
            }
            
            __device__ __forceinline__ real3 virial(const real4& pos,const real& radius,const real& epsilon, real& surface){
    	        
                return make_real3(0);
            }
    	    
    	    std::tuple<const real4 *,const real *,const real *,const real *> getArrays(ParticleData *pd){
    	    	
                real4* pos = pd->getPos(access::location::gpu, access::mode::read).raw();

                real* radius;
                if(useInnerRadius){
                    radius = pd->getInnerRadius(access::location::gpu, access::mode::readwrite).raw();     
                } else {
                    radius = pd->getRadius(access::location::gpu, access::mode::readwrite).raw();     
                }
                
                real*  epsilon = pd->getEpsilon(access::location::gpu, access::mode::read).raw();
                real*  surface = pd->getSurface(access::location::gpu, access::mode::read).raw();
    	    	
                return std::make_tuple(pos,radius,epsilon,surface);
    	    }
    };
    
}}}}

#endif
