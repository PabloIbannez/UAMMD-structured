#ifndef __SPHERICAL_SHELL__
#define __SPHERICAL_SHELL__
    
namespace uammd{
namespace structured{ 
namespace Potentials{
namespace Bounds{
    
    class SphericalShell: public ParameterUpdatable{
    	
        private:
            
            real epsilonShell;
            real sigmaShell;
            real Ashell;
            real Bshell;

            real3 shellCenter;
            real  shellRadius;
    	        
        public:

            struct Parameters{
            
                real epsilonShell;
                real sigmaShell;
                real Ashell;
                real Bshell;
    	        
                real3 shellCenter;
                real  shellRadius;
            };
    	    
    	    SphericalShell(Parameters par):epsilonShell(par.epsilonShell),
                                           sigmaShell(par.sigmaShell),
                                           Ashell(par.Ashell),
                                           Bshell(par.Bshell),
                                           shellCenter(par.shellCenter),
                                           shellRadius(par.shellRadius){}     
    	    
    	    __device__ __forceinline__ real3 force(const real4 &pos){
                
                const real3 rij = shellCenter-make_real3(pos);
                
                real  r2 = abs(sqrt(dot(rij,rij))-shellRadius);
                      r2 = r2*r2;
                
                real3 f12 = Ashell*CommonPotentials::Steric::Steric::force<12>(rij,r2,epsilonShell,sigmaShell);
                real3 f6  = Bshell*CommonPotentials::Steric::Steric::force<6>(rij,r2,epsilonShell,sigmaShell);
                
                real3 f = -(f12 + f6);

                return f;
            }
    	    
            __device__ __forceinline__ real  energy(const real4 &pos){
                
                const real3 rij = shellCenter-make_real3(pos);
                
                real  r2 = abs(sqrt(dot(rij,rij))-shellRadius);
                      r2 = r2*r2;
                
                real e = Ashell*CommonPotentials::Steric::Steric::energy<12>(rij,r2,epsilonShell,sigmaShell) + 
                         Bshell*CommonPotentials::Steric::Steric::energy<6>(rij,r2,epsilonShell,sigmaShell);

                return real(2.0)*e;
            
            }
            __device__ __forceinline__ real3 virial(const real4 &pos){return make_real3(0);}
    	    
    	    std::tuple<const real4 *> getArrays(ParticleData *pd){
    	    	
                auto pos = pd->getPos(access::location::gpu, access::mode::read);
    	    	
                return std::make_tuple(pos.raw());
    	    }

            void setShellRadius(real newShellRadius){
                shellRadius = newShellRadius;
            }
            
            real getShellRadius(){return shellRadius;}
    };
    
}}}}

#endif
