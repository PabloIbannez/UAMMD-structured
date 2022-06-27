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
                
                const real3 dr = shellCenter-make_real3(pos);
                
                const real r  = sqrt(dot(dr,dr));
                      real r2 = r-shellRadius;
                           r2 = r2*r2;
                
                const real sinvr2  = sigmaShell*sigmaShell/r2;
                const real sinvr6  = sinvr2*sinvr2*sinvr2;
                const real sinvr12 = sinvr6*sinvr6;
                
                real fmod = -epsilonShell*real(6.0)*(real(2.0)*Ashell*sinvr12+Bshell*sinvr6)/abs(r-shellRadius);

                real3 f = -fmod*(dr/r);
                
                //real  r2 = sqrt(dot(rij,rij))-shellRadius;
                //      r2 = r2*r2;
                //
                //real3 f12 = Ashell*CommonPotentials::Steric::Steric::force<12>(rij,r2,epsilonShell,sigmaShell);
                //real3 f6  = Bshell*CommonPotentials::Steric::Steric::force<6>(rij,r2,epsilonShell,sigmaShell);

                //real3 f = -(f12 + f6);
                //
                //if(sqrt(dot(f,f))>real(0.1)){
                //    printf("A:%f B:%f f: %f f6:%f f6x:%f f6y:%f f6z:%f f12:%f f12x:%f f12y:%f f12z:%f r:%f\n",Ashell,Bshell,sqrt(dot(f,f)),sqrt(dot(f6,f6)),f6.x,f6.y,f6.z,
                //                                                                                                   sqrt(dot(f12,f12)),f12.x,f12.y,f12.z,
                //                                                                                                   sqrt(dot(rij,rij)));
                //}

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
            __device__ __forceinline__ real virial(const real4 &pos){
                
                const real3 dr = shellCenter-make_real3(pos);
                
                const real r  = sqrt(dot(dr,dr));
                      real r2 = r-shellRadius;
                           r2 = r2*r2;
                
                const real sinvr2  = sigmaShell*sigmaShell/r2;
                const real sinvr6  = sinvr2*sinvr2*sinvr2;
                const real sinvr12 = sinvr6*sinvr6;
                
                real fmod = -epsilonShell*real(6.0)*(real(2.0)*Ashell*sinvr12+Bshell*sinvr6)/abs(r-shellRadius);

                real3 f = -fmod*(dr/r);

                return -dot(f,dr/r);
            }
            
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

            void setShellRadius(real newShellRadius){
                shellRadius = newShellRadius;
            }
            
            real getShellRadius(){return shellRadius;}
    };
    
}}}}

#endif
