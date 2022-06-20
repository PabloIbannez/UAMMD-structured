#ifndef __COMMON_POTENTIALS__
#define __COMMON_POTENTIALS__

namespace uammd{
namespace structured{ 
namespace Potentials{
namespace CommonPotentials{

    struct Harmonic{

        static inline __device__ real3 force(const real3& rij, const real& r2,
                                             const real& K,const real& r0){

            const real invr = rsqrt(r2);
            
            const real fmod = K*(real(1.0)-r0*invr);

            return fmod*rij;
        }
        
        static inline __device__ real energy(const real3& rij, const real& r2,
                                             const real& K,const real& r0){
            
            const real dr = sqrt(r2)-r0;
            
            return real(0.5)*K*dr*dr;
            
        }
    };
    
    struct HarmonicAnisotropic{

        static inline __device__ real3 force(const real3& rij,const real3& K,const real3& r0){

            return {K.x*(rij.x-r0.x),
                    K.y*(rij.y-r0.y),
                    K.z*(rij.z-r0.z)};
        }
        
        static inline __device__ real energy(const real3& rij,const real3& K,const real3& r0){
            
            return real(0.5)*K.x*(rij.x-r0.x)*(rij.x-r0.x)+
                   real(0.5)*K.y*(rij.y-r0.y)*(rij.y-r0.y)+
                   real(0.5)*K.z*(rij.z-r0.z)*(rij.z-r0.z);
            
        }
    };
    
    struct Morse{

        static inline __device__ real3 force(const real3& rij, const real& r2,
                                             const real& e,const real& r0,const real& D){

            const real  r = sqrt(r2);
            const real dr = r-r0;
            const real factor = exp(-dr/D);
            
            const real fmod = real(2.0)*(e/D)*(real(1.0)-factor)*factor;

            return fmod*rij/r;
        }
        
        static inline __device__ real energy(const real3& rij, const real& r2,
                                             const real& e,const real& r0,const real& D){
            
            const real dr = sqrt(r2)-r0;
            const real factor = real(1.0)-exp(-dr/D);
            
            return e*(factor*factor-real(1.0));
            
        }
    };

    struct GaussianWell{
        
        static inline __device__ real3 force(const real3& rij, const real& r2,
                                             const real& e,const real& r0,const real& D){

            const real r = sqrt(r2);
            const real dr = r-r0;
            
            const real fmod = (e/D)*(real(1.0)-r0/r)*exp(-dr*dr/(real(2.0)*D));

            return fmod*rij;
        }
        
        static inline __device__ real energy(const real3& rij, const real& r2,
                                             const real& e,const real& r0,const real& D){
            
            const real dr = sqrt(r2)-r0;
            
            return -e*exp(-dr*dr/(real(2.0)*D))/real(2.0);
            
        }

    };

    namespace LennardJones{

        //U(r)=4e((s/r)^12-(s/r)^6)
        struct Type1{
        
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                
                const real invr2   = real(1.0)/r2;
                const real sinvr2  = sigma*sigma*invr2;
                const real sinvr6  = sinvr2*sinvr2*sinvr2;
                const real sinvr12 = sinvr6*sinvr6;

                real fmod = -real(24.0)*epsilon*(real(2.0)*sinvr12-sinvr6)*invr2;
                
                return fmod*rij;
            }
            
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return computeVirial(rij,force(rij,r2,epsilon,sigma));
            }
            
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return computeStress(rij,force(rij,r2,epsilon,sigma));
            }
            
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                
                const real sinvr2  = sigma*sigma/r2;
                const real sinvr6  = sinvr2*sinvr2*sinvr2;
                const real sinvr12 = sinvr6*sinvr6;

                real e = real(4.0)*epsilon*(sinvr12-sinvr6);
                
                return e/real(2.0);
            }

        };
        
        //U(r)=e((s/r)^12-2*(s/r)^6)
        struct Type2{
        
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                
                const real invr2   = real(1.0)/r2;
                const real sinvr2  = sigma*sigma*invr2;
                const real sinvr6  = sinvr2*sinvr2*sinvr2;
                const real sinvr12 = sinvr6*sinvr6;

                real fmod = -real(12.0)*epsilon*(sinvr12-sinvr6)*invr2;

                return fmod*rij;
            }
            
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return computeVirial(rij,force(rij,r2,epsilon,sigma));
            }
            
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return computeStress(rij,force(rij,r2,epsilon,sigma));
            }
            
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                
                const real sinvr2  = sigma*sigma/r2;
                const real sinvr6  = sinvr2*sinvr2*sinvr2;
                const real sinvr12 = sinvr6*sinvr6;

                real e = epsilon*(sinvr12-real(2.0)*sinvr6);
                
                return e/real(2.0);
            }

        };
        
        //U(r)=e(5(s/r)^12-6(s/r)^10)
        struct Type3{
        
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                
                const real invr2   = real(1.0)/r2;
                const real sinvr2  = sigma*sigma*invr2;
                const real sinvr4  = sinvr2*sinvr2;
                const real sinvr6  = sinvr4*sinvr2;
                const real sinvr10 = sinvr6*sinvr4;
                const real sinvr12 = sinvr6*sinvr6;

                real fmod = -real(60.0)*epsilon*(sinvr12-sinvr10)*invr2;
                
                return fmod*rij;
            }
            
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return computeVirial(rij,force(rij,r2,epsilon,sigma));
            }
            
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return computeStress(rij,force(rij,r2,epsilon,sigma));
            }
            
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                
                const real invr2   = real(1.0)/r2;
                const real sinvr2  = sigma*sigma*invr2;
                const real sinvr4  = sinvr2*sinvr2;
                const real sinvr6  = sinvr4*sinvr2;
                const real sinvr10 = sinvr6*sinvr4;
                const real sinvr12 = sinvr6*sinvr6;

                real e = epsilon*(real(5.0)*sinvr12-real(6.0)*sinvr10);
                
                return e/real(2.0);
            }

        };
        
        //U(r)=e(13(s/r)^12-18(s/r)^10+4(s/r)^6)
        struct KaranicolasBrooks{
        
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                
                const real invr2   = real(1.0)/r2;
                const real sinvr2  = sigma*sigma*invr2;
                const real sinvr4  = sinvr2*sinvr2;
                const real sinvr6  = sinvr4*sinvr2;
                const real sinvr10 = sinvr6*sinvr4;
                const real sinvr12 = sinvr6*sinvr6;

                real fmod = -epsilon*(real(156.0)*sinvr12-real(180.0)*sinvr10+real(24.0)*sinvr6)*invr2;
                
                return fmod*rij;
            }
            
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return computeVirial(rij,force(rij,r2,epsilon,sigma));
            }
            
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return computeStress(rij,force(rij,r2,epsilon,sigma));
            }
            
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                
                const real invr2   = real(1.0)/r2;
                const real sinvr2  = sigma*sigma*invr2;
                const real sinvr4  = sinvr2*sinvr2;
                const real sinvr6  = sinvr4*sinvr2;
                const real sinvr10 = sinvr6*sinvr4;
                const real sinvr12 = sinvr6*sinvr6;

                real e = epsilon*(real(13.0)*sinvr12-real(18.0)*sinvr10+real(4.0)*sinvr6);
                
                return e/real(2.0);
            }

        };

    }

    namespace Steric{
    
        struct Steric{

            //Force
            template<int power>
            static inline __device__ real3 force(const real3& rij, const real& r2, 
                                          const real& epsilon,const real& sigma){
                static_assert(power==6 or power==12,"Steric power has to be 6 or 12");
                return make_real3(0);
            }
            
            //Virial
            template<int power>
            static inline __device__ real virial(const real3& rij, const real& r2, 
                                                 const real& epsilon,const real& sigma){
                static_assert(power==6 or power==12,"Steric power has to be 6 or 12");
                return real(0);
            }
            
            //Stress
            template<int power>
            static inline __device__ tensor3 stress(const real3& rij, const real& r2, 
                                             const real& epsilon,const real& sigma){
                static_assert(power==6 or power==12,"Steric power has to be 6 or 12");
                return tensor3(0);
            }
            
            
            //Energy
            template<int power>
            static inline __device__ real energy(const real3& rij, const real& r2, 
                                          const real& epsilon,const real& sigma){
                static_assert(power==6 or power==12,"Steric power has to be 6 or 12");
                return real(0);
            }
            

        };
            
        template<>
        inline __device__ real3 Steric::force<6>(const real3& rij, const real& r2, 
                                         const real& epsilon,const real& sigma){
            
            const real  invr2  = real(1.0)/r2;
            const real sinvr2  = sigma*sigma*invr2;
            const real sinvr6  = sinvr2*sinvr2*sinvr2;

            const real fmod = -real(6.0)*epsilon*sinvr6*invr2;

            return fmod*rij;
        }
        
        template<>
        inline __device__ real3 Steric::force<12>(const real3& rij, const real& r2, 
                                          const real& epsilon,const real& sigma){
            
            const real  invr2  = real(1.0)/r2;
            const real sinvr2  = sigma*sigma*invr2;
            const real sinvr6  = sinvr2*sinvr2*sinvr2;
            const real sinvr12 = sinvr6*sinvr6;

            const real fmod = -real(12.0)*epsilon*sinvr12*invr2;

            return fmod*rij;
        }
        
        
        template<>
        inline __device__ real Steric::virial<6>(const real3& rij, const real& r2, 
                                                 const real& epsilon,const real& sigma){
            
            return computeVirial(rij,force<6>(rij,r2,epsilon,sigma));
        }
        
        template<>
        inline __device__ tensor3 Steric::stress<6>(const real3& rij, const real& r2, 
                                            const real& epsilon,const real& sigma){
            
            return computeStress(rij,force<6>(rij,r2,epsilon,sigma));
        }
        
        template<>
        inline __device__ real Steric::virial<12>(const real3& rij, const real& r2, 
                                                  const real& epsilon,const real& sigma){
            
            return computeVirial(rij,force<12>(rij,r2,epsilon,sigma));
        }
        
        template<>
        inline __device__ tensor3 Steric::stress<12>(const real3& rij, const real& r2, 
                                             const real& epsilon,const real& sigma){
            
            return computeStress(rij,force<12>(rij,r2,epsilon,sigma));
        }
        
        template<>
        inline __device__ real Steric::energy<6>(const real3& rij, const real& r2, 
                                         const real& epsilon,const real& sigma){
            
            const real  invr2  = real(1.0)/r2;
            const real sinvr2  = sigma*sigma*invr2;
            const real sinvr6  = sinvr2*sinvr2*sinvr2;

            const real e = epsilon*sinvr6;

            return e/real(2.0);
        }
        
        template<>
        inline __device__ real Steric::energy<12>(const real3& rij, const real& r2, 
                                          const real& epsilon,const real& sigma){
            
            const real  invr2  = real(1.0)/r2;
            const real sinvr2  = sigma*sigma*invr2;
            const real sinvr6  = sinvr2*sinvr2*sinvr2;
            const real sinvr12 = sinvr6*sinvr6;

            const real e = epsilon*sinvr12;

            return e/real(2.0);
        }
        
        struct Steric6{
            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                
                return Steric::force<6>(rij,r2,epsilon,sigma);
            }
            
            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return Steric::virial<6>(rij,r2,epsilon,sigma);
            }
            
            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return Steric::stress<6>(rij,r2,epsilon,sigma);
            }
            
            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return Steric::energy<6>(rij,r2,epsilon,sigma);
            }
        };
        
        struct Steric12{
            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                
                return Steric::force<12>(rij,r2,epsilon,sigma);
            }
            
            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return Steric::virial<12>(rij,r2,epsilon,sigma);
            }
            
            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return Steric::stress<12>(rij,r2,epsilon,sigma);
            }
            
            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return Steric::energy<12>(rij,r2,epsilon,sigma);
            }
        };


    }
    
    namespace WCA{
        
        struct Type1{
            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                
                //r2 > (simga*2^(1/6))^2
                if(r2 > sigma*sigma*real(1.259921)) return make_real3(0);

                return LennardJones::Type1::force(rij,r2,epsilon,sigma);
            }
            
            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return real(0);
            }
            
            //Stress
            static inline __device__ tensor3 Stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return tensor3(0);
            }
            
            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                //r2 > (simga*2^(1/6))^2
                if(r2 > sigma*sigma*real(1.259921)) return real(0);

                return LennardJones::Type1::energy(rij,r2,epsilon,sigma)+epsilon;
            }
        };
        
        struct Type2{
            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                
                if(r2 > sigma*sigma) return make_real3(0);

                return LennardJones::Type2::force(rij,r2,epsilon,sigma);
            }
            
            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return real(0);
            }
            
            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return tensor3(0);
            }
            
            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                if(r2 > sigma*sigma) return real(0);

                return LennardJones::Type2::energy(rij,r2,epsilon,sigma)+epsilon;
            }
        };
        
        struct Type3{
            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                
                if(r2 > sigma*sigma) return make_real3(0);

                return LennardJones::Type3::force(rij,r2,epsilon,sigma);
            }
            
            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return real(0);
            }
            
            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return tensor3(0);
            }
            
            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                if(r2 > sigma*sigma) return real(0);

                return LennardJones::Type3::energy(rij,r2,epsilon,sigma)+epsilon;
            }
        };
    }
    
    namespace GeneralLennardJones{
        
        struct Type1{
            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                
                if(epsilon >= real(0.0)){
                    return WCA::Type1::force(rij,r2,epsilon,sigma);
                }
                    
                return LennardJones::Type1::force(rij,r2,fabs(epsilon),sigma);

            }
            
            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return real(0);
            }
            
            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return tensor3(0);
            }
            
            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                if(epsilon >= real(0.0)){
                    return WCA::Type1::energy(rij,r2,epsilon,sigma);
                }
                    
                return LennardJones::Type1::energy(rij,r2,fabs(epsilon),sigma);
            }
        };
        
        struct Type2{
            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                
                if(epsilon >= real(0.0)){
                    return WCA::Type2::force(rij,r2,epsilon,sigma);
                }
                    
                return LennardJones::Type2::force(rij,r2,fabs(epsilon),sigma);

            }
            
            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return real(0);
            }
            
            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return tensor3(0);
            }
            
            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                if(epsilon >= real(0.0)){
                    return WCA::Type2::energy(rij,r2,epsilon,sigma);
                }
                    
                return LennardJones::Type2::energy(rij,r2,fabs(epsilon),sigma);
            }
        };
        
        struct Type3{
            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                
                if(epsilon >= real(0.0)){
                    return WCA::Type3::force(rij,r2,epsilon,sigma);
                }
                    
                return LennardJones::Type3::force(rij,r2,fabs(epsilon),sigma);

            }
            
            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){
                return real(0);
            }
            
            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma){
                return tensor3(0);
            }
            
            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma){

                if(epsilon >= real(0.0)){
                    return WCA::Type3::energy(rij,r2,epsilon,sigma);
                }
                    
                return LennardJones::Type3::energy(rij,r2,fabs(epsilon),sigma);
            }
        };
    }
        
    namespace ModifiedLennardJones{

        struct Gaussian{
            
            //Force
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma,const real& D){
                
                return WCA::Type2::force(rij,r2,epsilon,sigma)+
                       GaussianWell::force(rij,r2,epsilon,sigma,D);
            }
            
            //Virial
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma,const real& D){
                return real(0);
            }
            
            //Stress
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& epsilon,const real& sigma,const real& D){
                return tensor3(0);
            }
            
            //Energy
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& epsilon,const real& sigma,const real& D){
                return WCA::Type2::energy(rij,r2,epsilon,sigma)+
                       GaussianWell::energy(rij,r2,epsilon,sigma,D);
            }

        };
    }
        
    struct NonPolar{

        //Force
        static inline __device__ real3 force(const real3& rij, const real& r2, 
                                             const real& epsilon,const real& sigma,const real& zeroEnergy){
            
            if(epsilon == real(0.0) ){
                return CommonPotentials::Steric::Steric::force<12>(rij,r2,zeroEnergy,sigma);
            }
        
            real3 f = CommonPotentials::LennardJones::Type1::force(rij,r2,abs(epsilon),sigma);
            //(2^(1/6)*sigma)^2
            const real r02 = (real(1.259921)*sigma*sigma);

            if(epsilon > real(0.0) and r2>=r02){
                f=-f;
            }

            return f;
        }
        
        //Virial
        static inline __device__ real virial(const real3& rij, const real& r2, 
                                             const real& epsilon,const real& sigma,const real& zeroEnergy){
            return real(0);
        }
        
        //Stress
        static inline __device__ tensor3 stress(const real3& rij, const real& r2, 
                                                const real& epsilon,const real& sigma,const real& zeroEnergy){
            return tensor3(0);
        }
        
        //Energy
        static inline __device__ real energy(const real3& rij, const real& r2, 
                                             const real& epsilon,const real& sigma,const real& zeroEnergy){
            
            if(epsilon == real(0.0) ){
                return CommonPotentials::Steric::Steric::energy<12>(rij,r2,zeroEnergy,sigma);
            }
            
            real e  = CommonPotentials::LennardJones::Type1::energy(rij,r2,abs(epsilon),sigma);
        
            //(2^(1/6)*sigma)^2
            const real r02 = (real(1.259921)*sigma*sigma);

            if(epsilon > real(0.0)){
                   
                if(r2<r02){
                    e+=epsilon; //Implicit 1/2
                } else {
                    e=-e;
                }
            }
        
            return e;
        }

    };

    namespace DebyeHuckel{
        
        struct DebyeHuckel{

            //Force
            template<typename Units_>
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& chgProduct, 
                                                 const real& dielectricConstant,const real& debyeLenght){
                
                const real      r  = sqrt(r2);
                const real  invr2  = real(1.0)/r2;

                const real efactor = Units_::ELECOEF*chgProduct/dielectricConstant; //ELECOEF = 1/(4*pi*e_0)

                //printf("%f\n",efactor);

                //printf("e %f %f %f %f %f %f\n",efactor,pCG::ELECOEF,dielectricConstant,rij.x,rij.y,rij.z);

                real fmod = -efactor*exp(-r/debyeLenght)*invr2*(real(1.0)/debyeLenght+real(1.0)/r);

                //printf("%f\n",fmod);

                return fmod*rij;
            }
            
            //Virial
            template<typename Units_>
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& chgProduct, 
                                                 const real& dielectricConstant,const real& debyeLenght){
                return real(0);
            }
            
            //Stress
            template<typename Units_>
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& chgProduct, 
                                                    const real& dielectricConstant,const real& debyeLenght){
                return tensor3(0);
            }
            
            //Energy
            template<typename Units_>
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& chgProduct, 
                                                 const real& dielectricConstant,const real& debyeLenght){

                const real      r  = sqrt(r2);

                const real efactor = Units_::ELECOEF*chgProduct/dielectricConstant; //ELECOEF = 1/(4*pi*e_0)

                //printf("e %f %f %f %f %f %f\n",efactor,pCG::ELECOEF,dielectricConstant,rij.x,rij.y,rij.z);

                real e = efactor*exp(-r/debyeLenght)/r;

                return e/real(2.0);
            }
        };
         
        struct DebyeHuckelSpheres{

            //Force
            template<typename Units_>
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& chgProduct, 
                                                 const real& radius1, const real& radius2,
                                                 const real& dielectricConstant,const real& debyeLenght){
                
                const real      r  = sqrt(r2);
                const real  invr2  = real(1.0)/r2;

                const real efactor = (Units_::ELECOEF*chgProduct/dielectricConstant)/
                                     ((real(1.0)+radius1/debyeLenght)*(real(1.0)+radius2/debyeLenght)); //ELECOEF = 1/(4*pi*e_0)

                real fmod = -efactor*exp(-(r-radius1-radius2)/debyeLenght)*invr2*(real(1.0)/debyeLenght+real(1.0)/r);

                return fmod*rij;
            }
            
            //Virial
            template<typename Units_>
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& chgProduct, 
                                                 const real& radius1, const real& radius2,
                                                 const real& dielectricConstant,const real& debyeLenght){
                return real(0);
            }
            
            //Stress
            template<typename Units_>
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& chgProduct, 
                                                    const real& radius1, const real& radius2,
                                                    const real& dielectricConstant,const real& debyeLenght){
                return tensor3(0);
            }
            
            //Energy
            template<typename Units_>
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& chgProduct, 
                                                 const real& radius1, const real& radius2,
                                                 const real& dielectricConstant,const real& debyeLenght){

                const real      r  = sqrt(r2);

                const real efactor = (Units_::ELECOEF*chgProduct/dielectricConstant)/
                                     ((real(1.0)+radius1/debyeLenght)*(real(1.0)+radius2/debyeLenght)); //ELECOEF = 1/(4*pi*e_0)

                real e = efactor*exp(-(r-radius1-radius2)/debyeLenght)/r;

                return e/real(2.0);
            }
        };
        
        inline __device__ real distanceDependentDielectric(const real& r,
                                                           const real& A, const real& B,
                                                           const real& k, const real lmd){
            return A+B/(real(1.0)+k*exp(-lmd*B*r));
        }
        
        inline __device__ real distanceDependentDielectricDerivate(const real& r,
                                                                   const real& A, const real& B,
                                                                   const real& k, const real lmd){

            real den = exp(lmd*B*r)+k;
            return B*B*k*lmd*exp(lmd*B*r)/(den*den);
        }
        
        struct DebyeHuckelDistanceDependentDielectric{

            //Force
            template<typename Units_>
            static inline __device__ real3 force(const real3& rij, const real& r2,
                                                 const real& chgProduct,const real& debyeLenght,
                                                 const real& A, const real& B,
                                                 const real& k, const real lmd){
                
                const real      r  = sqrt(r2);
                const real  invr2  = real(1.0)/r2;

                const real dielectricConstant = distanceDependentDielectric(r,A,B,k,lmd);

                const real efactor = Units_::ELECOEF*chgProduct/dielectricConstant; 

                real fmod = -efactor*exp(-r/debyeLenght)*invr2*(real(1.0)/debyeLenght+real(1.0)/r);
                     fmod = fmod - (efactor*exp(-r/debyeLenght)/r)*distanceDependentDielectricDerivate(r,A,B,k,lmd)/(r*dielectricConstant);

                return fmod*rij;
            }
            
            //Virial
            template<typename Units_>
            static inline __device__ real virial(const real3& rij, const real& r2,
                                                 const real& chgProduct,const real& debyeLenght,
                                                 const real& A, const real& B,
                                                 const real& k, const real lmd){

                return real(0);
            }
            
            //Stress
            template<typename Units_>
            static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                    const real& chgProduct,const real& debyeLenght,
                                                    const real& A, const real& B,
                                                    const real& k, const real lmd){

                return tensor3(0);
            }
            
            //Energy
            template<typename Units_>
            static inline __device__ real energy(const real3& rij, const real& r2,
                                                 const real& chgProduct,const real& debyeLenght,
                                                 const real& A, const real& B,
                                                 const real& k, const real lmd){

                const real      r  = sqrt(r2);

                const real dielectricConstant = distanceDependentDielectric(r,A,B,k,lmd);
                
                const real efactor = Units_::ELECOEF*chgProduct/dielectricConstant; //ELECOEF = 1/(4*pi*e_0)

                //printf("e %f %f %f %f %f %f\n",efactor,pCG::ELECOEF,dielectricConstant,rij.x,rij.y,rij.z);

                real e = efactor*exp(-r/debyeLenght)/r;

                return e/real(2.0);
            }
        };

    }

    struct Contact{
        
        //Force
        static inline __device__ real3 force(const real3& rij, const real& r2,
                                             const real& epsilon,const real& sigma,
                                             const real& n,const real& r0,
                                             const real& zeroEnergy){

            real eps;
            if(epsilon==real(0.0)){
                eps=  abs(zeroEnergy);
            } else {
                eps = abs(epsilon);
            }

            real3 fs = Steric::Steric::force<12>(rij,r2,eps,sigma);

            real r     = sqrt(r2);
            real sech2 = real(1.0)/cosh(n*(r0-r));
                 sech2 = sech2*sech2;
            
            real3 f  = -real(0.5)*epsilon*n*sech2*rij/r;
        
            return fs+f;
        }
        
        //Virial
        static inline __device__ real virial(const real3& rij, const real& r2,
                                             const real& epsilon,const real& sigma,
                                             const real& n,const real& r0,
                                             const real& zeroEnergy){

            return real(0);
        }
        
        //Stress
        static inline __device__ tensor3 stress(const real3& rij, const real& r2,
                                                const real& epsilon,const real& sigma,
                                                const real& n,const real& r0,
                                                const real& zeroEnergy){

            return tensor3(0);
        }
        
        //Energy
        static inline __device__ real energy(const real3& rij, const real& r2,
                                             const real& epsilon,const real& sigma,
                                             const real& n,const real& r0,
                                             const real& zeroEnergy){

            real eps;
            if(epsilon==real(0.0)){
                eps=zeroEnergy;
            } else {
                eps = epsilon;
            }

            real es = Steric::Steric::energy<12>(rij,r2,eps,sigma);
            real e  = real(0.5)*epsilon*(real(1.0)+tanh(n*(r0-sqrt(r2))));

            return es+e/real(2.0);
        }
    };
    
}}}}

#endif
