#ifndef __GAMD_INTERACTOR__
#define __GAMD_INTERACTOR__

namespace uammd{
namespace structured{

    namespace Interactor{

    namespace GAMD_ns{
    
        template<typename Tin1,typename Tin2,typename Tout>
        struct multiplies
        {
            typedef Tin1 first_argument_type;
            typedef Tin2 second_argument_type;
            typedef Tout result_type;
        
            __thrust_exec_check_disable__
                __host__ __device__ Tout operator()(const Tin1 &lhs, const Tin2 &rhs) const {return lhs * rhs;}
        }; // end multiplies
    
    }
    
    template<class Base>
    class GAMD: public Base{

        private:

            real k;

            real E;
            real k0;
            real Vmax;
            real Vmin;

        public:
            
            GAMD(std::shared_ptr<System>       sys,
                 std::shared_ptr<ParticleData>  pd,
                 std::shared_ptr<ParticleGroup> pg,
                 InputFile& in):Base(sys,pd,pg,in){
                
                    in.getOption("E",InputFile::Required)
                                  >>E;
                    in.getOption("k0",InputFile::Required)
                                  >>k0;
                    in.getOption("Vmax",InputFile::Required)
                                  >>Vmax;
                    in.getOption("Vmin",InputFile::Required)
                                  >>Vmin;
                
                    this->sys->template log<System::MESSAGE>("[GAMD] "
                                                             "Parameter E %f",E);
                    this->sys->template log<System::MESSAGE>("[GAMD] "
                                                              "Parameter k0 %f",k0);
                    this->sys->template log<System::MESSAGE>("[GAMD] "
                                                              "Parameter Vmax %f",Vmax);
                    this->sys->template log<System::MESSAGE>("[GAMD] "
                                                     "Parameter Vmin %f",Vmin);

                    k = k0/(Vmax-Vmin);
            }
            
            ~GAMD(){}
            
            void sum(typename Base::Computables comp,cudaStream_t st) override {

                if(comp.force == true){

                    {
                        auto energy = this->pd->getEnergy(access::location::gpu, access::mode::readwrite);     
                        thrust::fill(thrust::cuda::par.on(st), energy.begin(), energy.end(), real(0));
                    }

                    comp.energy = true;
                    Base::sum(comp,st);

                    real V = Measures::totalPotentialEnergy(this->sys,this->pd,this->pg,st);

                    //Add boost potential
                    if(V<E){

                        real factor = real(1.0) - k*(E-V);
                        auto force  = this->pd->getForce(access::location::gpu, access::mode::readwrite);     

                        thrust::transform(thrust::cuda::par.on(st),
                                          force.begin(),
                                          force.end(),
                                          thrust::make_constant_iterator(factor),
                                          force.begin(),
                                          GAMD_ns::multiplies<real4,real,real4>());
                    }
                }
            }
            
    };

}}}

#endif
