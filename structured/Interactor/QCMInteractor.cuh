#ifndef __QCM_INTERACTOR__
#define __QCM_INTERACTOR__

namespace uammd{
namespace structured{

    namespace Interactor{
    namespace QCM_ns{

    /* QCM interactor (sinusoidal forcing of the wall and vertical potential) */
    struct QCMForcing : public ParameterUpdatable {
            
        std::shared_ptr<System> sys;

        real time;
            
        real amplitude;
        real frequency;

        struct Parameters{
            real amplitude;
            real frequency;
        };
            
        Parameters inputFileToParam(InputFile& in){

            Parameters param;

            in.getOption("QCMamplitude",InputFile::Required)
                >>param.amplitude;
            in.getOption("QCMfrequency",InputFile::Required)
                >>param.frequency;

            return param;

        }

        QCMForcing(std::shared_ptr<System> sys,
                   Parameters param):sys(sys),
                                     amplitude(param.amplitude), frequency(param.frequency) {
                
                this->sys->log<System::MESSAGE>("[QCMForcing] "
                                                 "Parameter QCMaplitude %f",amplitude);
                this->sys->log<System::MESSAGE>("[QCMForcing] "
                                                 "Parameter frequency %f",frequency);

                time = 0; /* Local clock variable for the interactor */
            }
        
        QCMForcing(std::shared_ptr<System> sys,
                   InputFile& in):QCMForcing(sys,inputFileToParam(in)){}
    
        /* Forces acting on QCM wall particles */
        __device__ real3 force(real4 position) {
            return make_real3(amplitude*sin(real(2*M_PI)*frequency*time),
                              0.0,
                              0.0);
            /* In addition to these forces, QCM particles are also usually bonded to
               each other and/or to fixed points in space. */
        }
    
        auto getArrays(ParticleData *pd) {
            auto pos = pd->getPos(access::location::gpu, access::mode::read);
            return std::make_tuple(pos.raw());
        }
    
        /* Update interactor time */
        virtual void updateSimulationTime(real newTime) override {
            time = newTime;
        }
        
        template<class UNITS>
        void applyUnits(){

            frequency = frequency/UNITS::TO_INTERNAL_TIME;

            this->sys->log<System::MESSAGE>("[QCMForcing] "
                                            "Parameter frequency (after units) %f",frequency);

        }


    };
    
    }

    using QCMInteractor = ExternalForces<QCM_ns::QCMForcing>;

}}}

#endif
