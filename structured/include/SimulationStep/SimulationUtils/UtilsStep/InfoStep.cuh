#ifndef __INFO_STEP__
#define __INFO_STEP__

namespace uammd{
namespace structured{
namespace SimulationStep{
namespace SimulationUtils{

class InfoStep: public SimulationStepBase {

        ullint safeSum(ullint a, ullint b){
            ullint sum;
            if(a > std::numeric_limits<ullint>::max() - b){
                sum = std::numeric_limits<ullint>::max();
            } else {
                sum = a + b;
            }
            return sum;
        }

        Timer tim;

        ullint totalSteps;

    public:

        InfoStep(std::shared_ptr<ParticleGroup>              pg,
                 std::shared_ptr<IntegratorManager>  integrator,
                 std::shared_ptr<ForceField>                 ff,
                 DataEntry& data,
                 std::string name):SimulationStepBase(pg,integrator,ff,data,name){
        }

        void init(cudaStream_t st) override{

            auto integratorsSteps = integrator->getSortedIntegratorSteps();

            totalSteps = 0;
            for(auto intStep : integratorsSteps){
                totalSteps = safeSum(totalSteps, intStep.steps);
            }

            if(totalSteps == std::numeric_limits<ullint>::max()){
                //We consider to be infinite
                System::log<System::WARNING>("[InfoStep] (%s) Total steps: INFINITE", name.c_str());
            } else {
                System::log<System::MESSAGE>("[InfoStep] (%s) Total steps: %llu",name.c_str(),totalSteps);
            }


            tim = Timer();
            tim.tic();
        }

        void applyStep(ullint step, cudaStream_t st) override{

            if(step != 0){
                real time = tim.toc();

                std::chrono::seconds rtime_s(ullint((totalSteps-step)/(real(this->intervalStep)/time)));

                //rtime_s to days, hours, minutes, seconds
                uint days    = std::chrono::duration_cast<std::chrono::hours>(rtime_s).count()/24;
                uint hours   = std::chrono::duration_cast<std::chrono::hours>(rtime_s).count()%24;
                uint minutes = std::chrono::duration_cast<std::chrono::minutes>(rtime_s).count()%60;
                uint seconds = std::chrono::duration_cast<std::chrono::seconds>(rtime_s).count()%60;

                if(totalSteps == std::numeric_limits<ullint>::max()){
                    System::log<System::MESSAGE>("[InfoStep] (%s) Step %llu, ETA: unknown, mean FPS: %0.2f",
                                                  this->name.c_str(),step,real(this->intervalStep)/time);
                }else{
                    System::log<System::MESSAGE>("[InfoStep] (%s) Step %llu, ETA: d:%u h:%u m:%u s:%u, mean FPS: %0.2f",
                                                  this->name.c_str(),step,days,hours,minutes,seconds,real(this->intervalStep)/time);
                }
            }

            tim.tic();
        }

};

}}}}

#endif
