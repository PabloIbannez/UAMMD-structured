#include "Utils/Backup/BackupStep.cuh"

namespace uammd{
namespace structured{
namespace Backup{

    bool BackupStep::isSimulationStateCorrect(){
        //Check if the simulation state is correct
        bool correct = true;

        auto id = pd->getId(access::location::cpu, access::mode::read).raw();

        //Check if some particle position is NaN
        if(pd->isPosAllocated()){
            bool posCorrect = true;
            auto pos = pd->getPos(access::location::cpu, access::mode::read);
            for(int i = 0; i < pd->getNumParticles(); i++){
                if(std::isnan(pos.raw()[i].x) || std::isnan(pos.raw()[i].y) || std::isnan(pos.raw()[i].z)){
                    posCorrect = false;
                    System::log<System::DEBUG>("Particle %i has NaN position", id[i]);
                }
            }
            if(!posCorrect){
                System::log<System::WARNING>("[BackupStep] (%s) Found NaN in particle position",name.c_str());
                correct = false;
            }
        }

        //Check if some particle direction is NaN
        if(pd->isDirAllocated()){
            bool dirCorrect = true;
            auto dir = pd->getDir(access::location::cpu, access::mode::read);
            for(int i = 0; i < pd->getNumParticles(); i++){
                if(std::isnan(dir.raw()[i].x) || std::isnan(dir.raw()[i].y) || std::isnan(dir.raw()[i].z)){
                    dirCorrect = false;
                    System::log<System::DEBUG>("Particle %i has NaN direction", id[i]);
                }
            }
            if(!dirCorrect){
                System::log<System::WARNING>("[BackupStep] (%s) Found NaN in particle direction",name.c_str());
                correct = false;
            }
        }

        //Check if some particle velocity is NaN
        if(pd->isVelAllocated()){
            bool velCorrect = true;
            auto vel = pd->getVel(access::location::cpu, access::mode::read);
            for(int i = 0; i < pd->getNumParticles(); i++){
                if(std::isnan(vel.raw()[i].x) || std::isnan(vel.raw()[i].y) || std::isnan(vel.raw()[i].z)){
                    velCorrect = false;
                    System::log<System::DEBUG>("Particle %i has NaN velocity", id[i]);
                }
            }
            if(!velCorrect){
                System::log<System::WARNING>("[BackupStep] (%s) Found NaN in particle velocity",name.c_str());
                correct = false;
            }
        }

        //Check if some particle force is NaN
        if(pd->isForceAllocated()){
            bool forceCorrect = true;
            auto force = pd->getForce(access::location::cpu, access::mode::read);
            for(int i = 0; i < pd->getNumParticles(); i++){
                if(std::isnan(force.raw()[i].x) || std::isnan(force.raw()[i].y) || std::isnan(force.raw()[i].z)){
                    forceCorrect = false;
                    System::log<System::DEBUG>("Particle %i has NaN force", id[i]);
                }
            }
            if(!forceCorrect){
                System::log<System::WARNING>("[BackupStep] (%s) Found NaN in particle force",name.c_str());
                correct = false;
            }
        }

        //Check if some particle torque is NaN
        if(pd->isTorqueAllocated()){
            bool torqueCorrect = true;
            auto torque = pd->getTorque(access::location::cpu, access::mode::read);
            for(int i = 0; i < pd->getNumParticles(); i++){
                if(std::isnan(torque.raw()[i].x) || std::isnan(torque.raw()[i].y) || std::isnan(torque.raw()[i].z)){
                    torqueCorrect = false;
                    System::log<System::DEBUG>("Particle %i has NaN torque", id[i]);
                }
            }
            if(!torqueCorrect){
                System::log<System::WARNING>("[BackupStep] (%s) Found NaN in particle torque",name.c_str());
                correct = false;
            }
        }

        return correct;
    }

    void BackupStep::writeBackup(std::string backupFileName){

        if(!isSimulationStateCorrect()){
            System::log<System::WARNING>("[BackupStep] (%s) Simulation state is not correct",name.c_str());
            throw std::runtime_error("Simulation state is not correct");
        } else {
            System::log<System::DEBUG>("[BackupStep] (%s) Writing backup to %s",name.c_str(),backupFileName.c_str());
            this->sys->setLastBackupStep(this->gd->getFundamental()->getCurrentStep());

            bool restartedFromBackup = this->sys->getRestartedFromBackup();

            //We update restartedFromBackup to true, so that if the simulation crashes
            this->sys->setRestartedFromBackup(true);

            this->sys->updateInputSystem();
            this->gd->updateInputGlobalData();
            this->pd->updateInputState();
            this->input->write(backupFileName);

            //We back restartedFromBackup to its original value
            this->sys->setRestartedFromBackup(restartedFromBackup);
        }
    }

    BackupStep::BackupStep(std::shared_ptr<ParticleGroup>              pg,
                           std::shared_ptr<IntegratorManager>  integrator,
                           std::shared_ptr<ForceField>                 ff,
                           Parameters par,
                           std::string name):SimulationStepBase(pg,integrator,ff,par,name){
        input = sys->getInput();
    }

    BackupStep::~BackupStep(){
        //Write backup at the end of the simulation
        if(isSimulationStateCorrect()){
            this->writeBackup(this->sys->getBackupFilePath()+"End");
        } else {
            System::log<System::WARNING>("[BackupStep] (%s) Not writing end backup because simulation state is not correct",name.c_str());
        }
    }

    void BackupStep::init(cudaStream_t st){}

    void BackupStep::applyStep(ullint step, cudaStream_t st) {
        this->writeBackup(this->sys->getBackupFilePath());
    }



}}}

