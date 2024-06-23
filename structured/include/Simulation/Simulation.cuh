#ifndef __SIMULATION__
#define __SIMULATION__

#include<unistd.h>
#include<sys/wait.h>

namespace uammd{
namespace structured{

class Simulation{

    private:

        std::shared_ptr<ExtendedSystem>      sys;

        std::shared_ptr<GlobalData>           gd;
        std::shared_ptr<ExtendedParticleData> pd;

        /////////////////////////////////////////

        std::shared_ptr<Topology>           topology;
        std::shared_ptr<IntegratorManager>  integrators;

        std::shared_ptr<ForceField>         ff;

        std::shared_ptr<SimulationStepManager> simulationSteps;

    public:

        Simulation(std::shared_ptr<ExtendedSystem> sys):sys(sys){

            gd  = std::make_shared<GlobalData>(sys);
            pd  = std::make_shared<ExtendedParticleData>(sys);

            /////////////////////////////////////////

            if(sys->getRestartedFromBackup()){
                System::log<System::WARNING>("[Simulation] Simulation loaded from backup file.");

                System::log<System::WARNING>("[Simulation] Restarting from step: %llu", sys->getBackupRestartStep());

                //We change the seed to avoid the same random numbers in the next run
                System::log<System::WARNING>("[Simulation] Changing seed to avoid the same random numbers in the next run...");
                System::log<System::DEBUG>("[Simulation] Old seed: %llu", sys->getSeed());
                this->sys->setSeed(this->sys->getSeed()+ullint(time(NULL)));
                System::log<System::DEBUG>("[Simulation] New seed: %llu", sys->getSeed());
            }

            /////////////////////////////////////////

            //Load topology
            topology = std::make_shared<Topology>(sys, gd, pd);

            //Load force field
            ff = std::make_shared<ForceField>(topology);

            //Load integrators
            //Note that the integrators are loaded after the topology and the force field
            //this is because topology can set some particle properties that are needed
            //by integrators initialization.
            //For example, the particle mass and the particle radius.
            integrators = std::make_shared<IntegratorManager>(sys, gd, pd);

            /////////////////////////////////////////

            //Load simulations steps

            simulationSteps = std::make_shared<SimulationStepManager>(integrators,ff);

            //If saveBackup, add backupStep
            if(this->sys->getSaveBackup()){
                System::log<System::MESSAGE>("[Simulation] Adding backup step...");

                typename Backup::BackupStep::Parameters backupStepParam;

                backupStepParam.startStep    = this->sys->getBackupStartStep();
                backupStepParam.endStep      = this->sys->getBackupEndStep();
                backupStepParam.intervalStep = this->sys->getBackupIntervalStep();

                std::shared_ptr<Backup::BackupStep> backupStep =
                std::make_shared<Backup::BackupStep>(topology->getParticleGroup(),
                                                     integrators,
                                                     ff,
                                                     backupStepParam,
                                                     "backup");

                simulationSteps->addSimulationStep(backupStep,backupStep->getName());

            }
        }

        //Getters

        std::shared_ptr<ExtendedSystem>       getSystem(){ return this->sys;}
        std::shared_ptr<GlobalData>           getGlobalData(){ return this->gd;}
        std::shared_ptr<ExtendedParticleData> getParticleData(){return this->pd;}

        std::shared_ptr<Topology>   getTopology(){return this->topology;}
        std::shared_ptr<ForceField> getForceField(){return this->ff;}

        std::shared_ptr<IntegratorManager>   getIntegrators(){return this->integrators;}

        std::shared_ptr<SimulationStepManager> getSimulationSteps(){return this->simulationSteps;}

        //Run simulation

        int run(){

            Timer tim;
            tim.tic();

            std::map<std::string,
                     std::shared_ptr<SimulationStep::SimulationStepBase>> simSteps = simulationSteps->getSimulationSteps();

            try {

                if(this->getSystem()->getInput()->checkAllPathsUsed()){
                    System::log<System::DEBUG>("[Simulation] All input paths were used.");
                }else{
                    System::log<System::CRITICAL>("[Simulation] Some input paths were not used.");
                }

                System::log<System::MESSAGE>("[Simulation] Running simulation...");

                ullint previousIntegratorSteps = 0;
                for(auto& integratorInfo: integrators->getSortedIntegratorSteps()){


                    std::string name  = integratorInfo.name;
                    ullint      steps = integratorInfo.steps;

                    if(steps != std::numeric_limits<ullint>::max()){
                        System::log<System::MESSAGE>("[Simulation] Running integrator (%s), %llu steps...",
                                                     name.c_str(), steps);
                    } else {
                        System::log<System::MESSAGE>("[Simulation] Running integrator (%s) for an infinite number of steps...",
                                                     name.c_str(), steps);
                    }

                    std::shared_ptr<uammd::Integrator> currentIntegrator = integrators->getIntegrator(name);

                    //Load force field into integrators
                    currentIntegrator->addInteractor(ff);

                    ullint startStep = (gd->getFundamental()->getCurrentStep()-previousIntegratorSteps);
                    if(startStep > 0){
                        System::log<System::WARNING>("[Simulation] Integrator (%s) has already run %llu steps.",
                                                     name.c_str(), startStep);
                    }
                    if(startStep >= steps){
                        System::log<System::WARNING>("[Simulation] Integrator (%s) has already run more steps than total steps (%llu), skipping...",
                                                     name.c_str(), steps);
                    }
                    System::log<System::DEBUG>("[Simulation] Initializing simulation steps...");
                    for(auto sStep : simSteps){
                        sStep.second->tryInit();
                    }

                    for(ullint i = startStep; i < steps; i++){
                        #ifdef PYBIND11_MODULE
                        if (PyErr_CheckSignals() != 0) {
                            System::log<System::WARNING>("[Simulation] Python signal received. Exiting ...");
                            return 2;
                        }
                        #endif
                        for(auto sStep : simSteps){
                            sStep.second->tryApplyStep();
                        }
                        if(this->getSystem()->getState() == ExtendedSystem::SIMULATION_STATE::STOPPED){
                            System::log<System::MESSAGE>("[Simulation] Integrator (%s) stopped.",
                                                         name.c_str());
                            break;
                        }
                        currentIntegrator->forwardTime();
                    }
                    cudaDeviceSynchronize();

                    if(this->getSystem()->getState() == ExtendedSystem::SIMULATION_STATE::STOPPED){
                        System::log<System::MESSAGE>("[Simulation] Simulation stopped.");
                        break;
                    }

                    previousIntegratorSteps += steps;
                }

            } catch (std::exception& e) {

                System::log<System::WARNING>("[Simulation] Error running simulation: %s", e.what());
                System::log<System::WARNING>("[Simulation] Error at step %llu", gd->getFundamental()->getCurrentStep());

                if(sys->getSaveBackup()){

                    if(sys->getLastBackupStep() < sys->getBackupIntervalStep()){
                        System::log<System::WARNING>("[Simulation] Simulation failed before the first backup step");
                        return 2;
                    }

                    if(sys->getRestartedFromBackup()){
                        if(sys->getLastBackupStep() == sys->getBackupRestartStep()){
                            System::log<System::WARNING>("[Simulation] Simulation started from backup file and failed before a new backup step");
                            return 2;
                        } else {
                            System::log<System::WARNING>("[Simulation] Simulation started from backup file and failed after a new backup step");
                            return 1;
                        }
                    } else {
                        System::log<System::WARNING>("[Simulation] Simulation started from scratch and failed after a new backup step");
                        return 1;
                    }
                } else {
                    System::log<System::WARNING>("[Simulation] No backup detected.");
                    return 2;
                }
            }

            auto totalTime = tim.toc();

            System::log<System::MESSAGE>("[Simulation] Mean FPS: %f",
                                         real(gd->getFundamental()->getCurrentStep())/totalTime);

            return 0; //Success
        }
};

template<class inTyp>
void startSelfStartingSimulation(const inTyp& in){

    std::string backupFilePath = "";
    {
        bool saveBackup = false;
        auto input = std::make_shared<typename ExtendedSystem::InputType>(in);

        //Check if simulation has backup
        for(std::string entryName : input->getEntriesList({"system"})){
            std::vector<std::string> entryPath = {"system"};
            entryPath.push_back(entryName);

            auto dataEntry = input->getDataEntry(entryPath);

            std::string entryType    = dataEntry.getType();
            std::string entrySubType = dataEntry.getSubType();

            if(entryType == "Simulation"){
                if(entrySubType == "Backup"){
                    saveBackup = true;
                    backupFilePath = dataEntry.template getParameter<std::string>("backupFilePath");
                    if (backupFilePath.find(".") == std::string::npos){
                        backupFilePath += input->getFileExtension();
                    }
                }
            }
        }

        if(saveBackup){
            System::log<System::MESSAGE>("[Simulation] Backup enabled, backup file: %s", backupFilePath.c_str());
        } else {
            System::log<System::WARNING>("[Simulation] Simulation has no backup, will not restart if it fails.");
        }
    }

    pid_t pid;

    pid = fork();

    if(pid >= 0){
        //fork succeeded

        if(pid == 0){
            // child process
            int exitCode = 2; //Default exit code is 2 (simulation failed)
            {
                std::shared_ptr<ExtendedSystem> sys = std::make_shared<ExtendedSystem>(in);
                try {

                    {
                        std::shared_ptr<Simulation>  sim = std::make_shared<Simulation>(sys);
                        exitCode = sim->run();
                        uammd::System::log<uammd::System::MESSAGE>("[Simulation] Simulation returned exit code %d", exitCode);
                        uammd::System::log<uammd::System::MESSAGE>("[Simulation] Finalizing simulation...");
                    }

                    try {
                        sys->finish();
                    } catch (std::exception& e) {
                        uammd::System::log<uammd::System::WARNING>("[Simulation] Error finalizing simulation: %s", e.what());
                        std::exit(exitCode);
                    }
                } catch (std::exception& e) {
                    uammd::System::log<uammd::System::WARNING>("[Simulation] Error creating/deleting simulation: %s", e.what());
                    std::exit(exitCode);
                }

                if(exitCode != 0){
                    uammd::System::log<uammd::System::MESSAGE>("[Simulation] Simulation failed, returning exit code %d", exitCode);
                    std::exit(exitCode);
                }

            } //Need to destroy the simulation before exiting
            uammd::System::log<uammd::System::MESSAGE>("[Simulation] Simulation finished successfully");
            std::exit(exitCode);

        } else {
            //parent process
            int status;
            waitpid(pid,&status,0);
            if(WIFEXITED(status)){
                int childExitStatus = WEXITSTATUS(status);
                if        (childExitStatus == 0){
                    uammd::System::log<uammd::System::MESSAGE>("Simulation finished successfully");
                } else if (childExitStatus == 1){ //This code is returned when the simulation fails but a backup is available
                    uammd::System::log<uammd::System::WARNING>("Simulation finished with error code %d", childExitStatus);
                    uammd::System::log<uammd::System::WARNING>("Trying to restart simulation");

                    //Check the backup file exists
                    if(std::ifstream(backupFilePath)){
                        uammd::System::log<uammd::System::MESSAGE>("Backup file found, restarting simulation");
                        startSelfStartingSimulation(backupFilePath);
                    } else {
                        uammd::System::log<uammd::System::WARNING>("Backup file (%s) not found, aborting simulation", backupFilePath.c_str());
                    }
                } else {
                    uammd::System::log<uammd::System::WARNING>("Simulation finished with error code %d", childExitStatus);
                    uammd::System::log<uammd::System::CRITICAL>("Not trying to restart simulation. Aborting");
                }
            } else {
                uammd::System::log<uammd::System::CRITICAL>("Child process did not exit normally!");
            }
        }
    } else {
        //fork failed
        uammd::System::log<uammd::System::CRITICAL>("Fork failed!");
    }

}

void startSelfStartingSimulationFromFile(std::string inputFilePath){
    startSelfStartingSimulation(inputFilePath);
}

void startSelfStartingSimulationFromInput(const typename ExtendedSystem::InputType::DataType& in){
    startSelfStartingSimulation(in);
}

}}

#endif
