#pragma once

#include <fstream>

#include"System/System.h"

#include"InputOutput/Input/InputFormats/InputJSON.cuh"

namespace uammd{
namespace structured{

    template<class InputType_>
    class ExtendedSystem_ : public uammd::System {

        public:

            using InputType = InputType_;

            //Enum available states, running and stopped
            enum SIMULATION_STATE {RUNNING, STOPPED};

        private:

            //Simulation state
            SIMULATION_STATE state = RUNNING;

            bool cudaStreamCreated = false;
            cudaStream_t stream;

            std::shared_ptr<InputType> input;

            ///////////////////////////
            std::vector<std::string> path;

            std::vector<std::string> infoPath;
            std::vector<std::string> backupPath;

            std::string name;
            ullint seed;

            //Backup
            bool        saveBackup;
            std::string backupFilePath;

            bool    restartedFromBackup;
            ullint  backupRestartStep;

            ullint lastBackupStep;

            ullint backupStartStep;
            ullint backupEndStep;
            ullint backupIntervalStep;

            void loadSimulationInformation(std::string entryName);
            void loadSimulationBackup(std::string entryName);

            void init();

        public:

            ExtendedSystem_(std::string inputFilePath);
            ExtendedSystem_(int argc, char *argv[],std::string inputFilePath,std::vector<std::string> path);

            cudaStream_t getCudaStream(){
                if(!cudaStreamCreated){
                    System::log<System::MESSAGE>("[ExtendedSystem] (%s) Creating cuda stream",path.back().c_str());
                    CudaSafeCall(cudaStreamCreate(&stream));
                    cudaStreamCreated = true;
                }
                return stream;
            }

            //void synchronizeCudaStream(bool force=false){

            //    if(!cudaStreamCreated){
            //        System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to synchronize a non created stream",path.back().c_str());
            //    }

            //    if(force){
            //        cudaStreamSynchronize(stream);
            //    } else {
            //        cudaStreamSynchronize(stream);
            //    }
            //}

            SIMULATION_STATE getState(){
                return state;
            }

            void setState(SIMULATION_STATE state){
                this->state = state;
            }

            std::shared_ptr<InputType> getInput(){
                return input;
            }

            //Getters
            std::string getName(){ return name; }
            ullint getSeed(){ return seed; }

            //Setters
            void setSeed(ullint seed){ this->seed = seed; }

            //Backup getters
            bool getSaveBackup(){ return saveBackup; }

            std::string getBackupFilePath(){
                if(!saveBackup){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to get backup file path,"
                                                  " but the simulation is not handling a backup!",path.back().c_str());
                }
                return backupFilePath;
            }

            bool getRestartedFromBackup(){ return restartedFromBackup; }

            ullint getBackupRestartStep(){
                if(!restartedFromBackup){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to get backup restart step, "
                                                  "but the simulation is not a backup!",path.back().c_str());
                }
                return backupRestartStep;
            }

            ullint getLastBackupStep(){
                if(!saveBackup){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to get last backup step,"
                                                  " but the simulation is not handling a backup!",path.back().c_str());
                }
                return lastBackupStep;
            }

            ullint getBackupStartStep(){
                if(!saveBackup){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to get backup start step,"
                                                  " but the simulation is not handling a backup!",path.back().c_str());
                }
                return backupStartStep;
            }

            ullint getBackupEndStep(){
                if(!saveBackup){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to get backup end step,"
                                                  " but the simulation is not handling a backup!",path.back().c_str());
                }
                return backupEndStep;
            }

            ullint getBackupIntervalStep(){
                if(!saveBackup){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to get backup interval step,"
                                                  " but the simulation is not handling a backup!",path.back().c_str());
                }
                return backupIntervalStep;
            }

            //Backup setters

            void setRestartedFromBackup(bool restartedFromBackup){
                if(!saveBackup){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to set restartedFromBackup,"
                                                  " but the simulation is not handling a backup!",path.back().c_str());
                }
                this->restartedFromBackup = restartedFromBackup;
            }

            void setBackupRestartStep(ullint backupRestartStep){
                if(!restartedFromBackup){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to set backup restart step,"
                                                  " but the simulation is not a backup!",path.back().c_str());
                }
                this->backupRestartStep = backupRestartStep;
            }

            void setLastBackupStep(ullint lastBackupStep){
                if(!saveBackup){
                    System::log<System::CRITICAL>("[ExtendedSystem] (%s) Trying to set last backup step,"
                                                  " but the simulation is not handling a backup!",path.back().c_str());
                }
                this->lastBackupStep = lastBackupStep;
            }

            void updateInputSystem();

            ///////////////////////////////////////////

            void finish();
    };

    using ExtendedSystem = ExtendedSystem_<InputJSON::InputJSON>;

    using ParametersEntry = ExtendedSystem::InputType::parametersEntry;
    using TypeEntry       = ExtendedSystem::InputType::typeEntry;
    using DataEntry       = ExtendedSystem::InputType::dataEntry;

    std::shared_ptr<ExtendedSystem> getExtendedSystem(std::shared_ptr<uammd::System> sys);
}}
