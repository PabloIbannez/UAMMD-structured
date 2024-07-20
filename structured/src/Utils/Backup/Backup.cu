#include "Utils/Backup/Backup.cuh"

namespace uammd{
namespace structured{
namespace Backup{

    bool openFile(std::shared_ptr<ExtendedSystem> sys,
                  std::string outputFilePath,
                  std::ofstream& outputFile,
                  bool binary){

        bool isFileEmpty = false;

        if(sys->getRestartedFromBackup()){
            //Init from backup
            //Check if the outputFilePath exists.
            bool fileExists = false;
            {
                std::fstream file(outputFilePath);
                fileExists = file.good();
            }

            if(fileExists){
                //If file exists, open in append mode
                System::log<System::DEBUG>("[Backup] Backup detected and previos file found. Opening file (%s) in append mode", outputFilePath.c_str());
                if(binary){
                    outputFile.open(outputFilePath, std::ios::binary | std::ios::app);
                }else{
                    outputFile.open(outputFilePath, std::ios::app);
                }
            } else {
                System::log<System::DEBUG>("[Backup] Backup detected but not previous file found. Opening file (%s) in write mode", outputFilePath.c_str());
                if(binary){
                    outputFile.open(outputFilePath, std::ios::binary);
                }else{
                    outputFile.open(outputFilePath);
                }
                isFileEmpty = true;
            }
        } else {
            System::log<System::DEBUG>("[Backup] No backup detected. Opening file (%s) in write mode", outputFilePath.c_str());
            if(binary){
                outputFile = std::ofstream(outputFilePath, std::ios::binary);
            } else {
                outputFile = std::ofstream(outputFilePath);
            }
            isFileEmpty = true;
        }

        return isFileEmpty;
    }
}}}

