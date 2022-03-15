#ifndef __BLOCKS_FILE__
#define __BLOCKS_FILE__

#include <regex>

namespace uammd{
namespace structured{
namespace InputOutput{

bool checkCommented(std::string& line){

    if (line.empty() or 
        line.find_first_not_of (' ') == line.npos or
        line[line.find_first_not_of(' ')] == '#'){

        return true;

    } else {

        return false;
    }

}

class InputBlocksFile{

    protected:

        std::shared_ptr<System> sys;
        
        std::string blocksFilePath;
        std::shared_ptr<std::ifstream> blocksFile;

        struct labeledBlock{
            std::ifstream::pos_type blockBegin;
            std::ifstream::pos_type blockEnd;
        };
        
        std::map<std::string,labeledBlock> labeledBlocks;
        

    public:
        
        class fileBlockIterator{

            private: 

                std::shared_ptr<std::ifstream> file;

                std::ifstream::pos_type begin;
                std::ifstream::pos_type end;

            public:

                fileBlockIterator(std::shared_ptr<std::ifstream> file,
                                  std::ifstream::pos_type begin,
                                  std::ifstream::pos_type end):file(file),begin(begin),end(end){
                    this->reset();
                }

                void reset(){
                    file->clear();
                    file->seekg(begin);
                }

                bool next(std::string& line){
                    if(std::getline(*file,line) and (file->tellg() < end or end == -1)){
                        if (line.empty() or 
                            line.find_first_not_of (' ') == line.npos or
                            line[line.find_first_not_of(' ')] == '#'){
                            
                            return this->next(line);
                        }
                        return true;
                    }  
                    return false;
                }
        };


        InputBlocksFile(std::shared_ptr<System>      sys,
                        std::string blocksFilePath):sys(sys),
                                                    blocksFilePath(blocksFilePath){

            blocksFile = std::make_shared<std::ifstream>(blocksFilePath);
                    
            if (blocksFile->fail()) {
                    sys->log<System::CRITICAL>("[Blocks file] Blocks file %s cannot be opened.", blocksFilePath.c_str());
            } else {
                    sys->log<System::MESSAGE>("[Blocks file] Blocks file %s opened.", blocksFilePath.c_str());
            }

            std::regex rgx("\\[\\s*(\\w*)\\s*\\]");

            std::string previousLabel;
            for(std::string line ; std::getline(*blocksFile,line);){
                
                std::smatch match;

                if(std::regex_search(line,match,rgx)){

                    if(labeledBlocks.count(match[1])>0){
                        sys->log<System::CRITICAL>("[Blocks file] Block %s in blocks file %s is repeated.",
                                                    match[1].str().c_str(),blocksFilePath.c_str());
                    } else {
                        labeledBlocks[match[1]].blockBegin=blocksFile->tellg();

                        if(!previousLabel.empty()){
                            labeledBlocks[previousLabel].blockEnd=labeledBlocks[match[1]].blockBegin;
                        }
                        
                        previousLabel = match[1];
                    }
                }
            
            }
            
            labeledBlocks[previousLabel].blockEnd=blocksFile->tellg();
        }


        ~InputBlocksFile(){
            blocksFile->close();
        }

        bool isBlockPresent(std::string blockName){
            if(labeledBlocks.count(blockName)==0){
                return false;
            }
            return true;
        }
        
        fileBlockIterator getFileBlockIterator(std::string blockName){
            
            if(!isBlockPresent(blockName)){
                sys->log<System::CRITICAL>("[Blocks file] The block named %s is not present in the blocks file %s.", 
                                            blockName.c_str(),blocksFilePath.c_str());
            }

            return fileBlockIterator(blocksFile,labeledBlocks[blockName].blockBegin,labeledBlocks[blockName].blockEnd);
        }

        void printBlock(std::ostream& out,std::string blockLabel){
            
            fileBlockIterator block = this->getFileBlockIterator(blockLabel);

            std::string line;
            while(block.next(line)){
                out << line << std::endl;
            }
        }
};

}}}

#endif
