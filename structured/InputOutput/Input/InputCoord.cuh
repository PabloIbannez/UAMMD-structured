#ifndef __INPUT_COORD__
#define __INPUT_COORD__

namespace uammd{
namespace structured{
namespace InputOutput{

namespace InputCoord_ns{

template<class coordFormat>
coordFormat processCoordLine(std::string& line){

    std::stringstream ss;
    ss.str(line);

    std::string pName;
    
    coordFormat coordBuffer;

    coordBuffer.vel = make_real3(0);

    ss >> coordBuffer.id    >>
          coordBuffer.pos.x >>
          coordBuffer.pos.y >>
          coordBuffer.pos.z >>
          coordBuffer.vel.x >>
          coordBuffer.vel.y >>
          coordBuffer.vel.z >>
          coordBuffer.dir.x >>
          coordBuffer.dir.y >>
          coordBuffer.dir.z >>
          coordBuffer.dir.w ;

    return coordBuffer;
}

}

template<class coordFormat>
std::vector<coordFormat> loadCoordFromBlock(std::shared_ptr<System> sys,std::string blocksFilePath,std::string blockLabel){
    
    std::vector<coordFormat> coord;

    InputBlocksFile blocksFile =  InputBlocksFile(sys,blocksFilePath);

    auto coordBlock = blocksFile.getFileBlockIterator(blockLabel);

    std::string line;
    std::stringstream parser;

    while (coordBlock.next(line)){
        
        if (!checkCommented(line)){

            coordFormat coordBuffer = InputCoord_ns::processCoordLine<coordFormat>(line);

            coord.push_back(coordBuffer);
            
            sys->log<System::DEBUG>("[Load coord] Added particle: pId:%i, pos:%f,%f,%f, vel:%f,%f,%f, dir:%f,%f,%f,%f",coordBuffer.id,         
                                                                                                                       coordBuffer.pos.x,       
                                                                                                                       coordBuffer.pos.y,       
                                                                                                                       coordBuffer.pos.z, 
                                                                                                                       coordBuffer.vel.x,       
                                                                                                                       coordBuffer.vel.y,       
                                                                                                                       coordBuffer.vel.z, 
                                                                                                                       coordBuffer.dir.x,       
                                                                                                                       coordBuffer.dir.y,       
                                                                                                                       coordBuffer.dir.z, 
                                                                                                                       coordBuffer.dir.w); 
        }

    }

    return coord;
}

template<class coordFormat>
std::vector<coordFormat> loadCoordFromFile(std::shared_ptr<System> sys,std::string coordFilePath){
    
    std::vector<coordFormat> coord;

    std::fstream coordFile(coordFilePath);

    if (coordFile.fail()) {
        sys->log<System::CRITICAL>("[Load coord] System file %s cannot be opened.", coordFilePath.c_str());
    } else {
        sys->log<System::MESSAGE>("[Load coord] System file %s opened.", coordFilePath.c_str());
    }

    std::string line;

    while (std::getline(coordFile, line)){
        
        if (!checkCommented(line)){

            coordFormat coordBuffer = InputCoord_ns::processCoordLine<coordFormat>(line);

            coord.push_back(coordBuffer);
            
            sys->log<System::DEBUG>("[Load coord] Added particle: pId:%i, pos:%f,%f,%f, vel:%f,%f,%f, dir:%f,%f,%f,%f",coordBuffer.id,         
                                                                                                                       coordBuffer.pos.x,       
                                                                                                                       coordBuffer.pos.y,       
                                                                                                                       coordBuffer.pos.z, 
                                                                                                                       coordBuffer.vel.x,       
                                                                                                                       coordBuffer.vel.y,       
                                                                                                                       coordBuffer.vel.z, 
                                                                                                                       coordBuffer.dir.x,       
                                                                                                                       coordBuffer.dir.y,       
                                                                                                                       coordBuffer.dir.z, 
                                                                                                                       coordBuffer.dir.w); 
        }
    }
    
    sys->log<System::MESSAGE>("[Load coord] Coord file %s loaded into buffer (number particles = %i).", coordFilePath.c_str(),coord.size());

    return coord;
}

}}}

#endif
