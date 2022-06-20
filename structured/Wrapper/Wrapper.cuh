#ifndef __WRAPPER__
#define __WRAPPER__

namespace uammd{
namespace structured{
namespace Wrapper{

std::shared_ptr<ParticleData> setUpParticleData(std::shared_ptr<System> sys,
                                                InputFile& in){
    
    std::string inputCoordPath;

    in.getOption("inputCoordPath"   ,InputFile::Required) >> inputCoordPath;

    struct coordFormat{
        int   id;
        real3 pos;
        real3 vel;
        real4 dir;
    };
    
    std::vector<coordFormat> pdBuffer = InputOutput::loadCoordFromFile<coordFormat>(sys,inputCoordPath);
            
    std::shared_ptr<ParticleData> pd = std::make_shared<ParticleData>(pdBuffer.size(),sys);
            
    auto pId = pd->getId(access::location::cpu, access::mode::write);
    auto pos = pd->getPos(access::location::cpu, access::mode::write);
    auto vel = pd->getVel(access::location::cpu, access::mode::write);
    auto dir = pd->getDir(access::location::cpu, access::mode::write);

    fori(0,pdBuffer.size()){

        if(i==pdBuffer[i].id){
            pId[i]   = pdBuffer[i].id;
            pos[i].x = pdBuffer[i].pos.x;
            pos[i].y = pdBuffer[i].pos.y;
            pos[i].z = pdBuffer[i].pos.z;
            vel[i].x = pdBuffer[i].vel.x;
            vel[i].y = pdBuffer[i].vel.y;
            vel[i].z = pdBuffer[i].vel.z;
            dir[i].x = pdBuffer[i].dir.x;
            dir[i].y = pdBuffer[i].dir.y;
            dir[i].z = pdBuffer[i].dir.z;
            dir[i].w = pdBuffer[i].dir.w;

        } else {
            sys->log<System::CRITICAL>("[Wrapper] The internal id has to "
                    "match with the given by coord file. Inte: %i, File %i",i,pdBuffer[i].id);
        }
    }

    return pd;
}

std::shared_ptr<ParticleGroup> setUpParticleGroup(std::shared_ptr<ParticleData>  pd,
                                                  InputFile& in){
            
    std::shared_ptr<ParticleGroup> pg = std::make_shared<ParticleGroup>(pd,"All");

    return pg;
}


template<class ForceField>        
std::shared_ptr<ForceField> setUpForceField(std::shared_ptr<ParticleGroup> pg,
                                            InputFile& in){
    
    std::shared_ptr<ForceField> ff = std::make_shared<ForceField>(pg,in);
            
    auto top = ff->getTopology();

    top->loadStructureData(pg->getParticleData());
    top->loadTypes(pg->getParticleData());

    return ff;
}

}}}

#endif
