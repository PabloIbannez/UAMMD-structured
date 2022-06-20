#ifndef __GROUP_UTILS__
#define __GROUP_UTILS__

namespace uammd{
namespace structured{
namespace groupUtils{

std::set<int> getModelList(std::shared_ptr<ParticleData> pd){
    
    std::set<int> mdlList;
    {
        auto mdlId = pd->getModelId(uammd::access::location::cpu,uammd::access::mode::read);

        fori(0,pd->getNumParticles()){
            mdlList.emplace(mdlId[i]);
        }
    }

    return mdlList;
}

std::set<int> getSimList(std::shared_ptr<ParticleData> pd){
    
    std::set<int> simList;
    {
        auto simId = pd->getSimulationId(uammd::access::location::cpu,uammd::access::mode::read);

        fori(0,pd->getNumParticles()){
            simList.emplace(simId[i]);
        }
    }

    return simList;
}

std::vector<std::shared_ptr<uammd::ParticleGroup>> getSimGroups(std::shared_ptr<ParticleData> pd){

    std::vector<std::shared_ptr<uammd::ParticleGroup>> simGroups;

    std::set<int> simList = getSimList(pd);

    for(const int& s : simList){

        selectors::simulationId selector(s);

        auto pg = std::make_shared<uammd::ParticleGroup>(selector,
                                                         pd,
                                                         "simId_"+std::to_string(s));
        simGroups.push_back(pg);
        
    }

    return simGroups;
}

std::vector<std::shared_ptr<uammd::ParticleGroup>> getSimModelGroups(std::shared_ptr<ParticleData> pd,
                                                                     int mdlId){

    std::vector<std::shared_ptr<uammd::ParticleGroup>> simMdlGroups;

    std::set<int> simList = getSimList(pd);

    for(const int& s : simList){

        selectors::simulationIdModelId selector(s,mdlId);

        auto pg = std::make_shared<uammd::ParticleGroup>(selector,
                                                         pd,
                                                         "simId_"+std::to_string(s)+"_mdl_"+std::to_string(mdlId));
        simMdlGroups.push_back(pg);
        
    }

    return simMdlGroups;
}

bool checkSimGroupsEqual(std::shared_ptr<ParticleData> pd,
                         std::vector<std::shared_ptr<uammd::ParticleGroup>> simGroups){
    
    auto sys = pd->getSystem();

    int N = simGroups[0]->getNumberParticles();
    for(auto& sg : simGroups){
        int Nnew = sg->getNumberParticles();
        if(Nnew != N){
            sys->log<uammd::System::WARNING>("[CheckSimGroupsEqual] The number of particles for the groups \"%s\" and \"%s\" does not match",
                                              simGroups[0]->getName().c_str(),sg->getName().c_str());    
            return false;
        }    
    }
    
    bool equal=true;
    
    auto sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

    auto pos    = pd->getPos(access::location::cpu,     access::mode::read);
    auto resId  = pd->getResId(access::location::cpu,   access::mode::read);
    auto chnId  = pd->getChainId(access::location::cpu, access::mode::read);
    auto molId  = pd->getModelId(access::location::cpu, access::mode::read);

    auto refGroupIndex  = simGroups[0]->getIndexIterator(access::location::cpu);
    
    fori(0,N){

        int refIndex = refGroupIndex[i];
        
        int refType = int(pos[refIndex].w);
        int refRes  = resId[refIndex];
        int refCh   = chnId[refIndex];
        int refMol  = molId[refIndex];

        for(auto& sg : simGroups){
            
            auto groupIndex = sg->getIndexIterator(access::location::cpu);
        
            int index = groupIndex[i];
        
            int type = int(pos[index].w);
            int res  = resId[index];
            int ch   = chnId[index];
            int mol  = molId[index];

            if((type != refType) or 
               (res  != refRes ) or 
               (ch   != refCh  ) or
               (mol  != refMol )){
                    sys->log<uammd::System::WARNING>("[CheckSimGroupsEqual] The particle (%i,%i,%i,%i) (\"%s\") and (%i,%i,%i,%i) (\"%s\") do not match",
                                                      refType,refRes,refCh,refMol,simGroups[0]->getName().c_str(),
                                                      type,res,ch,mol,sg->getName().c_str());

                    equal = false;
            }
        }
    }

    return equal;
}

}}}

#endif
