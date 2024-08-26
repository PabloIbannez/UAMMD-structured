#include "Topology/Topology.cuh"

namespace uammd{
namespace structured{

    void Topology::loadStructure(){

        std::vector<std::string> structurePath = path;
        structurePath.push_back("structure");

        auto structureData = sys->
                             getInput()->
                             getDataEntry(structurePath);

        std::vector<int> id = structureData.getData<int>("id");

        int N = id.size();

        if (N != pd->getNumParticles()){
            System::log<System::CRITICAL>("[Topology] Number of particles in the structure file "
                                          "does not match the number of particles in the system.");
        }

        std::vector<std::string> type = structureData.getData<std::string>("type");

        ////////////////////////////////////////////////////////////////////////////

        std::vector<int> resBuffer;
        if(structureData.isDataAdded("resId")){
            resBuffer = structureData.getData<int>("resId");
        } else {
            resBuffer.resize(N);
            std::fill(resBuffer.begin(),resBuffer.end(),0);
        }

        std::vector<int> chainBuffer;
        if(structureData.isDataAdded("chainId")){
            chainBuffer = structureData.getData<int>("chainId");
        } else {
            chainBuffer.resize(N);
            std::fill(chainBuffer.begin(),chainBuffer.end(),0);
        }

        std::vector<int> mdlBuffer;
        if(structureData.isDataAdded("modelId")){
            mdlBuffer = structureData.getData<int>("modelId");
        } else {
            mdlBuffer.resize(N);
            std::fill(mdlBuffer.begin(),mdlBuffer.end(),0);
        }

        std::vector<int> batchBuffer;
        if(structureData.isDataAdded("batchId")){
            batchBuffer = structureData.getData<int>("batchId");
        } else {
            batchBuffer.resize(N);
            std::fill(batchBuffer.begin(),batchBuffer.end(),0);
        }

        ////////////////////////////////////////////////////////////////////////////

        try {

            auto typeParamHandler = gd->getTypes();

            const int * sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

            auto pos  = pd->getPos(access::location::cpu,     access::mode::write);
            auto res  = pd->getResId(access::location::cpu,   access::mode::write);
            auto chn  = pd->getChainId(access::location::cpu, access::mode::write);
            auto mdl  = pd->getModelId(access::location::cpu, access::mode::write);
            auto bat  = pd->getBatchId(access::location::cpu, access::mode::write);

            for(int i=0;i<N;i++){

                int id_ = id[i];

                pos[sortedIndex[id_]].w = int(typeParamHandler->getTypeId(type[i]));
                res[sortedIndex[id_]]   = resBuffer[i];
                chn[sortedIndex[id_]]   = chainBuffer[i];
                mdl[sortedIndex[id_]]   = mdlBuffer[i];
                bat[sortedIndex[id_]]   = batchBuffer[i];

                System::log<System::DEBUG2>("[Topology] Added structure for particle %i (index %i),"
                                            "type: %s, res:%i, chain:%i, mol:%i, batch:%i",
                                            id_,sortedIndex[id_],type[i].c_str(),
                                            resBuffer[i],chainBuffer[i],mdlBuffer[i],batchBuffer[i]);
            }

        } catch(...) {
            System::log<System::CRITICAL>("[Topology] Error reading structure.");
        }

        {
            System::log<System::MESSAGE>("[Topology] Loading types into particle data.");
            gd->getTypes()->loadTypesIntoParticleData(pd);
        }
    }

    //Force field
    void Topology::loadGroups(){
        groups = GroupUtils::loadGroupsListFromInputEntries(sys,gd,pd,this->forceFieldInfo);
    }

    void Topology::loadNeighbourLists(){
        VConListSet = VerletConditionalListSetUtils::loadVerletConditionalListSetsFromInputEntries(sys,gd,groups,this->forceFieldInfo);
    }

    void Topology::loadInteractors(){

        for(auto& entry : forceFieldInfo->getEntriesInfo()){

            if(Interactor::isInteractorAvailable(sys,entry.second.path)){
                std::shared_ptr<typename uammd::Interactor> inter = Interactor::loadInteractor(sys,gd,groups,VConListSet,entry.second.path);

                if(interactors.count(entry.second.name) == 0){
                    interactors[entry.second.name] = inter;
                    entry.second.used = true;
                } else {
                    System::log<System::CRITICAL>("[Topology] Error loading interactors,"
                                                  "interactor \"%s\" has already been added.",entry.second.name.c_str());
                }

            } else {

                DataEntry data = sys->getInput()->getDataEntry(entry.second.path);

                std::string interactorType    = data.getType();
                std::string interactorSubType = data.getSubType();

                System::log<System::CRITICAL>("[Topology] Error loading interactor: \"%s\". "
                                              "Interactor (type,subtype): (\"%s\",\"%s\") not recognized.",
                                              entry.second.name.c_str(),
                                              interactorType.c_str(),interactorSubType.c_str());
            }
        }

        //Print information about interactors
        for(auto i : interactors){
            System::log<System::MESSAGE>("[Topology] Added interactor \"%s\".",
                                         i.first.c_str());
        }
    }

    Topology::Topology(std::shared_ptr<ExtendedSystem>       sys,
                       std::shared_ptr<GlobalData>            gd,
                       std::shared_ptr<ExtendedParticleData>  pd,
                       std::vector<std::string> path):sys(sys),gd(gd),pd(pd),path(path){

        this->loadStructure();

        //Load components

        std::vector<std::string> pathForceField = path;
        pathForceField.push_back("forceField");
        forceFieldInfo = std::make_shared<InputEntryManager>(sys,pathForceField);

        this->loadGroups();
        this->loadNeighbourLists();
        this->loadInteractors();

        forceFieldInfo->checkEntriesUsed();
    }

    Topology::Topology(std::shared_ptr<ExtendedSystem>       sys,
                       std::shared_ptr<GlobalData>            gd,
                       std::shared_ptr<ExtendedParticleData>  pd):Topology(sys,gd,pd,{"topology"}){}

    //Add a new interactor to the system
    void Topology::addInteractor(std::shared_ptr<typename uammd::Interactor> interactor, std::string name){
        if(interactors.count(name) == 0){
            System::log<System::MESSAGE>("[Topology] Added interactor \"%s\" (manually).",
                                         name.c_str());
            interactors[name] = interactor;
        } else {
            System::log<System::CRITICAL>("[Topology] Error adding interactor manually,"
                                          "interactor \"%s\" has already been added.",name.c_str());
        }
    }

    //Get interactor
    std::map<std::string,std::shared_ptr<typename uammd::Interactor>>
    Topology::getInteractorsByClass(std::string intClass){
        //Warning: manually added interactors will not be included in this list
        System::log<System::WARNING>("[Topology] Getting interactors by class. "
                                     "Manually added interactors are not included in this list.");

        std::map<std::string,std::shared_ptr<typename uammd::Interactor>> interactorsByClass;
        for(std::string entryName : forceFieldInfo->getEntriesByClass(intClass)){
            if(interactors.count(entryName) == 0){
                System::log<System::WARNING>("[Topology] Interactor \"%s\" (class: \"%s\")"
                                             " is present if the force field but it was not loaded."
                                             " Not including it in the list.",
                                             entryName.c_str(),intClass.c_str());
            }
            interactorsByClass[entryName] = interactors[entryName];
        }
        if(interactorsByClass.size() == 0){
            System::log<System::WARNING>("[Topology] No interactors of class \"%s\" "
                                         "were found in the force field.",intClass.c_str());
        }
        return interactorsByClass;
    }

}}


