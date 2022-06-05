#ifndef __TOPOLOGY__
#define __TOPOLOGY__

namespace uammd{
namespace structured{

template<class Units_,
         class Types_>
class Topology: public InputOutput::InputBlocksFile{
    
    private:
        
        std::string labelDelimiter = ";";
        std::string paramDelimiter = ",";
        std::string assocDelimiter = "=";

    public:

        using Units            = Units_;
        using Types            = Types_;

        using TypeParameterHandler = TypeParameterHandler<Types>;

    private:

        struct entryInfo{
            
            std::string label;

            std::string alias;
            std::string type;
            std::string name;
            std::map<std::string,std::string> param;

        };
        
        std::map<std::string,entryInfo> entries;

        std::shared_ptr<TypeParameterHandler> typeParamH;
        
        std::shared_ptr<TypeParameterHandler> readTypes(std::string typesLabel){
            
            std::shared_ptr<TypeParameterHandler> types = std::make_shared<TypeParameterHandler>();
            
            fileBlockIterator typesBlock = this->getFileBlockIterator(typesLabel);

            typename TypeParameterHandler::InputTypeParameters typeParamBuffer;
            int typeId=0;

            //Auxiliar vector to check that a type is not added twice
            std::vector<std::string> typesNameList;
            
            std::string line;
            std::stringstream parser;

            while (typesBlock.next(line)){
                
                if (!InputOutput::checkCommented(line)){

                    typeParamBuffer = types->readTypeParameters(line);

                    if(std::count(typesNameList.begin(), typesNameList.end(), typeParamBuffer.name) != 0){
                        sys->log<System::CRITICAL>("[Topology] Error loading particles types. Particle type %s has been added previosly", typeParamBuffer.name.c_str());
                    } else {
                        typesNameList.push_back(typeParamBuffer.name);
                        types->add(typeId,typeParamBuffer);
                        sys->log<System::MESSAGE>("[Topology] Particle type %s (type id: %i) added.", typeParamBuffer.name.c_str(),typeId);
                    }
                    //Increase typeId counter
                    typeId ++;
                }
            }

            return types;

        }
        
        void loadTypeParameterHandler(std::string typesLabel){
            typeParamH = this->readTypes(typesLabel);
        }

    public:
        
        Topology(std::shared_ptr<System> sys,
                 std::string topologyFilePath):InputBlocksFile(sys,topologyFilePath){
            sys->log<System::MESSAGE>("[Topology] Parameter inputTopologyPath added: %s",topologyFilePath.c_str());

            ///

            std::vector<std::string> entriesList = this->getBlocksList();

            for(std::string label : entriesList){
                std::vector<std::string> labelBuffer = Miscellany::split(label,labelDelimiter);

                entries[label].label = label;

                //
                
                entries[label].alias = labelBuffer[0];
                
                if(labelBuffer.size()>1){
                    entries[label].type = labelBuffer[1];
                }

                if(labelBuffer.size()>2){
                    entries[label].name = labelBuffer[2];
                }

                if(labelBuffer.size()>3){
                    std::vector<std::string> paramBuffer = Miscellany::split(labelBuffer[3],paramDelimiter);
                    for(std::string p : paramBuffer){
                        std::vector<std::string> pBuffer = Miscellany::split(p,assocDelimiter);
                        entries[label].param[pBuffer[0]] = pBuffer[1];
                    }
                }
            }

            /*
            for(auto l : entries){
                std::cout << "key "     << l.first << 
                             " label: " << l.second.label << 
                             " alias: " << l.second.alias << 
                             " type: "  << l.second.type  << 
                             " name: "  << l.second.name  << " params: ";

                for(auto p : l.second.param){
                    std::cout << p.first << "," << p.second << "      ";
                }
                std::cout << std::endl;
            }
            
            std::cin.get();
            */

            ////
            
            std::set<std::string> aliasSet;
            for(auto& entry : entries){
                std::string alias=entry.second.alias;
                if(aliasSet.count(alias)){
                    sys->log<System::CRITICAL>("[Topology] "
                                               "Alias %s has been added twice.",alias.c_str());
                } else {
                    aliasSet.insert(alias);
                }
            }

            this->loadTypeParameterHandler("TYPES");
        }
        
        Topology(std::shared_ptr<System> sys,
                 uammd::InputFile& in):Topology(sys,in.getOption("inputTopologyPath",InputFile::Required).str()){}
        
        std::shared_ptr<TypeParameterHandler> getTypes(){return typeParamH;}

        bool isEntryPresent(std::string alias,std::string type,std::string name){
            for(auto& entry : entries){
                if(((entry.second.alias == alias) or (alias.empty())) and 
                   ((entry.second.type  == type ) or ( type.empty())) and 
                   ((entry.second.name  == name ) or ( name.empty())) ){
                    return true;
                }
            }
            return false;
        }
        
        bool isEntryPresent(std::string alias){
            return isEntryPresent(alias,std::string(),std::string());
        }
        
        bool isEntryPresent(std::string type,std::string name){
            return isEntryPresent(std::string(),type,name);
        }
        
        std::vector<entryInfo> getEntryInfo(std::string alias,std::string type,std::string name){
            std::vector<entryInfo> entriesInfo;
            for(auto& entry : entries){
                if(((entry.second.alias == alias) or (alias.empty())) and 
                   ((entry.second.type  == type ) or ( type.empty())) and 
                   ((entry.second.name  == name ) or ( name.empty())) ){
                    entriesInfo.push_back(entry.second);
                }
            }
            if(entries.size()>0){
                return entriesInfo;
            }
            sys->log<System::CRITICAL>("[Topology] "
                                       "There is not entry with alias:%s ,type:%s and name:%s.",alias.c_str(),type.c_str(),name.c_str());
            return std::vector<entryInfo>();
        }
        
        std::vector<entryInfo> getEntryInfo(std::string alias){
            return getEntryInfo(alias,
                                std::string(),
                                std::string());
        }
        
        std::vector<entryInfo> getEntryInfo(std::string type,std::string name){
            return getEntryInfo(std::string(),
                                type,
                                name);
        }

        template<class InteractionParameters>
        std::shared_ptr<typename structured::PairParameterHandler<InteractionParameters>> readPairs(std::string pairLabel){
            
            using PairParameterHandler  = typename structured::PairParameterHandler<InteractionParameters>;
    
            std::shared_ptr<PairParameterHandler> pairParameterHandler = std::make_shared<PairParameterHandler>();

            fileBlockIterator pairBlock = this->getFileBlockIterator(pairLabel);
                    
            //Auxiliar map which translate from name(string) to typeId(int)(internal representation)
            std::map<std::string,int> nameIdMap;

            for(int typeId : typeParamH->getTypeIdList()){
                nameIdMap[typeParamH->getTypeParameters(typeId).name]=typeId;
            }

            //Determinate particle number and check if all particles types are correct
            std::string line;
            std::stringstream parser;
    
            typename PairParameterHandler::InputPairParameters pairParameterBuffer;
            std::string type_i,type_j;

            while (pairBlock.next(line)){
                
                if (!InputOutput::checkCommented(line)){

                    parser.clear();
                    parser.str(line);

                    parser >> type_i >> type_j;
                    //parser.str().substr(parser.tellg()) remainder
                    pairParameterBuffer = pairParameterHandler->readPairParameters(parser.str().substr(parser.tellg()));

                    if(nameIdMap.count(type_i) == 0){
                        sys->log<System::CRITICAL>("[Topology] "
                                "Particle type %s has not been added.", type_i.c_str());
                    } 

                    if(nameIdMap.count(type_j) == 0){
                        sys->log<System::CRITICAL>("[Topology] "
                                "Particle type %s has not been added.", type_j.c_str());
                    } 

                    pairParameterHandler->add(nameIdMap[type_i],
                                              nameIdMap[type_j],
                                              pairParameterBuffer);

                    sys->log<System::DEBUG>("[Topology] Particle pair interaction"
                                            " parameters added for pair %s %s", 
                                            type_i.c_str(),type_j.c_str());
                }
            }

            return pairParameterHandler;
        }
        
        void loadTypes(std::shared_ptr<ParticleData> pd){
            
            auto pos = pd->getPos(access::location::cpu,access::mode::write);

            fori(0,pd->getNumParticles()){
                Types::load(pd,i,typeParamH->getTypeParameters(int(pos[i].w)));
            }
        }

        void loadStructureData(std::shared_ptr<ParticleData> pd){
            loadStructureData(pd,"STRUCTURE");
        }

        void loadStructureData(std::shared_ptr<ParticleData> pd,std::string structureLabel){
            
            fileBlockIterator structureBlock = this->getFileBlockIterator(structureLabel);
                    
            const int * sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

            auto pos    = pd->getPos(access::location::cpu,     access::mode::write);
            auto resId  = pd->getResId(access::location::cpu,   access::mode::write);
            auto chnId  = pd->getChainId(access::location::cpu, access::mode::write);
            auto molId  = pd->getModelId(access::location::cpu, access::mode::write);
            auto simId  = pd->getSimulationId(access::location::cpu, access::mode::write);

            //Auxiliar map which translate from name(string) to typeId(int)(internal representation)
            std::map<std::string,int> nameIdMap;

            for(int typeId : typeParamH->getTypeIdList()){
                nameIdMap[typeParamH->getTypeParameters(typeId).name]=typeId;
            }

            //Determinate particle number and check if all particles types are correct
            std::string line;
            std::stringstream parser;

            bool first=true;
            bool loadSimId=false;
            while (structureBlock.next(line)){
                
                if(first){
                    //Count number of columns
                    parser.clear();
                    parser.str(line);
                    
                    uint ncol=0;
                    std::string colBuffer;
                    while(parser >> colBuffer){ncol++;}
                    
                    if(ncol == 6){
                        loadSimId=true;
                    }
                    first=false;
                }

                parser.clear();
                parser.str(line);

                int id;
                std::string name;  
                int resIdBuffer;
                int chnIdBuffer;
                int molIdBuffer;
                int simIdBuffer;
                
                parser >> id           >>
                          name         >>
                          resIdBuffer  >>
                          chnIdBuffer  >>
                          molIdBuffer  ;

                if(!loadSimId){
                    simIdBuffer=0;
                } else {
                    parser>>simIdBuffer;
                }

                if(nameIdMap.count(name) == 0){
                    sys->log<System::CRITICAL>("Error while loading structure, particle type %s has not been added.", name.c_str());
                } 

                if(parser.fail()){
                    sys->log<System::CRITICAL>("Error processing line \"%s\", while loading structure %s.", line.c_str(),structureLabel.c_str());
                }

                pos[sortedIndex[id]].w = int(nameIdMap[name]);
                resId[sortedIndex[id]] = resIdBuffer;
                chnId[sortedIndex[id]] = chnIdBuffer;
                molId[sortedIndex[id]] = molIdBuffer;
                simId[sortedIndex[id]] = simIdBuffer;

                sys->log<System::DEBUG2>("Added structure for particle %i (%i),type: %s (%i), res:%i, chain:%i, mol:%i, sim:%i",
                                          id,sortedIndex[id],name.c_str(),int(nameIdMap[name]),resIdBuffer,
                                                                                               chnIdBuffer,
                                                                                               molIdBuffer,
                                                                                               simIdBuffer);
            }
        }
        
        template<class propertyIterator>
        void loadProperty(std::shared_ptr<ParticleData> pd,
                          std::string propertyName,
                          propertyIterator propertyItr){

            fileBlockIterator propertyBlock = this->getFileBlockIterator(propertyName);

            const int * sortedIndex = pd->getIdOrderedIndices(access::location::cpu);

            std::string line;
            std::stringstream parser;

            while(propertyBlock.next(line)){

                parser.clear();
                parser.str(line);

                int  id;
                typename std::remove_pointer<typename propertyIterator::Iterator>::type propertyBuffer;

                parser >> id            >>
                          propertyBuffer;

                if(parser.fail()){
                    sys->log<System::CRITICAL>("Error processing line \"%s\", while loading property %s.", line.c_str(),propertyName.c_str());
                }

                propertyItr[sortedIndex[id]]=propertyBuffer;

                sys->log<System::DEBUG2>("Added property %s for particle %i (%i), value: %f",
                                          propertyName.c_str(),id,sortedIndex[id],propertyBuffer);
            }
        }

        void loadNeighbourList(std::string listName, std::map<int,std::vector<int>>& neighbourList){
            
            fileBlockIterator propertyBlock = this->getFileBlockIterator(listName);

            std::string line;
            std::stringstream parser;

            while(propertyBlock.next(line)){
                
                parser.clear();
                parser.str(line);
                
                int id;
                parser >> id;
                
                if(neighbourList.count(id)>0){
                    sys->log<System::CRITICAL>("Neighbours for particle %i have been " 
                                                "added previously, line : %s",id, line.c_str());
                } else {
                    int neig;
                    while(parser >> neig){
                        neighbourList[id].push_back(neig);
                    }
                }
            }
        }
};
    
}}

#endif
