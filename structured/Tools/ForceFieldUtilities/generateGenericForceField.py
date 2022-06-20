import os,sys

import re

import CppHeaderParser

FORCE_FIELD_PATH = os.getcwd()+"/../../ForceFields/Generic/"

POTENTIAL_PATH    = os.getcwd()+"/../../Potentials/"
POTENTIAL_BONDS   = ["Bond"+str(i+1) for i in range(4)]
POTENTIAL_UNBOUND = ["UnBound"]

TEMPLATES = ["Units_","power","CommonPotentials::LennardJones::Type2","CommonPotentials::LennardJones::Type3","6","12"]
TEMPLATES_NEEDED = {"Units_":"Base::Units"}

def getBondPotentialInfoFromFile(filename):

    fileManager = CppHeaderParser.CppHeader(filename)
    
    potNames = [[j['name'],j['defaultValue']] for j in [i for i in fileManager.variables if i['type']=='using']]
    
    pot = {}
    for n in potNames:
        n[1] = n[1].replace(" ","")
        cl = n[1]

        pot[n[0]] = {"template":None}

        while ((templateStart := [m.start() for m in re.finditer('<',cl)]) and 
                 (templateEnd := [m.start() for m in re.finditer('>',cl)])):
            cl = cl[templateStart[0]+1:templateEnd[-1]]

        cl=cl.replace(" ","")
        #print(cl,cl in TEMPLATES)
        if cl in TEMPLATES:
            for i,info in enumerate(n[1].split("<")):
                if (cl in info):
                    pot[n[0]]['template'] = cl
                    cl=(n[1].split("<")[i-1]).replace(" ","")
                    #print("template detected",cl)

    
        pot[n[0]]['class']=cl.replace(" ","")
    
    for p in pot.keys():
        pot[p]['parameters']={}
    
        cl = pot[p]['class']
    
        #print(fileManager.classes.keys())
        parameters = [nested for nested in fileManager.classes[cl]['nested_classes'] if nested['name']=="Parameters"][0]
        while(True):
    
            for param in parameters['properties']['public']:
                if param['type'] != 'Box':
                    paramType = param['type']
                    if(paramType=="std::string"):
                        paramType="str"
                    pot[p]['parameters'][param['name']] = paramType
    
            inhe = parameters['inherits']
            if(inhe):
                baseName = inhe[0]['class']
        
                #print(baseName)
                if ((templateStart := [m.start() for m in re.finditer('<',baseName)]) and 
                    (templateEnd   := [m.start() for m in re.finditer('>',baseName)])):
                    sub = baseName[templateStart[0]:templateEnd[-1]+1]
                    baseName = baseName.replace(sub,'')
                    #print(baseName)
                
                parameters = [nested for nested in [fileManager.classes[baseName]] if nested['name']=="Parameters"][0]
            else:
                break

    return pot

def getUnBoundPotentialInfoFromFile(filename):

    fileManager = CppHeaderParser.CppHeader(filename)
    
    potNames = [[j['name'],j['defaultValue']] for j in [i for i in fileManager.variables if i['type']=='using']]
    
    pot = {}
    #Get class name
    for n in potNames:
        n[1] = n[1].replace(" ","")
        cl = n[1]

        cl=cl.split("<")[0]
        pot[n[0]] = {'class':cl.replace(" ","")}

        pot[n[0]]["template"]=None
    
    for p in pot.keys():
        pot[p]['parameters']={}
        cl = pot[p]['class']
    
        #print(fileManager.classes.keys())
        parameters = [nested for nested in fileManager.classes[cl]['nested_classes'] if nested['name']=="Parameters"][0]
        while(True):
    
            for param in parameters['properties']['public']:
                if param['type'] != 'Box':
                    paramType = param['type']
                    if(paramType=="std::string"):
                        paramType="str"
                    pot[p]['parameters'][param['name']] = paramType
    
            inhe = parameters['inherits']
            if(inhe):
                baseName = inhe[0]['class']
        
                #print(baseName)
                if ((templateStart := [m.start() for m in re.finditer('<',baseName)]) and 
                    (templateEnd   := [m.start() for m in re.finditer('>',baseName)])):
                    sub = baseName[templateStart[0]:templateEnd[-1]+1]
                    baseName = baseName.replace(sub,'')
                    #print(baseName)
                
                parameters = [nested for nested in [fileManager.classes[baseName]] if nested['name']=="Parameters"][0]
            else:
                break

    return pot

if __name__=="__main__":
    
    bond_types = []
    unbound_types = []
    
    interactors_bond_types = []
    interactors_unbound_types = []
    
    bond_load = []
    unbound_load = []

    for bnd in POTENTIAL_BONDS:
        for f in os.listdir(POTENTIAL_PATH+bnd):
            if(f==bnd+".cuh"):
                continue
            print(bnd,f)
            filename = POTENTIAL_PATH+bnd+"/"+f
            
            potInfo = getBondPotentialInfoFromFile(filename)
    
            for name,pI in potInfo.items():
                param = pI["parameters"]
                bndTypeName = bnd+name+"Type"

                print("     ",name,bndTypeName,param)
    
                if(pI["template"] in TEMPLATES_NEEDED.keys()):
                    bndType = "        using "+bndTypeName+" = Potentials::"+bnd+"::"+name+"<typename "+TEMPLATES_NEEDED[pI["template"]]+">;"
                else:
                    bndType = "        using "+bndTypeName+" = Potentials::"+bnd+"::"+name+";"
                
                #print(bndType)
    
                interactorBndType = """
        using Interactor{}   = Interactor::BondedInteractor<{},
        Interactor::BondedInteractor_ns::BondProcessor<{}>,
        Interactor::BondedInteractor_ns::BondReaderFromFile<typename Base::Topology,{}>>;""".format(bndTypeName,
                                                                                                    bndTypeName,
                                                                                                    bndTypeName,
                                                                                                    bndTypeName)
                
                #bndPtr = "        std::shared_ptr<Interactor{}>   {};".format(name+"Type",name.lower())
                
                bndParam = ""
                for pname,ptype in param.items():
                    bndParam+="bndParam.{}=Miscellany::str2{}(entryInfo[i].param[\"{}\"],this->sys);\n                        ".format(pname,ptype,pname)
                
                bndParamOut = ""
                for pname,ptype in param.items():
                    bndParamOut+="this->sys->template log<System::MESSAGE>(\"[Generic] Added parameter: {} with value %s, to interactor: {}\",\
                                  \n                                                             \
 entryInfo[i].param[\"{}\"].c_str());\n                    ".format(pname,bnd+"::"+name,pname)
    
                bndLoad = """
                if(this->top->isEntryPresent(\"{}\",\"{}\")){{

                    auto entryInfo = this->top->getEntryInfo(\"{}\",\"{}\");

                    for(uint i=0;i<entryInfo.size();i++){{
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : {}");
                        
                        typename {}::Parameters bndParam;
    
                        {}
                        {}
                        std::shared_ptr<{}> bnd = std::make_shared<{}>(this->pd,
                                                                       bndParam);
                        
                        typename Interactor{}::Parameters interactorBndParameters;
                        
                        interactorBndParameters.bondName = entryInfo[i].label;
                        
                        std::string interactorName = entryInfo[i].alias;
                        
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<Interactor{}>(this->pg, 
                                                                           this->top, bnd,
                                                                           interactorBndParameters);
                    }}
    
                }}""".format(bnd,name,
                             bnd,name,
                             name,
                             bndTypeName,
                             bndParam,
                             bndParamOut,
                             bndTypeName,
                             bndTypeName,
                             bndTypeName,
                             bndTypeName)
    
    
                #bond_ptrs.append(bndPtr)
                bond_types.append(bndType)
                interactors_bond_types.append(interactorBndType)
                bond_load.append(bndLoad)
                #pdb.set_trace()
    
    for ubnd in POTENTIAL_UNBOUND:
        for f in os.listdir(POTENTIAL_PATH+ubnd):
            print(ubnd,f)
            filename = POTENTIAL_PATH+ubnd+"/"+f
            
            potInfo = getUnBoundPotentialInfoFromFile(filename)
    
            for name,pI in potInfo.items():
                param = pI["parameters"]
                ubndTypeName = ubnd+name+"Type"

                print("     ",name,ubndTypeName,param)
    
                if(pI["template"] in TEMPLATES_NEEDED.keys()):
                    print("TODO")
                    sys.exit(1)
                else:
                    ubndType = "        using "+ubndTypeName+" = Potentials::"+ubnd+"::"+name+"<typename Base::Topology>;"

                interactorUbndType = """
        using Interactor{}   = Interactor::PairInteractor<{},NeighbourList>;""".format(ubndTypeName,
                                                                                                       ubndTypeName)
                
                ubndParam = ""
                for pname,ptype in param.items():
                    if pname == "label":
                        ubndParam+="ubndParam.{}=entryInfo[i].label;\n                        ".format(pname,ubnd,name)
                    else:
                        ubndParam+="ubndParam.{}=Miscellany::str2{}(entryInfo[i].param[\"{}\"],this->sys);\n                        ".format(pname,ptype,pname)
                
                ubndParamOut = ""
                for pname,ptype in param.items():
                    if pname == "label":
                        ubndParamOut+="this->sys->template log<System::MESSAGE>(\"[Generic] Added parameter: {} with value %s, to interactor: {}\",\
                                  \n                                                             \
 entryInfo[i].label.c_str());\n                    ".format(pname,ubnd+"::"+name)
                    else:
                        ubndParamOut+="this->sys->template log<System::MESSAGE>(\"[Generic] Added parameter: {} with value %s, to interactor: {}\",\
                                  \n                                                             \
 entryInfo[i].param[\"{}\"].c_str());\n                        ".format(pname,ubnd+"::"+name,pname)
    
                ubndLoad = """
                if(this->top->isEntryPresent(\"{}\",\"{}\")){{
                    
                    if(!isNeighbourListInit){{
                        this->initNeighbourList(in);
                    }}
                    
                    auto entryInfo = this->top->getEntryInfo(\"{}\",\"{}\");

                    for(uint i=0;i<entryInfo.size();i++){{
                
                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Detected interactor : {}");
                        
                        typename {}::Parameters ubndParam;
    
                        {}
                        {}
                        std::shared_ptr<{}> ubnd = std::make_shared<{}>(this->pg,
                                                                        this->top,
                                                                        ubndParam);
                        
                        if(ubnd->getCutOffDst() >= this->nl->getCutOffVerlet()){{
                            this->sys->template log<System::CRITICAL>(\"[Generic] Error in potential {}, cutOffDst (%f) \"
                                                                      \"has to be smaller than VerletListDst (%f)\",
                                                                       ubnd->getCutOffDst(),this->nl->getCutOffVerlet());
                        }}
                        
                        this->nl->setCutOff(std::max(this->nl->getCutOff(),ubnd->getCutOffDst()));
                
                        typename Interactor{}::Parameters interactorUbndParameters;
                        
                        interactorUbndParameters.name                     = \"{}\";
                        interactorUbndParameters.pot                      = ubnd;
                        interactorUbndParameters.nl                       = this->nl;
                        interactorUbndParameters.conditionInteractionName = Miscellany::str2str(entryInfo[i].param["condition"],this->sys);
                        
                        std::string interactorName = entryInfo[i].alias;

                        this->sys->template log<System::MESSAGE>("[Generic] "
                                                                 "Adding interactor : %s",interactorName.c_str());

                        interactors[interactorName]=std::make_shared<Interactor{}>(this->pg, 
                                                                                   interactorUbndParameters);
                    }}
    
                }}""".format(ubnd,name,
                             ubnd,name,
                             name,
                             ubndTypeName,
                             ubndParam,
                             ubndParamOut,
                             ubndTypeName,
                             ubndTypeName,
                             ubnd+name,
                             ubndTypeName,
                             ubnd+name,
                             ubndTypeName)
    
    
                unbound_types.append(ubndType)
                interactors_unbound_types.append(interactorUbndType)
                unbound_load.append(ubndLoad)
                        
    fout=open(FORCE_FIELD_PATH+"Generic.cuh","w")
    
    fout.write("""
#ifndef __GENERIC__
#define __GENERIC__

namespace uammd{
namespace structured{ 
namespace forceField{
namespace Generic{

template<class Units_,
         class Types_,
         template <class Topology_> class Condition_>
class Generic : public ForceFieldBase<Units_,Types_>{
        
    protected:

        using Base = ForceFieldBase<Units_,Types_>;

        using Condition = Condition_<typename Base::Topology>;
        
        using NeighbourList = ConditionedVerletListSet<Condition>;
        
        std::shared_ptr<Condition> condition;
        std::shared_ptr<NeighbourList>    nl; 
            
        real VerletListDst;\n\n""")
    
    for bt in bond_types:
        fout.write(bt+"\n")
    for ubt in unbound_types:
        fout.write(ubt+"\n")
    fout.write("\n")

    for bti in interactors_bond_types:
        fout.write(bti+"\n")
    for ubti in interactors_unbound_types:
        fout.write(ubti+"\n")
    fout.write("\n")
    
    fout.write("""
        std::map<std::string,std::shared_ptr<uammd::Interactor>> interactors;\n""")
    
    fout.write("""
        std::map<std::string,std::shared_ptr<uammd::Interactor>> stoppedInteractors;\n""")
        
    fout.write("""
        bool isNeighbourListInit = false;

        void initNeighbourList(InputFile& in){

            this->sys->template log<System::MESSAGE>(\"[Generic] Initializing neighbour list ...\");
            
            VerletListDst = std::stof(in.getOption(\"VerletListDst\",InputFile::Required).str());

            condition = std::make_shared<Condition>(this->pd,this->top,in);
            
            typename NeighbourList::Parameters NeighbourListParam;
            
            NeighbourListParam.cutOff       = real(0.0);
            NeighbourListParam.cutOffVerlet = VerletListDst;

            nl = std::make_shared<NeighbourList>(this->pg,
                                                 condition,
                                                 NeighbourListParam);
            isNeighbourListInit = true;
        }\n""")

    fout.write("""
    public:\n""")

    fout.write("""
        Generic(std::shared_ptr<ParticleGroup> pg,
                InputFile&                     in):Base(pg,in){

                this->sys->template log<System::MESSAGE>(\"[Generic] Start\");
                """)
    
    for bl in bond_load:
        fout.write(bl+"\n")
    for ubl in unbound_load:
        fout.write(ubl+"\n")
    fout.write("\n")
                    
    fout.write("""
                if(isNeighbourListInit){
                    this->sys->template log<System::MESSAGE>(\"[Generic] Neighbour list cutoff: %f, Neighbour list verlet cut off: %f\",
                                                              this->nl->getCutOff(),this->nl->getCutOffVerlet());
                }
        }\n""")
    
    fout.write("""
            
        std::vector<std::string> getComponentsList(){
            std::vector<std::string> compList;
            
            for(auto inter : interactors){
                compList.push_back(inter.first);
            }
            return compList;
        }
            
        void sum(std::string component,Computables comp,cudaStream_t st){
            for(auto inter : interactors){
                if(inter.first == component){
                    inter.second->sum(comp,st);
                }
            }
            this->sys->template log<System::CRITICAL>("[Generic] Requested potential %s to sum. "
                                                        "But %s has not been added.",
                                                        component.c_str(),component.c_str());
        }
        
        void stop(std::string alias){
            for(auto inter : interactors){
                if(inter.first == alias){
                    stoppedInteractors[alias]=interactors[alias];
                    interactors.erase(alias);
                
                    this->sys->template log<System::MESSAGE>("[Generic] "
                                                             "Stopped interactor : %s",alias.c_str());
                    return;
                }
            }
            this->sys->template log<System::CRITICAL>("[Generic] An attempt has been made to stop the %s interactor,"
                                                        "but is not present in the list of active interactors.",
                                                        alias.c_str());
        }
        
        void resume(std::string alias){
            for(auto inter : stoppedInteractors){
                if(inter.first == alias){
                    interactors[alias]=stoppedInteractors[alias];
                    stoppedInteractors.erase(alias);

                    this->sys->template log<System::MESSAGE>("[Generic] "
                                                             "Resumed interactor : %s",alias.c_str());
                    return;
                }
            }
            this->sys->template log<System::CRITICAL>("[Generic] An attempt has been made to resume the %s interactor,"
                                                        "but is not present in the list of stopped interactors.",
                                                        alias.c_str());
        }
    
        void sum(Computables comp,cudaStream_t st) override {
            Base::sum(comp,st);
            for(auto inter : interactors){
                inter.second->sum(comp,st);
            }
        }
        
        void updateBox(Box box){
            Base::updateBox(box);
            for(auto inter : interactors){
                inter.second->updateBox(box);
            }
        }

        """
    )
    
    fout.write(
"""};

}}}}
#endif\n
    """)
    
    fout.close()
