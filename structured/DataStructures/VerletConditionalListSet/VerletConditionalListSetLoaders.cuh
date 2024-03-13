#ifndef __VERLET_CONDITIONAL_LIST_SET_LOADERS__
#define __VERLET_CONDITIONAL_LIST_SET_LOADERS__

namespace uammd{
namespace structured{
namespace VerletConditionalListSetLoaders{

    std::shared_ptr<uammd::structured::VerletConditionalListSet<conditions::all>>
    loadVerletConditionalListSet_all(std::shared_ptr<ExtendedSystem> sys,
                                     std::shared_ptr<GlobalData>           gd,
                                     std::shared_ptr<ParticleGroup>        pg,
                                     std::vector<std::string>      path){

        std::shared_ptr<ExtendedParticleData> pd = getExtendedParticleData(pg->getParticleData());

        DataEntry data = sys->getInput()->getDataEntry(path);
        std::shared_ptr<conditions::all>  cond =
        std::make_shared<conditions::all>(gd,pd,data);

        return std::make_shared<uammd::structured::VerletConditionalListSet<conditions::all>>(gd,pg,data,cond,path.back());
    }

    std::shared_ptr<uammd::structured::VerletConditionalListSet<conditions::intra_inter>>
    loadVerletConditionalListSet_intra_inter(std::shared_ptr<ExtendedSystem> sys,
                                             std::shared_ptr<GlobalData>           gd,
                                             std::shared_ptr<ParticleGroup>        pg,
                                             std::vector<std::string>      path){

        std::shared_ptr<ExtendedParticleData> pd = getExtendedParticleData(pg->getParticleData());

        DataEntry data = sys->getInput()->getDataEntry(path);
        std::shared_ptr<conditions::intra_inter>  cond =
        std::make_shared<conditions::intra_inter>(gd,pd,data);

        return std::make_shared<uammd::structured::VerletConditionalListSet<conditions::intra_inter>>(gd,pg,data,cond,path.back());
    }

    std::shared_ptr<uammd::structured::VerletConditionalListSet<conditions::nonExcluded>>
    loadVerletConditionalListSet_nonExcluded(std::shared_ptr<ExtendedSystem> sys,
                                             std::shared_ptr<GlobalData>           gd,
                                             std::shared_ptr<ParticleGroup>        pg,
                                             std::vector<std::string>      path){

        std::shared_ptr<ExtendedParticleData> pd = getExtendedParticleData(pg->getParticleData());

        DataEntry data = sys->getInput()->getDataEntry(path);
        std::shared_ptr<conditions::nonExcluded>  cond =
        std::make_shared<conditions::nonExcluded>(gd,pd,data);

        return std::make_shared<uammd::structured::VerletConditionalListSet<conditions::nonExcluded>>(gd,pg,data,cond,path.back());
    }

    std::shared_ptr<uammd::structured::VerletConditionalListSet<conditions::nonExclIntra_nonExclInter>>
    loadVerletConditionalListSet_nonExclIntra_nonExclInter(std::shared_ptr<ExtendedSystem> sys,
                                                           std::shared_ptr<GlobalData>           gd,
                                                           std::shared_ptr<ParticleGroup>        pg,
                                                           std::vector<std::string>      path){

        std::shared_ptr<ExtendedParticleData> pd = getExtendedParticleData(pg->getParticleData());

        DataEntry data = sys->getInput()->getDataEntry(path);
        std::shared_ptr<conditions::nonExclIntra_nonExclInter>  cond =
        std::make_shared<conditions::nonExclIntra_nonExclInter>(gd,pd,data);

        return std::make_shared<uammd::structured::VerletConditionalListSet<conditions::nonExclIntra_nonExclInter>>(gd,pg,data,cond,path.back());
    }

    std::shared_ptr<uammd::structured::VerletConditionalListSet<conditions::nonExclTypeGroup1Intra_nonExclTypeGroup2Intra_nonExclInter_nonExclNoGroup>>
    loadVerletConditionalListSet_nonExclTypeGroup1Intra_nonExclTypeGroup2Intra_nonExclInter_nonExclNoGroup(std::shared_ptr<ExtendedSystem> sys,
                                                                                                           std::shared_ptr<GlobalData>           gd,
                                                                                                           std::shared_ptr<ParticleGroup>        pg,
                                                                                                           std::vector<std::string>      path){

        std::shared_ptr<ExtendedParticleData> pd = getExtendedParticleData(pg->getParticleData());

        DataEntry data = sys->getInput()->getDataEntry(path);
        std::shared_ptr<conditions::nonExclTypeGroup1Intra_nonExclTypeGroup2Intra_nonExclInter_nonExclNoGroup>  cond =
        std::make_shared<conditions::nonExclTypeGroup1Intra_nonExclTypeGroup2Intra_nonExclInter_nonExclNoGroup>(gd,pd,data);

        return std::make_shared<uammd::structured::VerletConditionalListSet<conditions::nonExclTypeGroup1Intra_nonExclTypeGroup2Intra_nonExclInter_nonExclNoGroup>>(gd,pg,data,cond,path.back());
    }

    std::shared_ptr<uammd::structured::VerletConditionalListSet<conditions::nonExclIntra_nonExclInter_nonExclCharged>>
    loadVerletConditionalListSet_nonExclIntra_nonExclInter_nonExclCharged(std::shared_ptr<ExtendedSystem> sys,
                                                                          std::shared_ptr<GlobalData>           gd,
                                                                          std::shared_ptr<ParticleGroup>        pg,
                                                                          std::vector<std::string>      path){

        std::shared_ptr<ExtendedParticleData> pd = getExtendedParticleData(pg->getParticleData());

        DataEntry data = sys->getInput()->getDataEntry(path);
        std::shared_ptr<conditions::nonExclIntra_nonExclInter_nonExclCharged>  cond =
        std::make_shared<conditions::nonExclIntra_nonExclInter_nonExclCharged>(gd,pd,data);

        return std::make_shared<uammd::structured::VerletConditionalListSet<conditions::nonExclIntra_nonExclInter_nonExclCharged>>(gd,pg,data,cond,path.back());
    }

    std::shared_ptr<uammd::structured::VerletConditionalListSet<conditions::nonExclIntra_nonExclInter_nonExclInterCharged>>
    loadVerletConditionalListSet_nonExclIntra_nonExclInter_nonExclInterCharged(std::shared_ptr<ExtendedSystem> sys,
                                                                               std::shared_ptr<GlobalData>           gd,
                                                                               std::shared_ptr<ParticleGroup>        pg,
                                                                               std::vector<std::string>      path){

        std::shared_ptr<ExtendedParticleData> pd = getExtendedParticleData(pg->getParticleData());

        DataEntry data = sys->getInput()->getDataEntry(path);
        std::shared_ptr<conditions::nonExclIntra_nonExclInter_nonExclInterCharged>  cond =
        std::make_shared<conditions::nonExclIntra_nonExclInter_nonExclInterCharged>(gd,pd,data);

        return std::make_shared<uammd::structured::VerletConditionalListSet<conditions::nonExclIntra_nonExclInter_nonExclInterCharged>>(gd,pg,data,cond,path.back());
    }

    std::shared_ptr<uammd::structured::VerletConditionalListSet<conditions::interDifferentType>>
    loadVerletConditionalListSet_interDifferentType(std::shared_ptr<ExtendedSystem> sys,
                                                    std::shared_ptr<GlobalData>           gd,
                                                    std::shared_ptr<ParticleGroup>        pg,
                                                    std::vector<std::string>      path){

        std::shared_ptr<ExtendedParticleData> pd = getExtendedParticleData(pg->getParticleData());

        DataEntry data = sys->getInput()->getDataEntry(path);
        std::shared_ptr<conditions::interDifferentType>  cond =
        std::make_shared<conditions::interDifferentType>(gd,pd,data);

        return std::make_shared<uammd::structured::VerletConditionalListSet<conditions::interDifferentType>>(gd,pg,data,cond,path.back());
    }

    std::shared_ptr<uammd::structured::VerletConditionalListSetBase>
    loadVerletConditionalListSet(std::shared_ptr<ExtendedSystem> sys,
                                 std::shared_ptr<GlobalData>    gd,
                                 std::map<std::string,std::shared_ptr<ParticleGroup>> groups,
                                 std::vector<std::string>       path){


        DataEntry data = sys->getInput()->getDataEntry(path);
        //Check data type is VerletConditionalListSet
        if(data.getType() != "VerletConditionalListSet"){
            System::log<System::CRITICAL>("[VerletConditionalListSetLoader] (%s) Data type is not VerletConditionalListSet!",path.back().c_str());
        }

        std::shared_ptr<ParticleGroup> pg = GroupUtils::getParticleGroupFromGroupsList(groups,data,"All");

        //////////////////////////////////////////////////////////////////////////////

        std::shared_ptr<uammd::structured::VerletConditionalListSetBase> vCondListSet;

        std::string condition = data.getSubType();
        bool found=false;

        if("all" == condition){
            System::log<System::MESSAGE>("[VerletConditionalListSetLoader] (%s) Detected VerletConditionalListSet all",path.back().c_str());
            vCondListSet = VerletConditionalListSetLoaders::loadVerletConditionalListSet_all(sys,gd,pg,path);
            found = true;
        }
        if("intra_inter" == condition){
            System::log<System::MESSAGE>("[VerletConditionalListSetLoader] (%s) Detected VerletConditionalListSet intra_inter",path.back().c_str());
            vCondListSet = VerletConditionalListSetLoaders::loadVerletConditionalListSet_intra_inter(sys,gd,pg,path);
            found = true;
        }
        if("nonExcluded" == condition){
            System::log<System::MESSAGE>("[VerletConditionalListSetLoader] (%s) Detected VerletConditionalListSet nonExcluded",path.back().c_str());
            vCondListSet = VerletConditionalListSetLoaders::loadVerletConditionalListSet_nonExcluded(sys,gd,pg,path);
            found = true;
        }
        if("nonExclIntra_nonExclInter" == condition){
            System::log<System::MESSAGE>("[VerletConditionalListSetLoader] (%s) Detected VerletConditionalListSet nonExclIntra_nonExclInter",path.back().c_str());
            vCondListSet = VerletConditionalListSetLoaders::loadVerletConditionalListSet_nonExclIntra_nonExclInter(sys,gd,pg,path);
            found = true;
        }
        if("nonExclTypeGroup1Intra_nonExclTypeGroup2Intra_nonExclInter_nonExclNoGroup" == condition){
            System::log<System::MESSAGE>("[VerletConditionalListSetLoader] (%s) Detected VerletConditionalListSet nonExclTypeGroup1Intra_nonExclTypeGroup2Intra_nonExclInter_nonExclNoGroup",path.back().c_str());
            vCondListSet = VerletConditionalListSetLoaders::loadVerletConditionalListSet_nonExclTypeGroup1Intra_nonExclTypeGroup2Intra_nonExclInter_nonExclNoGroup(sys,gd,pg,path);
            found = true;
        }
        if("nonExclIntra_nonExclInter_nonExclCharged" == condition){
            System::log<System::MESSAGE>("[VerletConditionalListSetLoader] (%s) Detected VerletConditionalListSet nonExclIntra_nonExclInter_nonExclCharged",path.back().c_str());
            vCondListSet = VerletConditionalListSetLoaders::loadVerletConditionalListSet_nonExclIntra_nonExclInter_nonExclCharged(sys,gd,pg,path);
            found = true;
        }
        if("nonExclIntra_nonExclInter_nonExclInterCharged" == condition){
            System::log<System::MESSAGE>("[VerletConditionalListSetLoader] (%s) Detected VerletConditionalListSet nonExclIntra_nonExclInter_nonExclInterCharged",path.back().c_str());
            vCondListSet = VerletConditionalListSetLoaders::loadVerletConditionalListSet_nonExclIntra_nonExclInter_nonExclInterCharged(sys,gd,pg,path);
            found = true;
        }
        if("interDifferentType" == condition){
            System::log<System::MESSAGE>("[VerletConditionalListSetLoader] (%s) Detected VerletConditionalListSet interDifferentType",path.back().c_str());
            vCondListSet = VerletConditionalListSetLoaders::loadVerletConditionalListSet_interDifferentType(sys,gd,pg,path);
            found = true;
        }

        if(not found){
            std::string available = "all, intra_inter, nonExcluded, nonExclIntra_nonExclInter, nonExclTypeGroup1Intra_nonExclTypeGroup2Intra_nonExclInter_nonExclNoGroup, nonExclIntra_nonExclInter_nonExclCharged, nonExclIntra_nonExclInter_nonExclInterCharged, interDifferentType";
            System::log<System::CRITICAL>("[VerletConditionalListSetLoader] (%s) Condition '%s' not available. Available conditions are: %s",path.back().c_str(),condition.c_str(),available.c_str());
        }

        return vCondListSet;

    }
}}}
#endif
