import json
import logging

def getIsInteractor(interactors):

    #Iterate over all potentials (merging all lists) and generate the isInteractor function
    isInteractor = ""
    for pot in interactors:
        for subClass,potType,potSubType,f in pot:
            isInteractor += f"""
        if(\"{potType}\" == potType and \"{potSubType}\" == potSubType){{
            return true;
        }}
            """

    isInteractorTemplate=f"""
    bool isInteractorAvailable(std::shared_ptr<ExtendedSystem> sys,
                               std::vector<std::string>       path){{

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string potType    = data.getType();
        std::string potSubType = data.getSubType();
        {isInteractor}
        return false;

    }}\n
    """
    return isInteractorTemplate


def getGenericLoader(basics,nls,pp):

    loaders =""
    for pot in basics:
        for subClass,potType,potSubType,f in pot:
            loader = subClass
            name   = f"load{potType}_{potSubType}"

            loaderTemplate = f"""
        if("{potType}" == potType and "{potSubType}" == potSubType){{
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected {potType}::{potSubType} potential",path.back().c_str());

            std::shared_ptr<{potType}::{potSubType}> pot =
            std::make_shared<{potType}::{potSubType}>(gd,pg,data);

            interactor = std::make_shared<typename Interactor::{subClass}Interactor<{potType}::{potSubType}>>(gd,pg,data,pot,path.back());
            found = true;
        }}"""
            loaders += loaderTemplate

    for pot in nls:
        for subClass,potType,potSubType,f in pot:
            loader = subClass
            name   = f"load{potType}_{potSubType}"

            loaderTemplate = f"""
        if("{potType}" == potType and "{potSubType}" == potSubType){{
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected {potType}::{potSubType} potential",path.back().c_str());

            std::shared_ptr<{potType}::{potSubType}> pot =
            std::make_shared<{potType}::{potSubType}>(gd,pg,data);

            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);

            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));

            interactor = std::make_shared<typename Interactor::{subClass}Interactor<{potType}::{potSubType},VerletConditionalListSetBase>>(gd,pg,data,pot,nl,path.back());

            found = true;
        }}"""
            loaders += loaderTemplate

    ppLoader=""
    for pot in pp:
        for subClass,potType,potSubType,f in pot:
            loader = subClass
            name   = f"load{potType}_{potSubType}"

            loaderTemplate = f"""
        if("{potType}" == potType and "{potSubType}" == potSubType){{
            System::log<System::MESSAGE>("[GenericLoader] (%s) Detected {potType}::{potSubType} potential",path.back().c_str());
            interactor = std::make_shared<typename Interactor::{potType}::{potSubType}>(gd,pg,path,path.back());
            found = true;
        }}"""
            loaders += loaderTemplate

    loaderTemplate=f"""
    std::shared_ptr<typename uammd::Interactor>
    loadGeneric(std::shared_ptr<ExtendedSystem> sys,
                std::shared_ptr<GlobalData>     gd,
                std::map<std::string,std::shared_ptr<ParticleGroup>>& groups,
                std::map<std::string,std::shared_ptr<VerletConditionalListSetBase>>& nls,
                std::vector<std::string>       path){{

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::shared_ptr<ParticleGroup> pg = GroupUtils::getParticleGroupFromGroupsList(groups,data,"All");

        std::string potType    = data.getType();
        std::string potSubType = data.getSubType();

        std::shared_ptr<typename uammd::Interactor> interactor;
        bool found = false;
{loaders}

        if(not found){{
            System::log<System::CRITICAL>("[GenericLoader] (%s) Potential of class %s and subclass %s not found",
                                           path.back().c_str(),data.getType().c_str(),data.getSubType().c_str());
        }}

        return interactor;

    }}\n
    """
    return loaderTemplate



def generateGenericPotentialLoader(componentsPath,potsComponentPath,uammdStructuredBaseFolder):

    logger = logging.getLogger("USCM")

    output = "/".join(potsComponentPath)+"/GenericPotentialLoader.cuh"

    with open(componentsPath,"r") as f:
        potsComp = json.load(f)
        for path in potsComponentPath:
            potsComp = potsComp[path]

    ignoredKeys = ["Patches"]
    nlsKeys     = ["Pair"]
    ppKeys      = ["PatchyParticles"]

    basics = []
    nls    = []
    pp     = []
    for key in potsComp.keys():
        if key not in ignoredKeys:
            tmp = [[key]+p for p in potsComp[key]]
            if(key in nlsKeys):
                nls.append(tmp)
            elif(key in ppKeys):
                pp.append(tmp)
            else:
                basics.append(tmp)

    with open(uammdStructuredBaseFolder+"/"+output,"w") as fout:
        fout.write("#ifndef __GENERIC_LOADER__\n")
        fout.write("#define __GENERIC_LOADER__\n")

        fout.write("namespace uammd{\n")
        fout.write("namespace structured{\n")
        fout.write("namespace Potentials{\n")
        fout.write("namespace GenericLoader{\n")

        logger.debug("Generating isInteractor function")
        logger.debug("Generating generic potential loader")

        isInteractor = getIsInteractor(basics+nls+pp)
        loader       = getGenericLoader(basics,nls,pp)

        fout.write(isInteractor)
        fout.write(loader)

        fout.write("}}}}\n")
        fout.write("#endif\n")

    return f"#include\"{output}\""
