import sys
import logging

from colorama import init, Fore, Back, Style
init()

import json

def getIsInteractor(interactors):

    #Iterate over all potentials (merging all lists) and generate the isInteractor function
    isInteractor = ""
    for potType,potSubType in interactors:
        isInteractor += f"""
        if(\"{potType}\" == potType and \"{potSubType}\" == potSubType){{
            return true;
        }}
            """

    isInteractorTemplate=f"""
    bool isPatchyParticlesInteractorAvailable(std::shared_ptr<ExtendedSystem> sys,
                                              std::vector<std::string>       path){{

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string potType    = data.getType();
        std::string potSubType = data.getSubType();
        {isInteractor}
        return false;

    }}\n
    """
    return isInteractorTemplate

def getGenericLoader(bounds,nonbonded):

    boundsLoader=""
    for potType,potSubType in bounds:
        loader = potType.strip("0123456789")

        loaderTemplate = f"""
        if("{potType}" == potType and "{potSubType}" == potSubType){{
            System::log<System::MESSAGE>("[GenericPatchyParticlesLoader] (%s) Detected {potType}::{potSubType} potential",path.back().c_str());
            std::shared_ptr<Potentials::SurfacePatches::{potSubType}> pot = std::make_shared<Potentials::SurfacePatches::{potSubType}>(gd,pg,patchesGd,patchesPg,data);
            interactor = std::make_shared<typename Interactor::SingleInteractor<Potentials::SurfacePatches::{potSubType}>>(patchesGd,patchesPg,data,pot,path.back());
            found = true;
        }}"""
        boundsLoader += loaderTemplate

    nonbondedLoader=""
    for potType,potSubType in nonbonded:
        loader = potType.strip("0123456789")

        loaderTemplate = f"""
        if("{potType}" == potType and "{potSubType}" == potSubType){{
            System::log<System::MESSAGE>("[GenericPatchyParticlesLoader] (%s) Detected {potType}::{potSubType} potential",path.back().c_str());
            std::shared_ptr<VerletConditionalListSetBase> nl = VerletConditionalListSetUtils::getNeighbourListFromNeighbourListsList(nls,data);
            std::shared_ptr<Potentials::NonBondedPatches::{potSubType}> pot = std::make_shared<Potentials::NonBondedPatches::{potSubType}>(gd,pg,patchesGd,patchesPg,data);
            nl->setCutOff(std::max(nl->getCutOff(),pot->getCutOff()));
            interactor = std::make_shared<typename Interactor::PairInteractor<Potentials::NonBondedPatches::{potSubType},VerletConditionalListSetBase>>(patchesGd,patchesPg,data,pot,nl,path.back());
            found = true;
        }}"""
        nonbondedLoader += loaderTemplate

    loaderTemplate=f"""
    std::shared_ptr<typename uammd::Interactor>
    loadGenericPatchyParticles(std::shared_ptr<ExtendedSystem> sys,
                               std::shared_ptr<GlobalData>  gd,std::shared_ptr<ParticleGroup>  pg,
                               std::shared_ptr<GlobalData>  patchesGd,std::map<std::string,std::shared_ptr<ParticleGroup>>& groups,
                               std::map<std::string,std::shared_ptr<VerletConditionalListSetBase>>& nls,
                               std::vector<std::string>       path){{

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::shared_ptr<ParticleGroup> patchesPg = GroupUtils::getParticleGroupFromGroupsList(groups,data,"All");

        std::string potType    = data.getType();
        std::string potSubType = data.getSubType();

        std::shared_ptr<typename uammd::Interactor> interactor;
        bool found = false;
{nonbondedLoader}
{boundsLoader}

        if(not found){{
            System::log<System::CRITICAL>("[GenericPatchyParticles] (%s) Potential of type %s and subType %s not found",
                                           path.back().c_str(),data.getType().c_str(),data.getSubType().c_str());
        }}

        return interactor;

    }}\n
    """
    return loaderTemplate


def generateGenericPatchesPotentialLoader(componentsPath,genPatPotComponentPath,uammdStructuredBaseFolder):

    logger = logging.getLogger("USCM")

    output = "/".join(genPatPotComponentPath)+"/GenericPatchesPotentialLoader.cuh"

    with open(componentsPath,"r") as f:
        genPatPotComp = json.load(f)
        for path in genPatPotComponentPath:
            genPatPotComp = genPatPotComp[path]

    bounds    = []
    nonbonded = []
    for cls,subCls,file in genPatPotComp:
        if   (cls == "SurfacePatches"):
            bounds.append([cls,subCls])
        elif (cls == "NonBondedPatches"):
            nonbonded.append([cls,subCls])
        else:
            logger.error(f"Unknown potential class: {cls}")
            sys.exit(1)

    with open(uammdStructuredBaseFolder+"/"+output,"w") as fout:
        fout.write("#ifndef __GENERIC_PATCHY_PARTICLES_POTENTIALS_LOADER__\n")
        fout.write("#define __GENERIC_PATCHY_PARTICLES_POTENTIALS_LOADER__\n")

        fout.write("namespace uammd{\n")
        fout.write("namespace structured{\n")
        fout.write("namespace Potentials{\n")
        fout.write("namespace GenericPatchyParticlesLoader{\n")

        logger.debug("Generating isInteractor function")
        logger.debug("Generating generic patchy loader")

        isInteractor = getIsInteractor(bounds+nonbonded)
        loader       = getGenericLoader(bounds,nonbonded)

        fout.write(isInteractor)
        fout.write(loader)

        fout.write("}}}}\n")
        fout.write("#endif\n")

    return f"#include\"{output}\""


