import json
import logging

def getLoader(tpy,comp):

    #First letter to upper case
    Tpy = tpy[0].upper() + tpy[1:]

    tpyLoader=""
    for c in comp:
        cType,cSubType,file = c

        tpyLoader += f"""
        if("{cType}" == {tpy}Type and "{cSubType}" == {tpy}SubType){{
            System::log<System::MESSAGE>("[{Tpy}Loader] (%s) Detected {cSubType} {tpy}",path.back().c_str());
            {tpy} = std::make_shared<{cType}::{cSubType}>(data);
            found = true;
        }}"""


    loaderTemplate=f"""
    std::shared_ptr<typename {Tpy}::{Tpy}Handler>
    inline
    load{Tpy}(std::shared_ptr<ExtendedSystem> sys,
              std::vector<std::string>       path){{

        DataEntry data = sys->getInput()->getDataEntry(path);

        std::string {tpy}Type    = data.getType();
        std::string {tpy}SubType = data.getSubType();

        std::shared_ptr<typename {Tpy}::{Tpy}Handler> {tpy};
        bool found = false;

        {tpyLoader}

        if(not found){{
            System::log<System::CRITICAL>("[{Tpy}Loader] (%s) Could not find {tpy} %s::%s",path.back().c_str(),
                                           {tpy}Type.c_str(),{tpy}SubType.c_str());
        }}

        return {tpy};

    }}\n
    """
    return loaderTemplate

    loaderTemplate=f"""
    """
    return loaderTemplate

def generateGenericLoader(typ,componentsPath,componentPath,uammdStructuredBaseFolder):

    logger = logging.getLogger("USCM")

    # First letter to upper case
    Typ = typ[0].upper() + typ[1:]

    output = "/".join(componentPath)+"/"+componentPath[-1]+"Loaders.cuh"

    with open(componentsPath,"r") as f:
        comp = json.load(f)
        for path in componentPath:
            comp = comp[path]

    with open(uammdStructuredBaseFolder+"/"+output,"w") as fout:
        fout.write(f"#ifndef __{typ.upper()}_LOADER__\n")
        fout.write(f"#define __{typ.upper()}_LOADER__\n")

        fout.write( "namespace uammd{\n")
        fout.write( "namespace structured{\n")
        fout.write(f"namespace {Typ}Loader{{\n")

        logger.debug(f"Generating generic {typ} loader")

        loader  = getLoader(typ,comp)

        fout.write(loader)
        fout.write("}}}\n")
        fout.write("#endif\n")

    return f"#include\"{output}\""

