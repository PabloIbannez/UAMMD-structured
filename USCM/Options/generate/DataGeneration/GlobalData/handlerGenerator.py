import json

def generateHandler(typ,componentsPath,particleDataComponentPath,outputFolder,outputName,genGetter=True,genSetter=True):

    # First letter uppercase
    Typ = typ[0].upper() + typ[1:]

    with open(componentsPath,"r") as f:
        variableParameters = json.load(f)
        variableParameters = variableParameters[particleDataComponentPath]

    with open(outputFolder+"/"+outputName,"w") as fout:

        whitespace = "        "

        getters = ""
        setters = ""
        for i,variable in enumerate(variableParameters):
            name,getter,vtyp = variable

            if genGetter:
                getters += whitespace
                getters += f"""virtual {vtyp} get{getter}(){{
            System::log<System::CRITICAL>("[{Typ}] {getter} not defined for {typ} \\"%s\\".",
                                          subType.c_str());
        }}\n"""

            if genSetter:
                setters += whitespace
                setters += f"""virtual void set{getter}({vtyp} value){{
            System::log<System::CRITICAL>("[{Typ}] {getter} not defined for {typ} \\"%s\\".",
                                          subType.c_str());
        }}\n"""

        update = ""
        if genSetter:
            update = whitespace + f"""virtual void updateDataEntry(DataEntry data) = 0;\n"""

        template = f"""
#ifndef __{typ.upper()}_HANDLER__
#define __{typ.upper()}_HANDLER__

namespace uammd{{
namespace structured{{
namespace {Typ}{{

class {Typ}Handler{{

    protected:

        std::string subType;

    public:

        {Typ}Handler(DataEntry& data){{
            subType = data.getSubType();
        }}

        std::string getSubType() const {{
            return subType;
        }}

        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wreturn-type"

{getters}
{setters}

        #pragma GCC diagnostic pop

{update}

}};

}}}}}}
#endif
"""

        fout.write(template)
