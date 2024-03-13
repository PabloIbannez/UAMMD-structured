import json

def generateParticleData(componentsPath,particleDataComponentPath,outputFolder,outputName):

    with open(componentsPath,"r") as f:
        variableParameters = json.load(f)
        variableParameters = variableParameters[particleDataComponentPath]

    with open(outputFolder+"/"+outputName,"w") as fout:

        whitespace = " "*len("#define EXTRA_PARTICLE_PROPERTIES")

        declarations = ""
        for i,variable in enumerate(variableParameters):
            name,getter,tpy = variable
            if i == 0:
                declarations += f"(({getter},{name},{tpy}))\\\n"
            else:
                declarations += whitespace+f" (({getter},{name},{tpy}))\\\n"

        template = f"""
#ifndef UAMMD_EXTENSIONS_PREAMBLE_H
#define UAMMD_EXTENSIONS_PREAMBLE_H

#define EXTRA_PARTICLE_PROPERTIES {declarations}



#define EXTRA_COMPUTABLES (magneticField)(transitionProbability)(hessian)(lambdaDerivative)(pairwiseForce)

#endif
"""

        fout.write(template)


