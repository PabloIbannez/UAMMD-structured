optGeneric_template = '''
nSteps               {nSteps}
nStepsInfoInterval   {nStepsInterval}
nStepsWriteInterval  {nStepsWrite}
nStepsBackupInterval {nStepsBackupInterval}

outPutFilePath {outPutFilePath}
outPutFormat   {outPutFormat}

boxSize {boxX} {boxY} {boxZ}
T {temperature}

h {steepestDecentScale}

nStepsSteepestDescent                 {nStepsSteepestDescent}
nStepsSteepestDescentProgressInterval {nStepsSteepestDescentProgressInterval}
maxObjectiveForce                     {maxObjectiveForce}

dt {timeStep}
frictionConstant {frictionConstant}

cutOffDst     {cutOffDst}
VerletListDst {VerletListDst}
    
inputCoordPath    {inputCoordFile}
inputTopologyPath {inputTopFile}'''

optGenericWithElec_template = optGeneric_template+'''
dielectricConstant {dielectricConstant}
debyeLength {debyeLength}
'''

optGenericWithSurface_template = optGenericWithElec_template+'''

epsilonSurf {epsilonSurf}
sigmaSurf   {sigmaSurf}

surfacePosition {surfacePosition}
'''

optGenericClash_template = optGeneric_template+'''

lambda {lambd}
gamma  {gamma}
'''

optGenericClashWithCompression_template = optGenericClash_template+'''

initialSphereRadius {initialSphereRadius}
minimalSphereRadius {minimalSphereRadius}

compressionVelocity {compressionVelocity}
'''

optAFM_template = optGenericWithSurface_template+'''

frictionConstantTip {frictionConstantTip}

initialTipSampleDst {initialTipSampleDst}
descentVelocity     {descentVelocity}

minimalChipHeight {minimalChipHeight}

Mtip {Mtip}
Rtip {Rtip}
Kxytip {Kxytip}
Ktip {Ktip}

epsilonTip {epsilonTip}
sigmaTip {sigmaTip}

Atip {Atip}
Btip {Btip}

epsilonTipSurf {epsilonTipSurf}
sigmaTipSurf   {sigmaTipSurf}

ATipSurf {ATipSurf}
BTipSurf {BTipSurf}

nStepsIndentMeasure {nStepsIndentMeasure}
outputIndentationMeasureFilePath {outputIndentationMeasureFilePath}
'''

optGenericUmbrella_template = optGenericWithElec_template+'''

umbrellaK    {umbrellaK}
umbrellaInit {umbrellaInit}
umbrellaEnd  {umbrellaEnd}
umbrellaWindowsNumber {umbrellaWindowsNumber}
umbrellaCopies        {umbrellaCopies}

nStepsUmbrellaMeasure         {nStepsUmbrellaMeasure}
outputUmbrellaMeasureFilePath {outputUmbrellaMeasureFilePath}
'''

def writeOptionsGeneric(path:str,
                        nSteps:int,
                        nStepsInterval:int,
                        nStepsWrite:int,
                        nStepsBackupInterval:int,
                        outPutFilePath:str,
                        outPutFormat:str,
                        boxX:float,boxY:float,boxZ:float,
                        temperature:float,
                        steepestDecentScale:float,
                        nStepsSteepestDescent:int,
                        nStepsSteepestDescentProgressInterval:int,
                        maxObjectiveForce:float,
                        timeStep:float,
                        frictionConstant:float,
                        cutOffDst:float,
                        VerletListDst:float,
                        dielectricConstant:float,
                        debyeLength:float,
                        inputCoordFile:str,
                        inputTopFile:str):

    opt = optGenericWithElec_template.format(nSteps=nSteps,
                                             nStepsInterval=nStepsInterval,
                                             nStepsWrite=nStepsWrite,
                                             nStepsBackupInterval=nStepsBackupInterval,
                                             outPutFilePath=outPutFilePath,
                                             outPutFormat=outPutFormat,
                                             boxX=boxX,boxY=boxY,boxZ=boxZ,
                                             temperature=temperature,
                                             steepestDecentScale=steepestDecentScale,
                                             nStepsSteepestDescent=nStepsSteepestDescent,
                                             nStepsSteepestDescentProgressInterval=nStepsSteepestDescentProgressInterval,
                                             maxObjectiveForce=maxObjectiveForce,
                                             timeStep=timeStep,
                                             frictionConstant=frictionConstant,
                                             cutOffDst=cutOffDst,
                                             VerletListDst=VerletListDst,
                                             dielectricConstant=dielectricConstant,
                                             debyeLength=debyeLength,
                                             inputCoordFile=inputCoordFile,
                                             inputTopFile=inputTopFile)
    
    with open(path,"w") as f:
        f.write(opt)

def writeOptionsGenericFromDict(path,optionsDict):
    
    nSteps=optionsDict["nSteps"]
    nStepsInterval=optionsDict["nStepsInterval"]
    nStepsWrite=optionsDict["nStepsWrite"]
    nStepsBackupInterval=optionsDict["nStepsBackupInterval"]
    outPutFilePath=optionsDict["outPutFilePath"]
    outPutFormat=optionsDict["outPutFormat"]
    boxX=optionsDict["boxX"]
    boxY=optionsDict["boxY"]
    boxZ=optionsDict["boxZ"]
    temperature=optionsDict["temperature"]
    steepestDecentScale=optionsDict["steepestDecentScale"]
    nStepsSteepestDescent=optionsDict["nStepsSteepestDescent"]
    nStepsSteepestDescentProgressInterval=optionsDict["nStepsSteepestDescentProgressInterval"]
    maxObjectiveForce=optionsDict["maxObjectiveForce"]
    timeStep=optionsDict["timeStep"]
    frictionConstant=optionsDict["frictionConstant"]
    cutOffDst=optionsDict["cutOffDst"]
    VerletListDst=optionsDict["VerletListDst"]
    dielectricConstant=optionsDict["dielectricConstant"]
    debyeLength=optionsDict["debyeLength"]
    inputCoordFile=optionsDict["inputCoordFile"]
    inputTopFile=optionsDict["inputTopFile"]

    writeOptionsGeneric(path,
                        nSteps,
                        nStepsInterval,
                        nStepsWrite,
                        nStepsBackupInterval,
                        outPutFilePath,
                        outPutFormat,
                        boxX,boxY,boxZ,
                        temperature,
                        steepestDecentScale,
                        nStepsSteepestDescent,
                        nStepsSteepestDescentProgressInterval,
                        maxObjectiveForce,
                        timeStep,
                        frictionConstant,
                        cutOffDst,
                        VerletListDst,
                        dielectricConstant,
                        debyeLength,
                        inputCoordFile,
                        inputTopFile)

def writeOptionsGenericWithSurface(path:str,
                                   nSteps:int,
                                   nStepsInterval:int,
                                   nStepsWrite:int,
                                   nStepsBackupInterval:int,
                                   outPutFilePath:str,
                                   outPutFormat:str,
                                   boxX:float,boxY:float,boxZ:float,
                                   temperature:float,
                                   steepestDecentScale:float,
                                   nStepsSteepestDescent:int,
                                   nStepsSteepestDescentProgressInterval:int,
                                   maxObjectiveForce:float,
                                   timeStep:float,
                                   frictionConstant:float,
                                   cutOffDst:float,
                                   VerletListDst:float,
                                   dielectricConstant:float,
                                   debyeLength:float,
                                   inputCoordFile:str,
                                   inputTopFile:str,
                                   epsilonSurf:float,
                                   sigmaSurf:float,
                                   surfacePosition:float):

    opt = optGenericWithSurface_template.format(nSteps=nSteps,
                                                nStepsInterval=nStepsInterval,
                                                nStepsWrite=nStepsWrite,
                                                nStepsBackupInterval=nStepsBackupInterval,
                                                outPutFilePath=outPutFilePath,
                                                outPutFormat=outPutFormat,
                                                boxX=boxX,boxY=boxY,boxZ=boxZ,
                                                temperature=temperature,
                                                steepestDecentScale=steepestDecentScale,
                                                nStepsSteepestDescent=nStepsSteepestDescent,
                                                nStepsSteepestDescentProgressInterval=nStepsSteepestDescentProgressInterval,
                                                maxObjectiveForce=maxObjectiveForce,
                                                timeStep=timeStep,
                                                frictionConstant=frictionConstant,
                                                cutOffDst=cutOffDst,
                                                VerletListDst=VerletListDst,
                                                dielectricConstant=dielectricConstant,
                                                debyeLength=debyeLength,
                                                inputCoordFile=inputCoordFile,
                                                inputTopFile=inputTopFile,
                                                epsilonSurf=epsilonSurf,
                                                sigmaSurf=sigmaSurf,
                                                surfacePosition=surfacePosition)
    
    with open(path,"w") as f:
        f.write(opt)

def writeOptionsGenericWithSurfaceFromDict(path,optionsDict):

    nSteps=optionsDict["nSteps"]
    nStepsInterval=optionsDict["nStepsInterval"]
    nStepsWrite=optionsDict["nStepsWrite"]
    nStepsBackupInterval=optionsDict["nStepsBackupInterval"]
    outPutFilePath=optionsDict["outPutFilePath"]
    outPutFormat=optionsDict["outPutFormat"]
    boxX=optionsDict["boxX"]
    boxY=optionsDict["boxY"]
    boxZ=optionsDict["boxZ"]
    temperature=optionsDict["temperature"]
    steepestDecentScale=optionsDict["steepestDecentScale"]
    nStepsSteepestDescent=optionsDict["nStepsSteepestDescent"]
    nStepsSteepestDescentProgressInterval=optionsDict["nStepsSteepestDescentProgressInterval"]
    maxObjectiveForce=optionsDict["maxObjectiveForce"]
    timeStep=optionsDict["timeStep"]
    frictionConstant=optionsDict["frictionConstant"]
    cutOffDst=optionsDict["cutOffDst"]
    VerletListDst=optionsDict["VerletListDst"]
    inputCoordFile=optionsDict["inputCoordFile"]
    dielectricConstant=optionsDict["dielectricConstant"]
    debyeLength=optionsDict["debyeLength"]
    inputTopFile=optionsDict["inputTopFile"]
    epsilonSurf=optionsDict["epsilonSurf"]
    sigmaSurf=optionsDict["sigmaSurf"]
    surfacePosition=optionsDict["surfacePosition"]

    writeOptionsGenericWithSurface(path,
                                   nSteps,
                                   nStepsInterval,
                                   nStepsWrite,
                                   nStepsBackupInterval,
                                   outPutFilePath,
                                   outPutFormat,
                                   boxX,boxY,boxZ,
                                   temperature,
                                   steepestDecentScale,
                                   nStepsSteepestDescent,
                                   nStepsSteepestDescentProgressInterval,
                                   maxObjectiveForce,
                                   timeStep,
                                   frictionConstant,
                                   cutOffDst,
                                   VerletListDst,
                                   dielectricConstant,
                                   debyeLength,
                                   inputCoordFile,
                                   inputTopFile,
                                   epsilonSurf,
                                   sigmaSurf,                                 
                                   surfacePosition)

def writeOptionsGenericClash(path:str,
                             nSteps:int,
                             nStepsInterval:int,
                             nStepsWrite:int,
                             nStepsBackupInterval:int,
                             outPutFilePath:str,
                             outPutFormat:str,
                             boxX:float,boxY:float,boxZ:float,
                             temperature:float,
                             steepestDecentScale:float,
                             nStepsSteepestDescent:int,
                             nStepsSteepestDescentProgressInterval:int,
                             maxObjectiveForce:float,
                             timeStep:float,
                             frictionConstant:float,
                             cutOffDst:float,
                             VerletListDst:float,
                             inputCoordFile:str,
                             inputTopFile:str,
                             lambd:float,
                             gamma:float):

    opt = optGenericClash_template.format(nSteps=nSteps,
                                          nStepsInterval=nStepsInterval,
                                          nStepsWrite=nStepsWrite,
                                          nStepsBackupInterval=nStepsBackupInterval,
                                          outPutFilePath=outPutFilePath,
                                          outPutFormat=outPutFormat,
                                          boxX=boxX,boxY=boxY,boxZ=boxZ,
                                          temperature=temperature,
                                          steepestDecentScale=steepestDecentScale,
                                          nStepsSteepestDescent=nStepsSteepestDescent,
                                          nStepsSteepestDescentProgressInterval=nStepsSteepestDescentProgressInterval,
                                          maxObjectiveForce=maxObjectiveForce,
                                          timeStep=timeStep,
                                          frictionConstant=frictionConstant,
                                          cutOffDst=cutOffDst,
                                          VerletListDst=VerletListDst,
                                          inputCoordFile=inputCoordFile,
                                          inputTopFile=inputTopFile,
                                          lambd=lambd,
                                          gamma=gamma)
    
    with open(path,"w") as f:
        f.write(opt)

def writeOptionsGenericClashFromDict(path,optionsDict):
    
    nSteps=optionsDict["nSteps"]
    nStepsInterval=optionsDict["nStepsInterval"]
    nStepsWrite=optionsDict["nStepsWrite"]
    nStepsBackupInterval=optionsDict["nStepsBackupInterval"]
    outPutFilePath=optionsDict["outPutFilePath"]
    outPutFormat=optionsDict["outPutFormat"]
    boxX=optionsDict["boxX"]
    boxY=optionsDict["boxY"]
    boxZ=optionsDict["boxZ"]
    temperature=optionsDict["temperature"]
    steepestDecentScale=optionsDict["steepestDecentScale"]
    nStepsSteepestDescent=optionsDict["nStepsSteepestDescent"]
    nStepsSteepestDescentProgressInterval=optionsDict["nStepsSteepestDescentProgressInterval"]
    maxObjectiveForce=optionsDict["maxObjectiveForce"]
    timeStep=optionsDict["timeStep"]
    frictionConstant=optionsDict["frictionConstant"]
    cutOffDst=optionsDict["cutOffDst"]
    VerletListDst=optionsDict["VerletListDst"]
    inputCoordFile=optionsDict["inputCoordFile"]
    inputTopFile=optionsDict["inputTopFile"]
    lambd=optionsDict["lambd"]
    gamma=optionsDict["gamma"]
    
    writeOptionsGenericClash(path,
                             nSteps,
                             nStepsInterval,
                             nStepsWrite,
                             nStepsBackupInterval,
                             outPutFilePath,
                             outPutFormat,
                             boxX,boxY,boxZ,
                             temperature,
                             steepestDecentScale,
                             nStepsSteepestDescent,
                             nStepsSteepestDescentProgressInterval,
                             maxObjectiveForce,
                             timeStep,
                             frictionConstant,
                             cutOffDst,
                             VerletListDst,
                             inputCoordFile,
                             inputTopFile,
                             lambd,
                             gamma)


def writeOptionsGenericClashWithCompression(path:str,
                                            nSteps:int,
                                            nStepsInterval:int,
                                            nStepsWrite:int,
                                            nStepsBackupInterval:int,
                                            outPutFilePath:str,
                                            outPutFormat:str,
                                            boxX:float,boxY:float,boxZ:float,
                                            temperature:float,
                                            steepestDecentScale:float,
                                            nStepsSteepestDescent:int,
                                            nStepsSteepestDescentProgressInterval:int,
                                            maxObjectiveForce:float,
                                            timeStep:float,
                                            frictionConstant:float,
                                            cutOffDst:float,
                                            VerletListDst:float,
                                            inputCoordFile:str,
                                            inputTopFile:str,
                                            lambd:float,
                                            gamma:float,
                                            initialSphereRadius:float,
                                            minimalSphereRadius:float,
                                            compressionVelocity:float):
    
    opt = optGenericClashWithCompression_template.format(nSteps=nSteps,
                                                         nStepsInterval=nStepsInterval,
                                                         nStepsWrite=nStepsWrite,
                                                         nStepsBackupInterval=nStepsBackupInterval,
                                                         outPutFilePath=outPutFilePath,
                                                         outPutFormat=outPutFormat,
                                                         boxX=boxX,boxY=boxY,boxZ=boxZ,
                                                         temperature=temperature,
                                                         steepestDecentScale=steepestDecentScale,
                                                         nStepsSteepestDescent=nStepsSteepestDescent,
                                                         nStepsSteepestDescentProgressInterval=nStepsSteepestDescentProgressInterval,
                                                         maxObjectiveForce=maxObjectiveForce,
                                                         timeStep=timeStep,
                                                         frictionConstant=frictionConstant,
                                                         cutOffDst=cutOffDst,
                                                         VerletListDst=VerletListDst,
                                                         inputCoordFile=inputCoordFile,
                                                         inputTopFile=inputTopFile,
                                                         lambd=lambd,
                                                         gamma=gamma,
                                                         initialSphereRadius=initialSphereRadius,
                                                         minimalSphereRadius=minimalSphereRadius,
                                                         compressionVelocity=compressionVelocity)
    
    with open(path,"w") as f:
        f.write(opt)

def writeOptionsGenericClashWithCompressionFromDict(path,optionsDict):
    
    nSteps=optionsDict["nSteps"]
    nStepsInterval=optionsDict["nStepsInterval"]
    nStepsWrite=optionsDict["nStepsWrite"]
    nStepsBackupInterval=optionsDict["nStepsBackupInterval"]
    outPutFilePath=optionsDict["outPutFilePath"]
    outPutFormat=optionsDict["outPutFormat"]
    boxX=optionsDict["boxX"]
    boxY=optionsDict["boxY"]
    boxZ=optionsDict["boxZ"]
    temperature=optionsDict["temperature"]
    steepestDecentScale=optionsDict["steepestDecentScale"]
    nStepsSteepestDescent=optionsDict["nStepsSteepestDescent"]
    nStepsSteepestDescentProgressInterval=optionsDict["nStepsSteepestDescentProgressInterval"]
    maxObjectiveForce=optionsDict["maxObjectiveForce"]
    timeStep=optionsDict["timeStep"]
    frictionConstant=optionsDict["frictionConstant"]
    cutOffDst=optionsDict["cutOffDst"]
    VerletListDst=optionsDict["VerletListDst"]
    inputCoordFile=optionsDict["inputCoordFile"]
    inputTopFile=optionsDict["inputTopFile"]
    lambd=optionsDict["lambd"]
    gamma=optionsDict["gamma"]
    initialSphereRadius=optionsDict["initialSphereRadius"]
    minimalSphereRadius=optionsDict["minimalSphereRadius"]
    compressionVelocity=optionsDict["compressionVelocity"]
    
    writeOptionsGenericClashWithCompression(path,
                                            nSteps,
                                            nStepsInterval,
                                            nStepsWrite,
                                            nStepsBackupInterval,
                                            outPutFilePath,
                                            outPutFormat,
                                            boxX,boxY,boxZ,
                                            temperature,
                                            steepestDecentScale,
                                            nStepsSteepestDescent,
                                            nStepsSteepestDescentProgressInterval,
                                            maxObjectiveForce,
                                            timeStep,
                                            frictionConstant,
                                            cutOffDst,
                                            VerletListDst,
                                            inputCoordFile,
                                            inputTopFile,
                                            lambd,
                                            gamma,
                                            initialSphereRadius,
                                            minimalSphereRadius,
                                            compressionVelocity)


def writeOptionsAFM(path:str,
                    nSteps:int,
                    nStepsInterval:int,
                    nStepsWrite:int,
                    nStepsBackupInterval:int,
                    outPutFilePath:str,
                    outPutFormat:str,
                    boxX:float,boxY:float,boxZ:float,
                    temperature:float,
                    steepestDecentScale:float,
                    nStepsSteepestDescent:int,
                    nStepsSteepestDescentProgressInterval:int,
                    maxObjectiveForce:float,
                    timeStep:float,
                    frictionConstant:float,
                    frictionConstantTip:float,
                    cutOffDst:float,
                    VerletListDst:float,
                    dielectricConstant:float,
                    debyeLength:float,
                    inputCoordFile:str,
                    inputTopFile:str,
                    epsilonSurf:float,
                    sigmaSurf:float,
                    surfacePosition:float,
                    initialTipSampleDst:float,
                    descentVelocity:float,
                    minimalChipHeight:float,
                    Mtip:float,
                    Rtip:float,
                    Kxytip:float,
                    Ktip:float,
                    epsilonTip:float,
                    sigmaTip:float,
                    Atip:float,
                    Btip:float,
                    epsilonTipSurf:float,
                    sigmaTipSurf  :float,
                    ATipSurf:float,
                    BTipSurf:float,
                    nStepsIndentMeasure:str,
                    outputIndentationMeasureFilePath:str):

    opt = optAFM_template.format(nSteps=nSteps,
                                 nStepsInterval=nStepsInterval,
                                 nStepsWrite=nStepsWrite,
                                 nStepsBackupInterval=nStepsBackupInterval,
                                 outPutFilePath=outPutFilePath,
                                 outPutFormat=outPutFormat,
                                 boxX=boxX,boxY=boxY,boxZ=boxZ,
                                 temperature=temperature,
                                 steepestDecentScale=steepestDecentScale,
                                 nStepsSteepestDescent=nStepsSteepestDescent,
                                 nStepsSteepestDescentProgressInterval=nStepsSteepestDescentProgressInterval,
                                 maxObjectiveForce=maxObjectiveForce,
                                 timeStep=timeStep,
                                 frictionConstant=frictionConstant,
                                 frictionConstantTip=frictionConstantTip,
                                 cutOffDst=cutOffDst,
                                 VerletListDst=VerletListDst,
                                 dielectricConstant=dielectricConstant,
                                 debyeLength=debyeLength,
                                 inputCoordFile=inputCoordFile,
                                 inputTopFile=inputTopFile,
                                 epsilonSurf=epsilonSurf,
                                 sigmaSurf=sigmaSurf,
                                 surfacePosition=surfacePosition,
                                 initialTipSampleDst=initialTipSampleDst,
                                 descentVelocity=descentVelocity,
                                 minimalChipHeight=minimalChipHeight,
                                 Mtip=Mtip,
                                 Rtip=Rtip,
                                 Kxytip=Kxytip,
                                 Ktip=Ktip,
                                 epsilonTip=epsilonTip,
                                 sigmaTip=sigmaTip,
                                 Atip=Atip,
                                 Btip=Btip,
                                 epsilonTipSurf=epsilonTipSurf,
                                 sigmaTipSurf  =sigmaTipSurf,
                                 ATipSurf=ATipSurf,
                                 BTipSurf=BTipSurf,
                                 nStepsIndentMeasure=nStepsIndentMeasure,
                                 outputIndentationMeasureFilePath=outputIndentationMeasureFilePath)
    
    with open(path,"w") as f:
        f.write(opt)

def writeOptionsAFMFromDict(path,optionsDict):
    
    nSteps=optionsDict["nSteps"]
    nStepsInterval=optionsDict["nStepsInterval"]
    nStepsWrite=optionsDict["nStepsWrite"]
    nStepsBackupInterval=optionsDict["nStepsBackupInterval"]
    outPutFilePath=optionsDict["outPutFilePath"]
    outPutFormat=optionsDict["outPutFormat"]
    boxX=optionsDict["boxX"]
    boxY=optionsDict["boxY"]
    boxZ=optionsDict["boxZ"]
    temperature=optionsDict["temperature"]
    steepestDecentScale=optionsDict["steepestDecentScale"]
    nStepsSteepestDescent=optionsDict["nStepsSteepestDescent"]
    nStepsSteepestDescentProgressInterval=optionsDict["nStepsSteepestDescentProgressInterval"]
    maxObjectiveForce=optionsDict["maxObjectiveForce"]
    timeStep=optionsDict["timeStep"]
    frictionConstant=optionsDict["frictionConstant"]
    frictionConstantTip=optionsDict["frictionConstantTip"]
    cutOffDst=optionsDict["cutOffDst"]
    VerletListDst=optionsDict["VerletListDst"]
    dielectricConstant=optionsDict["dielectricConstant"]
    debyeLength=optionsDict["debyeLength"]
    inputCoordFile=optionsDict["inputCoordFile"]
    inputTopFile=optionsDict["inputTopFile"]
    epsilonSurf=optionsDict["epsilonSurf"]
    sigmaSurf=optionsDict["sigmaSurf"]
    surfacePosition=optionsDict["surfacePosition"]
    initialTipSampleDst=optionsDict["initialTipSampleDst"]
    descentVelocity=optionsDict["descentVelocity"]
    minimalChipHeight=optionsDict["minimalChipHeight"]
    Mtip=optionsDict["Mtip"]
    Rtip=optionsDict["Rtip"]
    Kxytip=optionsDict["Kxytip"]
    Ktip=optionsDict["Ktip"]
    epsilonTip=optionsDict["epsilonTip"]
    sigmaTip=optionsDict["sigmaTip"]
    Atip=optionsDict["Atip"]
    Btip=optionsDict["Btip"]
    epsilonTipSurf=optionsDict["epsilonTipSurf"]
    sigmaTipSurf  =optionsDict["sigmaTipSurf"]
    ATipSurf=optionsDict["ATipSurf"]
    BTipSurf=optionsDict["BTipSurf"]
    nStepsIndentMeasure=optionsDict["nStepsIndentMeasure"]
    outputIndentationMeasureFilePath=optionsDict["outputIndentationMeasureFilePath"]

    writeOptionsAFM(path,
                    nSteps,
                    nStepsInterval,
                    nStepsWrite,
                    nStepsBackupInterval,
                    outPutFilePath,
                    outPutFormat,
                    boxX,boxY,boxZ,
                    temperature,
                    steepestDecentScale,
                    nStepsSteepestDescent,
                    nStepsSteepestDescentProgressInterval,
                    maxObjectiveForce,
                    timeStep,
                    frictionConstant,
                    frictionConstantTip,
                    cutOffDst,
                    VerletListDst,
                    dielectricConstant,
                    debyeLength,
                    inputCoordFile,
                    inputTopFile,
                    epsilonSurf,
                    sigmaSurf,
                    surfacePosition,
                    initialTipSampleDst,
                    descentVelocity,
                    minimalChipHeight,
                    Mtip,
                    Rtip,
                    Kxytip,
                    Ktip,
                    epsilonTip,
                    sigmaTip,
                    Atip,
                    Btip,
                    epsilonTipSurf,
                    sigmaTipSurf,
                    ATipSurf,
                    BTipSurf,
                    nStepsIndentMeasure,
                    outputIndentationMeasureFilePath)

def writeOptionsGenericUmbrella(path:str,
                                nSteps:int,
                                nStepsInterval:int,
                                nStepsWrite:int,
                                nStepsBackupInterval:int,
                                outPutFilePath:str,
                                outPutFormat:str,
                                boxX:float,boxY:float,boxZ:float,
                                temperature:float,
                                steepestDecentScale:float,
                                nStepsSteepestDescent:int,
                                nStepsSteepestDescentProgressInterval:int,
                                maxObjectiveForce:float,
                                timeStep:float,
                                frictionConstant:float,
                                cutOffDst:float,
                                VerletListDst:float,
                                dielectricConstant:float,
                                debyeLength:float,
                                inputCoordFile:str,
                                inputTopFile:str,
                                umbrellaK:float,
                                umbrellaInit:float,
                                umbrellaEnd:float,
                                umbrellaWindowsNumber:int, 
                                umbrellaCopies:int,
                                nStepsUmbrellaMeasure:int,
                                outputUmbrellaMeasureFilePath:str):

    opt = optGenericUmbrella_template.format(nSteps=nSteps,
                                             nStepsInterval=nStepsInterval,
                                             nStepsWrite=nStepsWrite,
                                             nStepsBackupInterval=nStepsBackupInterval,
                                             outPutFilePath=outPutFilePath,
                                             outPutFormat=outPutFormat,
                                             boxX=boxX,boxY=boxY,boxZ=boxZ,
                                             temperature=temperature,
                                             steepestDecentScale=steepestDecentScale,
                                             nStepsSteepestDescent=nStepsSteepestDescent,
                                             nStepsSteepestDescentProgressInterval=nStepsSteepestDescentProgressInterval,
                                             maxObjectiveForce=maxObjectiveForce,
                                             timeStep=timeStep,
                                             frictionConstant=frictionConstant,
                                             cutOffDst=cutOffDst,
                                             VerletListDst=VerletListDst,
                                             dielectricConstant=dielectricConstant,
                                             debyeLength=debyeLength,
                                             inputCoordFile=inputCoordFile,
                                             inputTopFile=inputTopFile,
                                             umbrellaK=umbrellaK,
                                             umbrellaInit=umbrellaInit,
                                             umbrellaEnd=umbrellaEnd,
                                             umbrellaWindowsNumber=umbrellaWindowsNumber, 
                                             umbrellaCopies=umbrellaCopies,
                                             nStepsUmbrellaMeasure=nStepsUmbrellaMeasure,
                                             outputUmbrellaMeasureFilePath=outputUmbrellaMeasureFilePath)
    
    with open(path,"w") as f:
        f.write(opt)

def writeOptionsGenericUmbrellaFromDict(path,optionsDict):
    
    nSteps=optionsDict["nSteps"]
    nStepsInterval=optionsDict["nStepsInterval"]
    nStepsWrite=optionsDict["nStepsWrite"]
    nStepsBackupInterval=optionsDict["nStepsBackupInterval"]
    outPutFilePath=optionsDict["outPutFilePath"]
    outPutFormat=optionsDict["outPutFormat"]
    boxX=optionsDict["boxX"]
    boxY=optionsDict["boxY"]
    boxZ=optionsDict["boxZ"]
    temperature=optionsDict["temperature"]
    steepestDecentScale=optionsDict["steepestDecentScale"]
    nStepsSteepestDescent=optionsDict["nStepsSteepestDescent"]
    nStepsSteepestDescentProgressInterval=optionsDict["nStepsSteepestDescentProgressInterval"]
    maxObjectiveForce=optionsDict["maxObjectiveForce"]
    timeStep=optionsDict["timeStep"]
    frictionConstant=optionsDict["frictionConstant"]
    cutOffDst=optionsDict["cutOffDst"]
    VerletListDst=optionsDict["VerletListDst"]
    dielectricConstant=optionsDict["dielectricConstant"]
    debyeLength=optionsDict["debyeLength"]
    inputCoordFile=optionsDict["inputCoordFile"]
    inputTopFile=optionsDict["inputTopFile"]
    umbrellaK=optionsDict["umbrellaK"]
    umbrellaInit=optionsDict["umbrellaInit"]
    umbrellaEnd=optionsDict["umbrellaEnd"]
    umbrellaWindowsNumber=optionsDict["umbrellaWindowsNumber"]
    umbrellaCopies=optionsDict["umbrellaCopies"]
    nStepsUmbrellaMeasure=optionsDict["nStepsUmbrellaMeasure"]
    outputUmbrellaMeasureFilePath=optionsDict["outputUmbrellaMeasureFilePath"]

    writeOptionsGenericUmbrella(path,
                                nSteps,
                                nStepsInterval,
                                nStepsWrite,
                                nStepsBackupInterval,
                                outPutFilePath,
                                outPutFormat,
                                boxX,boxY,boxZ,
                                temperature,
                                steepestDecentScale,
                                nStepsSteepestDescent,
                                nStepsSteepestDescentProgressInterval,
                                maxObjectiveForce,
                                timeStep,
                                frictionConstant,
                                cutOffDst,
                                VerletListDst,
                                dielectricConstant,
                                debyeLength,
                                inputCoordFile,
                                inputTopFile,
                                umbrellaK,
                                umbrellaInit,
                                umbrellaEnd,
                                umbrellaWindowsNumber,
                                umbrellaCopies,
                                nStepsUmbrellaMeasure,
                                outputUmbrellaMeasureFilePath)

################################################


#writeOptionsGeneric("optionsTestGeneric.dat",
#                    1000,
#                    1000,
#                    1000,
#                    1000,
#                    "kk",
#                    "asdf",
#                    1.23,2.13,3.12,
#                    1.0,
#                    1.5,
#                    55555,
#                    4444,
#                    -1.0,
#                    0.123456,
#                    2.0,
#                    100.0,
#                    150.0,
#                    "asdf.coord",
#                    "asdf.top")
#
#writeOptionsGenericWithSurface("optionsTestGenericSurface.dat",
#                               1000,
#                               1000,
#                               1000,
#                               1000,
#                               "kk",
#                               "asdf",
#                               1.23,2.13,3.12,
#                               1.0,
#                               1.5,
#                               55555,
#                               4444,
#                               -1.0,
#                               0.123456,
#                               2.0,
#                               100.0,
#                               150.0,
#                               "asdf.coord",
#                               "asdf.top",
#                               1.0,
#                               2.0,
#                               -300)
#
#writeOptionsGenericClash("optionsTestClash.dat",
#                          1000,
#                          1000,
#                          1000,
#                          1000,
#                          "kk",
#                          "asdf",
#                          1.23,2.13,3.12,
#                          1.0,
#                          1.5,
#                          55555,
#                          4444,
#                          -1.0,
#                          0.123456,
#                          2.0,
#                          100.0,
#                          150.0,
#                          "asdf.coord",
#                          "asdf.top",
#                          2.0,
#                          4.0)
#        
#writeOptionsGenericClashWithCompression("optionsTestClashWithCompression.dat",
#                                         1000,
#                                         1000,
#                                         1000,
#                                         1000,
#                                         "kk",
#                                         "asdf",
#                                         1.23,2.13,3.12,
#                                         1.0,
#                                         1.5,
#                                         55555,
#                                         4444,
#                                         -1.0,
#                                         0.123456,
#                                         2.0,
#                                         100.0,
#                                         150.0,
#                                         "asdf.coord",
#                                         "asdf.top",
#                                         2.0,
#                                         4.0,
#                                         300,
#                                         400,
#                                         0.01001)
#
#writeOptionsAFM("optionsAFM.dat",
#                 1000,
#                 1000,
#                 1000,
#                 1000,
#                 "kk",
#                 "asdf",
#                 1.23,2.13,3.12,
#                 1.0,
#                 1.5,
#                 55555,
#                 4444,
#                 -1.0,
#                 0.123456,
#                 2.0,
#                 100.0,
#                 150.0,
#                 "asdf.coord",
#                 "asdf.top",
#                 1.0,
#                 2.0,
#                 -300,
#                1.123,
#                0.001,
#                2.323,
#                1.02,
#                5.02,
#                5698,
#                456,
#                5897,
#                57,
#                142566,
#                146,
#                1486,
#                98765,
#                12453,
#                1289)


