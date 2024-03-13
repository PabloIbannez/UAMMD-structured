
def get_r0(cls,subCls,bondsList):
    if   cls == "Bond1":
        r0 = 0.0
    elif cls == "Bond2":
        if "r0" in bondsList[subCls]["commonParams"]:
            r0 = bondsList[subCls]["commonParams"]["r0"]
        elif "r0" in bondsList[subCls]["bondParams"]:
            r0 = bondsList[subCls]["bondParams"]["r0"]
        elif "sigma" in bondsList[subCls]["commonParams"]:
            r0 = bondsList[subCls]["commonParams"]["sigma"]*1.5
        elif "sigma" in bondsList[subCls]["bondParams"]:
            r0 = bondsList[subCls]["bondParams"]["sigma"]*1.5
        elif "debyeLength" in bondsList[subCls]["commonParams"]:
            r0 = bondsList[subCls]["commonParams"]["debyeLength"]*1.5
        elif "debyeLength" in bondsList[subCls]["bondParams"]:
            r0 = bondsList[subCls]["bondParams"]["debyeLength"]*1.5
        elif "maxDistance" in bondsList[subCls]["commonParams"]:
            r0 = bondsList[subCls]["commonParams"]["maxDistance"]
        elif "maxDistance" in bondsList[subCls]["bondParams"]:
            r0 = bondsList[subCls]["bondParams"]["maxDistance"]
        else:
            print("Error: r0 or sigma not found for bond type %s" % subCls)
            exit(1)
    elif cls == "Bond3":
        r0 = bondsList["HarmonicCommon_K_r0"]["commonParams"]["r0"]
    elif cls == "Bond4":
        r0 = bondsList["HarmonicCommon_K_r0"]["commonParams"]["r0"]
    else:
        print("Error: Unknown bond class %s" % cls)
        exit(1)

    return r0
