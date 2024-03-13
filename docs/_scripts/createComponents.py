import json

componentsFile = "/home/pablo/UAMMD/USCM/./Components/Components.json"

forceFieldComponents = ["Interactor","DataStructures"]

with open(componentsFile) as json_file:
    components = json.load(json_file)

compSorted = {}
for tpy,tpyList in components.items():
    if tpy not in forceFieldComponents:
        compSorted.setdefault(tpy,{})
        for subTpy,subTpyList in tpyList.items():
            for cls,subCls,_ in subTpyList:
                if subTpy != cls:
                    compSorted[tpy].setdefault(subTpy,{})
                    compSorted[tpy][subTpy].setdefault(cls,[])
                    compSorted[tpy][subTpy][cls].append(subCls)
                else:
                    compSorted[tpy].setdefault(cls,[])
                    compSorted[tpy][cls].append(subCls)
    else:
        compSorted.setdefault("Topology",{})
        compSorted["Topology"].setdefault(tpy,{})

        for subTpy,subTpyList in tpyList.items():
            for cls,subCls,_ in subTpyList:
                if subTpy != cls:
                    compSorted["Topology"][tpy].setdefault(subTpy,{})
                    compSorted["Topology"][tpy][subTpy].setdefault(cls,[])
                    compSorted["Topology"][tpy][subTpy][cls].append(subCls)
                else:
                    compSorted["Topology"][tpy].setdefault(cls,[])
                    compSorted["Topology"][tpy][cls].append(subCls)

#Next part commented for safety

#for tpy,tpyList in compSorted.items():
#    with open(tpy+".rst","w") as fout:
#        fout.write("#"*len(tpy)+"\n")
#        fout.write(tpy+"\n")
#        fout.write("#"*len(tpy)+"\n\n")
#        for subTpy,subTpyList in tpyList.items():
#            print("\t"+subTpy)
#            fout.write("*"*len(subTpy)+"\n")
#            fout.write(subTpy+"\n")
#            fout.write("*"*len(subTpy)+"\n\n")
#            if type(subTpyList) is dict:
#                for cls,clsList in subTpyList.items():
#                    print("\t\t"+cls)
#                    fout.write(cls+"\n")
#                    fout.write("="*len(cls)+"\n\n")
#                    for subCls in clsList:
#                        print("\t\t\t"+subCls)
#                        fout.write(subCls+"\n")
#                        fout.write("-"*len(subCls)+"\n\n")
#                        try:
#                            for subSubCls in clsList[subCls]:
#                                print("\t\t\t\t"+subSubCls)
#                                fout.write(subSubCls+"\n")
#                                fout.write("^"*len(subSubCls)+"\n\n")
#                        except:
#                            pass
#            elif type(subTpyList) is list:
#                for cls in subTpyList:
#                    print("\t\t"+cls)
#                    fout.write(cls+"\n")
#                    fout.write("="*len(cls)+"\n\n")




