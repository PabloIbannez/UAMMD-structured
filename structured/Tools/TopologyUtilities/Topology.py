import sys
import copy

import re

from scanf import scanf
from collections import OrderedDict


class Topology:

    def __isEmptyOrCommented(self,line:str):

        if(line.isspace()):
            return True

        first = line[0]
        if(first=="#"):
            return True
        if(first==";"):
            return True
        return False
    
    def __nextItemLine(self,item:str):
        if((list(self.d.keys()).index(item)+1)==(len(list(self.d)))):
            return self.nLines
        return self.d[list(self.d)[list(self.d.keys()).index(item) + 1]]

    def __loadLabels(self):
        self.d=OrderedDict()
        for i,line in enumerate(self.TopologyFile):
            a, b = line.find('['), line.find(']')
            if a != -1 and b != -1:
                label = line[a+1:b]
                label = label.strip()
                self.d[label]=i
        self.TopologyFile.seek(0)
    
    def __getCoord(self):
       
        c = []

        self.CoordFile.seek(0)
        
        bre = re.compile(r'^(\d+)\s+(.+)$')
        
        for i,line in enumerate(self.CoordFile):
            cb = []
            m = bre.match(line)
            cb.append(int(m.group(1)))
            cb.append(m.group(2))
            c.append(cb)

        return c
    
    def __getCommon(self,label:str):
       
        p = []

        self.TopologyFile.seek(0)
        
        start = self.d[label]
        end   = self.__nextItemLine(label)

        #print(start,end)

        for i,line in enumerate(self.TopologyFile):
            if(i>start and i<end and self.__isEmptyOrCommented(line)==False):
                pb = []
                pb.append(line)
                p.append(pb)

        return p
    
    def __getStructure(self,label:str):
       
        p = []

        self.TopologyFile.seek(0)
        
        start = self.d[label]
        end   = self.__nextItemLine(label)

        #print(start,end)

        bre = re.compile(r'^(\d+)\s+(\w+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d*)\s*$')

        for i,line in enumerate(self.TopologyFile):
            if(i>start and i<end and self.__isEmptyOrCommented(line)==False):
                m = bre.match(line)
                #print(line,"  :  ",int(m.group(1)),m.group(2),int(m.group(3)),int(m.group(4)),int(m.group(5)),int(m.group(6)))
                pb = [int(m.group(1)),
                          m.group(2) ,
                      int(m.group(3)),
                      int(m.group(4)),
                      int(m.group(5))]

                if(not m.group(6)):
                    pb.append(0)
                else:
                    pb.append(int(m.group(6)))

                p.append(pb)

        return p
    
    def __getNeigList(self,label:str):

        p = []

        self.TopologyFile.seek(0)
        
        start = self.d[label]
        end   = self.__nextItemLine(label)
        
        #print(start,end)

        for i,line in enumerate(self.TopologyFile):
            if(i>start and i<end and self.__isEmptyOrCommented(line)==False):
                pb = []
                for e in line.split():
                    pb.append(int(e))
                p.append(pb)

        return p

    
    def __getProperty(self,label:str,n:int):
       
        p = []

        self.TopologyFile.seek(0)
        
        start = self.d[label]
        end   = self.__nextItemLine(label)

        #print(start,end)

        rgx = r'^'
        for i in range(n):
            rgx=rgx+r'\s*(\d+)\s+'
        rgx=rgx+r'(.*)$'

        bre = re.compile(rgx)

        for i,line in enumerate(self.TopologyFile):
            if(i>start and i<end and self.__isEmptyOrCommented(line)==False):
                m = bre.match(line)
                #print(line,m)
                pb = []
                for i in range(n):
                    pb.append(int(m.group(i+1)))
                pb.append(m.group(n+1))
                p.append(pb)

        return p


    def __init__(self,
                 CoordFilePath    : str,
                 TopologyFilePath : str):

        self.propertiesIgnored = ["defaults",
                                  "atomtypes",
                                  "moleculetype",
                                  "system",
                                  "molecules"]
        
        self.propertiesCommon = ["TYPES",
                                 "MJ",
                                 "PDNA",
                                 "LJ_WCA",
                                 "UNBOUND",
                                 "CLASH",
                                 "POLAR",
                                 "STERIC",
                                 "STERIC_INTRA",
                                 "STERIC_INTER",
                                 "NONPOLAR"]

        self.propertyTypesDict   = {"Bond1":1,
                                    "Bond2":2,
                                    "Bond3":3,
                                    "Bond4":4}

        self.propertiesDict   = {"SASA":1,
                                 "BONDS":2,
                                 "BONDS_PROT":2,
                                 "BONDS_WLC":2,
                                 "BONDS_DNA":2,
                                 "BONDS_RNA":2,
                                 "ANGLES":3,
                                 "ANGLES_WLC":3,
                                 "ANGLES_DNA":3,
                                 "ANGLES_RNA":3,
                                 "DIHEDRALS":4,
                                 "FIXED":1,
                                 "BONDS_DH":2,
                                 "SOP_COVALENT":2,
                                 "PAIRS":2,
                                 "PAIRS_PROT":2,
                                 "SOP_NATIVE_CONTACT":2,
                                 "ENM_BONDS":2,
                                 "EXCLUSIONS":-1,
                                 "atoms":1,
                                 "bonds":2,
                                 "angles":3,
                                 "dihedrals":4,
                                 "pairs":2,
                                 "exclusions":1}

    
        self.propertiesLoaded = {}
        
        self.CoordFilePath = CoordFilePath
        self.CoordFile     = open(self.CoordFilePath,'r')

        self.coordLoaded = self.__getCoord()

        self.TopologyFilePath = TopologyFilePath
        self.TopologyFile     = open(self.TopologyFilePath,'r')

        self.nLines = sum(1 for line in self.TopologyFile)
        self.TopologyFile.seek(0)
        
        self.__loadLabels()
        
        #print(self.d)

        for i in self.d:
            
            tpy=""
            if len(i.split(";"))>1:
                tpy=i.split(";")[1]

            if   (i in self.propertiesCommon):
                #print("Detected common property:", i)
                self.propertiesLoaded[i]=self.__getCommon(i)
            elif (i == "STRUCTURE"):
                #print("Detected structure:", i)
                self.propertiesLoaded[i]=self.__getStructure(i)
            elif (tpy in self.propertyTypesDict.keys()):
                self.propertiesDict[i]=self.propertyTypesDict[tpy]
            elif (tpy == "UnBound"):
                self.propertiesCommon.append(i)
                self.propertiesLoaded[i]=self.__getCommon(i)
            elif (i in self.propertiesDict):
                #print("Detected property:", i)
                n = self.propertiesDict[i]
                if(n<0):
                    self.propertiesLoaded[i]=self.__getNeigList(i)
                else:
                    self.propertiesLoaded[i]=self.__getProperty(i,n)
            elif (i in self.propertiesIgnored):
                continue
            else:
                print("Topology property", i, "not listed anywhere")
                sys.exit(0)

        #for i in self.propertiesLoaded["EXCLUSIONS"]: 
        #    print(i)
        
        #for i in self.coordLoaded: 
        #    print(i)
    
    def addProperty(self,name:str,info):
        
        if (name not in self.propertiesCommon) and (name not in self.propertiesDict.keys()): 
            tpy=""
            if len(name.split(";"))>1:
                tpy=name.split(";")[1]

            if tpy in self.propertyTypesDict.keys():
                self.propertiesDict[name]=self.propertyTypesDict[name.split(";")[0]]
            elif tpy == "UnBound":
                self.propertiesCommon.append(name)
            else:
                print("ERROR property",name," is not in any property list")
                sys.exit(1)
        
        if info==None:
            self.propertiesLoaded[name] = []
        else:
            self.propertiesLoaded[name] = info

    def renamePropertyLoaded(self,oldName:str,newName:str):

        if newName not in self.propertiesLoaded.keys():

            self.propertiesLoaded[newName] = self.propertiesLoaded.pop(oldName)
            
            if oldName in self.propertiesCommon:
                self.propertiesCommon.append(newName)
            else:
                self.propertiesDict[newName]   = self.propertiesDict[oldName]

        else:
            for i in self.propertiesLoaded[oldName]:
                self.propertiesLoaded[newName].append(copy.deepcopy(i))
            
            del self.propertiesLoaded[oldName]


    def renamePropertiesLoadedUsingDict(self,old2newNameDict:dict):
        for oldName in old2newNameDict.keys():
            if oldName in self.propertiesLoaded.keys():
                self.renamePropertyLoaded(oldName,old2newNameDict[oldName])

    def append(self,top,mode:str):
            
        offset = 0
        for c in self.coordLoaded:
            pass
        offset = c[0]+1

        for c in top.coordLoaded:
            self.coordLoaded.append([c[0]+offset,c[1]])

        for common in top.propertiesCommon:
            if common in top.propertiesLoaded.keys():
                if common in self.propertiesLoaded.keys(): 
                    for p in  top.propertiesLoaded[common]:
                        if p not in self.propertiesLoaded[common]:
                            self.propertiesLoaded[common].append(copy.deepcopy(p))
                else:
                    self.propertiesLoaded[common] = {}
                    for p in  top.propertiesLoaded[common]:
                        self.propertiesLoaded[common].append(copy.deepcopy(p))

        if  (mode=="simId"):
            maxSimId=0
            for c in self.propertiesLoaded["STRUCTURE"]:
                if(c[5]):
                    maxSimId=max([maxSimId,c[5]])
            #print("MaxSimId:",maxSimId)

            for c in top.propertiesLoaded["STRUCTURE"]:
                tmp = [c[0]+offset,str(c[1]),c[2],c[3],c[4],c[5]+maxSimId]
                self.propertiesLoaded["STRUCTURE"].append(copy.deepcopy(tmp))

        elif(mode=="mdl"):
            maxMdl=0
            for c in self.propertiesLoaded["STRUCTURE"]:
                maxMdl=max([maxMdl,c[4]])
            if maxMdl==0:
                maxMdl=1

            for c in top.propertiesLoaded["STRUCTURE"]:
                tmp = [c[0]+offset,str(c[1]),c[2],c[3],c[4]+maxMdl,c[5]]
                self.propertiesLoaded["STRUCTURE"].append(copy.deepcopy(tmp))
        elif(mode=="none"):
            #Write structure
            for c in top.propertiesLoaded["STRUCTURE"]:
                tmp = [c[0]+offset,str(c[1]),c[2],c[3],c[4],c[5]]
                self.propertiesLoaded["STRUCTURE"].append(copy.deepcopy(tmp))
        else:
            print("ERROR in Topology.append(top,mode), not valid mode:", mode)
            sys.exit()
        
        #Load properties
        for prop in top.propertiesLoaded.keys():
            if prop in self.propertiesDict.keys():
                if prop in self.propertiesLoaded.keys():
                    ni = self.propertiesDict[prop]
                    if ni < 0:
                        for p in top.propertiesLoaded[prop]:
                            tmp = []
                            for i in p:
                                tmp.append(i+offset)
                            self.propertiesLoaded[prop].append(copy.deepcopy(tmp))
                    else:
                        for p in top.propertiesLoaded[prop]:
                            tmp = []
                            for i in range(ni):
                                tmp.append(p[i]+offset)
                            tmp.append(p[ni])
                            self.propertiesLoaded[prop].append(copy.deepcopy(tmp))

                else:
                    self.addProperty(prop,None)
                    
                    ni = self.propertiesDict[prop]
                    if ni < 0:
                        for p in top.propertiesLoaded[prop]:
                            tmp = []
                            for i in p:
                                tmp.append(i+offset)
                            self.propertiesLoaded[prop].append(copy.deepcopy(tmp))
                    else:
                        for p in top.propertiesLoaded[prop]:
                            tmp = []
                            for i in range(ni):
                                tmp.append(p[i]+offset)
                            tmp.append(p[ni])
                            self.propertiesLoaded[prop].append(copy.deepcopy(tmp))
                    

    def setSimId(self,simId:int):
        for c in self.propertiesLoaded["STRUCTURE"]:
            c[5]=simId

    def write(self,outputFile):
        with open(outputFile+'.coord','w') as cf:
            for c in self.coordLoaded:
                cf.write("%i %s\n" % (c[0],c[1]))
    
        with open(outputFile+'.top','w') as tf:
            #Write commons
            for common in self.propertiesCommon:
                if(common in self.propertiesLoaded.keys()):
                    tf.write("["+common+"]\n")
                    for i in self.propertiesLoaded[common]:
                        tf.write(i[0])
        
            tf.write("[STRUCTURE]\n")
            for c in self.propertiesLoaded["STRUCTURE"]:
                if(c[5] or c[5]==0):
                    tf.write("%i %s %i %i %i %i\n" % (c[0],str(c[1]),c[2],c[3],c[4],c[5]))
                else:
                    tf.write("%i %s %i %i %i\n" % (c[0],str(c[1]),c[2],c[3],c[4]))
            
            #Write properties
            for prop in self.propertiesDict.keys():
                if(prop in self.propertiesLoaded.keys()):
                    tf.write("["+prop+"]\n")
                    
                    ni = self.propertiesDict[prop]
                    if(ni<0):
                        for i in self.propertiesLoaded[prop]:
                            for j in i:
                                tf.write(str(j)+" ")
                            tf.write("\n")
                    else:
                        for i in self.propertiesLoaded[prop]:
                            for j in range(ni):
                                tf.write(str(i[j])+" ")
                            tf.write(i[ni]+"\n")
        return

def nCopies2File(top:Topology,nCopies:int,output:str,mode:str = "simId"):

    offsetDict = {}
    with open(output+'.coord','w') as cf:
        offset = 0
        for copy in range(nCopies):
            offsetDict[copy]=offset
            for c in top.coordLoaded:
                cf.write("%i %s\n" % (c[0]+offset,c[1]))
            offset = offset+c[0]+1
    
    with open(output+'.top','w') as tf:
        #Write commons
        for common in top.propertiesCommon:
            if(common in top.propertiesLoaded.keys()):
                tf.write("["+common+"]\n")
                for i in top.propertiesLoaded[common]:
                    tf.write(i[0])
        
        if  (mode=="simId"):
            maxSimId=0
            for c in top.propertiesLoaded["STRUCTURE"]:
                if(c[5]):
                    maxSimId=max([maxSimId,c[5]])
            #print("MaxSimId:",maxSimId)

            #Write structure
            tf.write("[STRUCTURE]\n")
            for copy in range(nCopies):
                offset = offsetDict[copy]
                for c in top.propertiesLoaded["STRUCTURE"]:
                    tf.write("%i %s %i %i %i %i\n" % (c[0]+offset,str(c[1]),c[2],c[3],c[4],c[5]+(maxSimId+1)*copy))
        
        elif(mode=="mdl"):
            maxMdl=0
            for c in top.propertiesLoaded["STRUCTURE"]:
                maxMdl=max([maxMdl,c[4]])
            if maxMdl == 0:
                maxMdl=1
            #print("maxMdl:",maxMdl)
            
            #Write structure
            tf.write("[STRUCTURE]\n")
            for copy in range(nCopies):
                offset = offsetDict[copy]
                for c in top.propertiesLoaded["STRUCTURE"]:
                    tf.write("%i %s %i %i %i %i\n" % (c[0]+offset,str(c[1]),c[2],c[3],c[4]+maxMdl*copy,c[5]))
        else:
            print("ERROR in nCopies2File, not valid mode:", mode)
            sys.exit()
        
        #Write properties
        for prop in top.propertiesDict.keys():
            if(prop in top.propertiesLoaded.keys()):
                tf.write("["+prop+"]\n")
                for copy in range(nCopies):
                    
                    offset = offsetDict[copy]
                    
                    ni = top.propertiesDict[prop]
                    if(ni<0):
                        for i in top.propertiesLoaded[prop]:
                            for j in i:
                                tf.write(str(j+offset)+" ")
                            tf.write("\n")
                    else:
                        for i in top.propertiesLoaded[prop]:
                            for j in range(ni):
                                tf.write(str(i[j]+offset)+" ")
                            tf.write(i[ni]+"\n")

def merging2File(top1:Topology,top2:Topology,output:str,mode:str = "simId"):
            
    offset = 0
    with open(output+'.coord','w') as cf:
        for c in top1.coordLoaded:
            cf.write("%i %s\n" % (c[0],c[1]))
        offset = c[0]+1
        for c in top2.coordLoaded:
            cf.write("%i %s\n" % (c[0]+offset,c[1]))
    
    commonAdded = set()
    for common in top1.propertiesCommon:
        if common in top1.propertiesLoaded.keys():
            if common in top2.propertiesCommon:
                if common in top2.propertiesLoaded.keys():
                    commonAdded.add(common)

    with open(output+'.top','w') as tf:
        #Write commons
        for common in commonAdded:
            commonEntries = []
            tf.write("["+common+"]\n")
            for i in top1.propertiesLoaded[common]:
                tf.write(i[0])
                commonEntries.append(i[0])
            for i in top2.propertiesLoaded[common]:
                if i[0] not in commonEntries:
                    tf.write(i[0])
        
        for common in top1.propertiesCommon:
            if(common in top1.propertiesLoaded.keys()):
                if(common not in commonAdded):
                    tf.write("["+common+"]\n")
                    for i in top1.propertiesLoaded[common]:
                        tf.write(i[0])
        
        for common in top2.propertiesCommon:
            if(common in top2.propertiesLoaded.keys()):
                if(common not in commonAdded):
                    tf.write("["+common+"]\n")
                    for i in top2.propertiesLoaded[common]:
                        tf.write(i[0])
        
        if  (mode=="simId"):
            maxSimId=0
            for c in top1.propertiesLoaded["STRUCTURE"]:
                if(c[5]):
                    maxSimId=max([maxSimId,c[5]])
            #print("MaxSimId:",maxSimId)

            #Write structure
            tf.write("[STRUCTURE]\n")
            for c in top1.propertiesLoaded["STRUCTURE"]:
                tf.write("%i %s %i %i %i %i\n" % (c[0]       ,str(c[1]),c[2],c[3],c[4],c[5]))
            for c in top2.propertiesLoaded["STRUCTURE"]:
                tf.write("%i %s %i %i %i %i\n" % (c[0]+offset,str(c[1]),c[2],c[3],c[4],c[5]+maxSimId+1))
        elif(mode=="mdl"):
            maxMdl=0
            for c in top1.propertiesLoaded["STRUCTURE"]:
                maxMdl=max([maxMdl,c[4]])
            if maxMdl==0:
                maxMdl=1

            #Write structure
            tf.write("[STRUCTURE]\n")
            for c in top1.propertiesLoaded["STRUCTURE"]:
                tf.write("%i %s %i %i %i %i\n" % (c[0]       ,str(c[1]),c[2],c[3],c[4],c[5]))
            for c in top2.propertiesLoaded["STRUCTURE"]:
                tf.write("%i %s %i %i %i %i\n" % (c[0]+offset,str(c[1]),c[2],c[3],c[4]+maxMdl,c[5]))
        elif(mode=="none"):
            #Write structure
            tf.write("[STRUCTURE]\n")
            for c in top1.propertiesLoaded["STRUCTURE"]:
                tf.write("%i %s %i %i %i %i\n" % (c[0]       ,str(c[1]),c[2],c[3],c[4],c[5]))
            for c in top2.propertiesLoaded["STRUCTURE"]:
                tf.write("%i %s %i %i %i %i\n" % (c[0]+offset,str(c[1]),c[2],c[3],c[4],c[5]))
        else:
            print("ERROR in merging2File, not valid mode:", mode)
            sys.exit()
        
        #Write properties
        for prop in top1.propertiesDict.keys():
            if  (prop in top1.propertiesLoaded.keys()):
                tf.write("["+prop+"]\n")
                ni = top1.propertiesDict[prop]
                if(ni<0):
                    for i in top1.propertiesLoaded[prop]:
                        for j in i:
                            tf.write(str(j)+" ")
                        tf.write("\n")
                    if(prop in top2.propertiesLoaded.keys()):
                        for i in top2.propertiesLoaded[prop]:
                            for j in i:
                                tf.write(str(j+offset)+" ")
                            tf.write("\n")
                else:
                    for i in top1.propertiesLoaded[prop]:
                        for j in range(ni):
                            tf.write(str(i[j])+" ")
                        tf.write(i[ni]+"\n")
                    if(prop in top2.propertiesLoaded.keys()):
                        for i in top2.propertiesLoaded[prop]:
                            for j in range(ni):
                                tf.write(str(i[j]+offset)+" ")
                            tf.write(i[ni]+"\n")
            elif(prop in top2.propertiesLoaded.keys()):
                tf.write("["+prop+"]\n")
                ni = top2.propertiesDict[prop]
                if(ni<0):
                    for i in top2.propertiesLoaded[prop]:
                        for j in i:
                            tf.write(str(j+offset)+" ")
                        tf.write("\n")
                else:
                    for i in top2.propertiesLoaded[prop]:
                        for j in range(ni):
                            tf.write(str(i[j]+offset)+" ")
                        tf.write(i[ni]+"\n")

