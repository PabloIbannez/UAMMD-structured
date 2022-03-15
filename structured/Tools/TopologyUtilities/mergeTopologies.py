import sys

from Topology import *

top1 = Topology(sys.argv[1],sys.argv[2])
top2 = Topology(sys.argv[3],sys.argv[4])

if(len(sys.argv)>6):
    merging2File(top1,top2,sys.argv[5],sys.argv[6])
else:
    merging2File(top1,top2,sys.argv[5])
