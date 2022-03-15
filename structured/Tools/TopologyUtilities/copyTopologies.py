import sys

from Topology import *

top = Topology(sys.argv[1],sys.argv[2])

if(len(sys.argv)>5):
    nCopies2File(top,int(sys.argv[3]),sys.argv[4],sys.argv[5])
else:
    nCopies2File(top,int(sys.argv[3]),sys.argv[4])
