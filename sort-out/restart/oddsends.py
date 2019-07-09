#!/usr/bin/env python
import sys, string,os;
from util import *;

try:
    (inputfile,jobdir,notdonefile) = sys.argv[1:];
except:
    print >> sys.stderr, "Syntax: %s <input-file> <job-directory> <not-done-file>\n" % sys.argv[0];
    exit(-1);

with open(notdonefile,"r") as f:
    notdone = [[int(s) for s in l.split()] for l in f.readlines()];

text   = readfile(inputfile);
lltext = readfile("llsubmit.sh");

def prepare(dirname,s,e):
    mkdirp(dirname);

    with open(dirname+"/info.txt","w") as f:
        print >> f, "Index: ",i;
        print >> f, "Interval: "+ str(s)+"-"+str(e);

    with open(dirname+"/input.inp","w") as f:
        f.write(replace_input(text,{"FROM":s,"TO":e}));
    
for i in range(len(notdone)):
    [s,e] = notdone[i];
    prepare(jobdir+"/"+str(i),s,e);

N  = len(notdone);
iN = N/8;
ifrom = [i*8 for i in range(iN)]
ito   = [(i+1)*8-1 for i in range(iN)];
ito[-1] = N;

for i in range(iN):
    with open(jobdir+"/llsubmit-%d-%d.sh" % (ifrom[i],ito[i]),"w") as f:
        f.write(replace_input(lltext,{"ROOT":os.environ["PWD"]+"/..","JOBDIR":os.environ["PWD"]+"/"+jobdir,
                                      "iFROM":ifrom[i],"iTO":ito[i],"I":i}));


