#!/usr/bin/env python
import sys, string,os;
from util import *;

try:
    (inputfile,jobdir,start,end,N) = sys.argv[1:];
    (start,end,N) = (int(start),int(end),int(N));
except:
    print >> sys.stderr, "Syntax: %s <input-file> <job-directory> <start> <end> <number-of-chunks>\n" % sys.argv[0];
    exit(-1);


length    = end-start+1;
chunksize = length/N;
chunkrem  = length % N;

print "(length,chunksize,remainder) = ",(length,chunksize,chunkrem);

text   = readfile(inputfile);
lltext = readfile("llsubmit.sh");

def prepare(dirname,s,e):
    mkdirp(dirname);

    with open(dirname+"/info.txt","w") as f:
        print >> f, "Parameters: ", (inputfile,jobdir,start,end,N);
        print >> f, "Index: ",i;
        print >> f, "Interval: "+ str(s)+"-"+str(e);

    with open(dirname+"/input.inp","w") as f:
        f.write(replace_input(text,{"FROM":s,"TO":e}));
    
for i in range(N):
    s = i*chunksize+start;
    e = (i+1)*chunksize+start-1;

    prepare(jobdir+"/"+str(i),s,e);
prepare(jobdir+"/"+str(N),e+1,end);

iN = N/8;
ifrom = [i*8 for i in range(iN)]
ito   = [(i+1)*8-1 for i in range(iN)];
ito[-1] = N;

for i in range(iN):
    with open(jobdir+"/llsubmit-%d-%d.sh" % (ifrom[i],ito[i]),"w") as f:
        f.write(replace_input(lltext,{"ROOT":os.environ["PWD"]+"/..","JOBDIR":os.environ["PWD"]+"/"+jobdir,
                                      "iFROM":ifrom[i],"iTO":ito[i],"I":i}));


