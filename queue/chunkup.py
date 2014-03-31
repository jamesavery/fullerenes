#!/usr/bin/env python
import sys, string,os;
from util import *;

try:
    (execname,paramstring,jobdir,start,end,N,NCPUs) = sys.argv[1:];
    (start,end,N,NCPUs) = (int(start),int(end),int(N),int(NCPUs));
except:
    print >> sys.stderr, "Syntax: %s <exec-name> [file:<input-file>|<input-string>] <job-directory> <start> <end> <number-of-chunks> <cpus-per-task>\n" % sys.argv[0];
    exit(-1);


length    = end-start+1;
chunksize = length/N;
chunkrem  = length % N;

print "(length,chunksize,remainder) = ",(length,chunksize,chunkrem);

if(paramstring[0:4] == "file:"):
    inputfile = paramstring[4:-1];
    text   = readfile(inputfile);
else:
    text = paramstring;

lltext = readfile("llsubmit.sh");
slurmtext = readfile("slurmsubmit.sh");

def prepare(dirname,s,e):
    mkdirp(dirname);

    with open(dirname+"/info.txt","w") as f:
        print >> f, "Parameters: ", (paramstring,jobdir,start,end,N);
        print >> f, "Index: ",i;
        print >> f, "Interval: "+ str(s)+"-"+str(e);

    with open(dirname+"/input.inp","w") as f:
        f.write(replace_input(text,{"FROM":s,"TO":e}));
    
for i in range(N):
    s = i*chunksize+start;
    e = (i+1)*chunksize+start-1;

    prepare(jobdir+"/"+str(i),s,e);
prepare(jobdir+"/"+str(N),e+1,end);

iN = N/NCPUs;
ifrom = [i*NCPUs for i in range(iN)]
ito   = [(i+1)*NCPUs-1 for i in range(iN)];
ito[-1] = N;

for i in range(iN):
    with open(jobdir+"/llsubmit-%d-%d.sh" % (ifrom[i],ito[i]),"w") as f:
        f.write(replace_input(lltext,{"EXECNAME":execname,
					"ROOT":os.environ["PWD"]+"/..","JOBDIR":os.environ["PWD"]+"/"+jobdir,
                                      "iFROM":ifrom[i],"iTO":ito[i],"I":i,"NCPUS":NCPUs}));

    with open(jobdir+"/slurmsubmit-%d-%d.sh" % (ifrom[i],ito[i]),"w") as f:
        f.write(replace_input(slurmtext,{"EXECNAME":execname,
					"ROOT":os.environ["PWD"]+"/..","JOBDIR":os.environ["PWD"]+"/"+jobdir,
                                         "iFROM":ifrom[i],"iTO":ito[i],"I":i,"NCPUS":NCPUs}));


