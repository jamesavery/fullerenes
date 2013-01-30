#!/usr/bin/python
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
    
for i in range(N):
    s = i*chunksize+start;
    e = (i+1)*chunksize+start-1;

    dirname = jobdir+"/"+str(i);
    mkdirp(dirname);

    with open(dirname+"/input-"+str(s)+"-"+str(e)+".inp","w") as f:
        f.write(replace_input(text,{"FROM":s,"TO":e}));

    with open(dirname+"/llsubmit.sh","w") as f:
        f.write(replace_input(lltext,{"PWD":os.environ["PWD"],"WORKDIR":os.environ["PWD"]+"/"+dirname}));


dirname = jobdir+"/"+str(N);
mkdirp(dirname);

with open(dirname+"/input-"+str(e+1)+"-"+str(end)+".inp","w") as f:
    f.write(replace_input(text,{"FROM":e+1,"TO":end}));
