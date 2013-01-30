#!/usr/bin/python
import sys, string;

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

def replace_input(text,keyvalues):
    for k in keyvalues:
        text = text.replace("$%s"%k,str(keyvalues[k]))
    return text;

text = ();
with open(inputfile,'r') as f:
    text = f.read();
    
for i in range(N):
    s = i*chunksize+start;
    e = (i+1)*chunksize+start-1;

    with open(jobdir+"/"+str(s)+"-"+str(e)+".inp","w") as f:
        f.write(replace_input(text,{"FROM":s,"TO":e}));

with open(jobdir+"/"+str(e+1)+"-"+str(end)+".inp","w") as f:
    f.write(replace_input(text,{"FROM":e+1,"TO":end}));



    
