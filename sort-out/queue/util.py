#!/usr/bin/env python
import os;

def replace_input(text,keyvalues):
    for k in keyvalues:
        text = text.replace("$%s"%k,str(keyvalues[k]))
    return text;

def mkdirp(path):
    try:
        os.makedirs(path);
    except:
        pass;

def readfile(filename):
    with open(filename,"r") as f:
        return f.read();
