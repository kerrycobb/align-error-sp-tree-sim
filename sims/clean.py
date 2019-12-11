#!/usr/bin/env python3

import os 
import glob
import shutil
import fire

def check(dir, file):
    path = os.path.join(dir, file)
    if os.path.exists(path):
        state = True
    else:
        state = False
        print("{} missing".format(path))
    return state

def clean(dir):
    states = [] 
    states.append(check(dir, "coverage-theta.csv"))
    states.append(check(dir, "coverage-time.csv"))
    states.append(check(dir, "summary-theta.csv"))
    states.append(check(dir, "summary-time.csv"))
    if False not in states:
        for i in glob.glob(os.path.join(dir, "rep-*")):
            shutil.rmtree(i) 
    print("{} cleanup complete".format(dir))
       
if __name__ == "__main__":
    fire.Fire(clean)