#!/usr/bin/env python3

import os 
import glob
import shutil
import re
import yaml
import fire


def clean(dir):
    for i in glob.glob(os.path.join(dir, "rep-*")):
        try:
            os.remove(os.path.join(i, "alignment.nex"))
        except:
            pass
        try:
            os.remove(os.path.join(i, "starbeast.xml"))
        except:
            pass
        for j in glob.glob(os.path.join(i, "ecoevo-chain-*")):
            try:
                shutil.rmtree(j)
            except:
                pass
        for j in glob.glob(os.path.join(i, "starbeast-chain-*")):
            try:
                shutil.rmtree(j)
            except:
                pass

if __name__ == "__main__":
    fire.Fire(clean)