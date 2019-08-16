#!/usr/bin/env python3

import os 
import glob
import shutil
import re
import yaml

out_dirs = glob.glob("out-*")



for i in out_dirs:
    for j in range(0, 100):
        rep = str(j).zfill(3)
        rep_dir = os.path.join(i, "seed1-reps100", "replicate-" + rep)
        shutil.copy("configs/eco-config.yml", rep_dir)