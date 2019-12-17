#! /usr/bin/python3
import glob
import subprocess
import os

files = os.popen("git status | grep modified | awk '{print $2}'").read().splitlines()

for filename in files:
    if (filename.endswith(".cpp") or filename.endswith(".hpp") or filename.endswith(".h")):
        try:
            print("Linting file " + filename)
            out = subprocess.check_output(["./cpplint.py", "--filter=-build/c++11,-build/namespaces", filename])
            print("No errors!")
        except:
            print("ERRORS!")
