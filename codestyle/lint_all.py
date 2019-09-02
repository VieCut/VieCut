#! /usr/bin/python3
import glob
import subprocess

for folder in ['../app','../lib','../tests']:
    files = glob.glob(folder+'/**/*', recursive=True)

    for filename in files:
        if (filename.endswith(".cpp") or filename.endswith(".hpp") or filename.endswith(".h")):
            try:
                print("Linting file " + filename)
                out = subprocess.check_output(["./cpplint.py", "--filter=-build/c++11,-build/namespaces", filename])
                print("No errors!")
            except:
                print("ERRORS!")
                exit(1)
