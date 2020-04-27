#! /usr/bin/python3
import glob
import os

cmdgit = "git status | grep modified | awk '{print $2}'"
files = os.popen(cmdgit).read().splitlines()

for filename in files:
    if (filename.endswith(".cpp") or filename.endswith(".hpp") or filename.endswith(".h")):
        try:
            print("Linting file " + filename)
            cmdlint = "./cpplint.py --filter=-build/c++11,-build/namespaces " + filename
            process = os.popen(cmdlint)
            output = process.read()
            print("No errors!")
        except:
            print("ERRORS!")
