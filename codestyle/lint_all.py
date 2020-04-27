#! /usr/bin/python3
import glob
import os

for folder in ['../app','../lib','../tests']:
    files = glob.glob(folder+'/**/*', recursive=True)

    for filename in files:
        if (filename.endswith(".cpp") or filename.endswith(".hpp") or filename.endswith(".h")):
            print("Linting file " + filename)
            process = os.popen("./cpplint.py --filter=-build/c++11,-build/namespaces " + filename)
            output = process.read()
