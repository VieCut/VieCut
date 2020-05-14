#! /usr/bin/python3

import glob
import os

for folder in ['../app','../lib','../tests']:
    files = glob.glob(folder+'/**/*', recursive=True)

    for filename in files:
        if (filename.endswith(".cpp") or filename.endswith(".hpp") or filename.endswith(".h")):
            uncrustcmd = "uncrustify -c uncrustify.cfg --no-backup --replace " + filename
            process = os.popen(uncrustcmd)
            output = process.read()
