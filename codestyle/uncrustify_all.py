#! /usr/bin/python3

import glob
import subprocess

for folder in ['../app','../lib','../tests']:
    files = glob.glob(folder+'/**/*', recursive=True)

    for filename in files:
        if (filename.endswith(".cpp") or filename.endswith(".hpp") or filename.endswith(".h")):
            subprocess.run(["uncrustify", "-c", "uncrustify.cfg", "--no-backup", filename])