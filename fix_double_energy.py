# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 17:46:12 2021

@author: chernir
"""
import os
 
#dirname = r"D:\BioXAS\testdata"
dirname = "."  # current dir
fixedname = "dEfixed"
dEmin = 0.1  # minimal E step in eV

for root, dirs, files in os.walk(dirname, topdown=False):
    for name in files:
        filepath = os.path.join(root, name)

        filename = os.path.basename(filepath)
        baseDirName = os.path.dirname(filepath)
        if str(baseDirName).endswith(fixedname):
            continue
             
        file1 = open(filepath, 'r')
        Lines = file1.readlines()
        file1.close()

        if not Lines[0].startswith("#"):  # not an XAS scan
            continue

        newLines = []

        Ep = None
        En = None
        for tline in Lines:
            if len(tline) < 2:
                continue
            if not str(tline).startswith("#"):
                fp = tline.split()[0]
                try:
                    En = float(fp)
                except:
                    print("Can't convert to float", fp)
                    continue
                if all([En, Ep]):
                    dE = abs(En-Ep)
                    if dE < dEmin:
                        print("Dropping repeating energy point at", En, "eV")
                        continue
                    else:
                        Ep = En
                else:
                    Ep = En
            newLines.append(tline)
                
#        print(len(Lines), type(Lines))
        fixedDirName = os.path.join(baseDirName, fixedname)
        if not os.path.exists(fixedDirName):
            os.mkdir(fixedDirName)
        file2 = open(os.path.join(fixedDirName, filename), 'w')
        file2.writelines(newLines)
        file2.close()
#        break
