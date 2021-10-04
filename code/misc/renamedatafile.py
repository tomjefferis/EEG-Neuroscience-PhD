import os
from shutil import copyfile

mainPath = "A:\PhD\PatternGlareData\participants"
part = ""
for i in range(1,40):
    part = str(i) 
    partpath = mainPath + "\participant_" + part

    if(i<10):
        part = "0"+ part
    

    if(os.path.exists(partpath)):
        copyfile(partpath + "\spmeeg_P"+part+"_075_80Hz.dat", partpath + "\spmeeg_P"+part+"_075_80Hz_rejected.dat")
        copyfile(partpath + "\spmeeg_P"+part+"_075_80Hz.mat", partpath + "\spmeeg_P"+part+"_075_80Hz_rejected.mat")
    
