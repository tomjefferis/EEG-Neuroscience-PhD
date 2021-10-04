import os

mainPath = "W:\PhD\PatternGlareData\participants"
part = ""
for i in range(1,40):
    part = str(i) 
    partpath = mainPath + "\participant_" + part   
    os.mkdir(partpath)
    os.rename(mainPath + "\spmeeg_P"+part+"_075_80Hz.dat", partpath + "\spmeeg_P"+part+"_075_80Hz.dat")
    os.rename(mainPath + "\spmeeg_P"+part+"_075_80Hz.mat", partpath + "\spmeeg_P"+part+"_075_80Hz.mat")
    
