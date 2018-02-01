import os


def forAllLineInFile(fileInfo, doworkwithline):
    # fileInfo=(fileIndex, inputFolder, inputfile)
    inputFolder = fileInfo[1]
    inputfile = fileInfo[2]

    filepath = os.path.join(inputFolder, inputfile)
    with open(filepath, "r") as f:
        for i, line in enumerate(f):
            lineInfo = (i, line)
            doworkwithline(fileInfo, lineInfo)


def forAllFiles(doworkwithfile, folder):
    inputfiles = os.listdir(folder)

    for i, inputfile in enumerate(inputfiles):
        fileInfo = (i, folder, inputfile)
        doworkwithfile(fileInfo)


def forAllLinesInAllFiles(doworkonline, folder):
    forAllInputFiles(lambda fileInfo: forAllLineInFile(fileInfo, doworkonline), folder)
