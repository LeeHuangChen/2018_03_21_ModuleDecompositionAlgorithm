import os
import stat
import sys
import configurations as conf


# generate all the directories needed for the given path (helper function)
def generateDirectories(path):
    folders = path.split("/")
    curdir = ""
    for folder in folders:
        curdir = os.path.join(curdir, folder)
        if not os.path.exists(curdir):
            os.mkdir(curdir)


def generateDirectoriesMult(paths):
    for path in paths:
        generateDirectories(path)


def printL(string):
    generateDirectories(conf.logFolder)
    with open(conf.logFile, "a") as f:
        f.write(string)
        sys.stdout.write(string)
        sys.stdout.flush()


def progressbar(i, length, numberNotification):
    scale = length / numberNotification
    if scale > 0:
        if i % scale == 0:
            sys.stdout.write('*')
            if length - i < scale:
                sys.stdout.write('\n')
            sys.stdout.flush()


def progressbarGuide(length):
    sys.stdout.write('|')
    for i in range(0, length - 1):
        sys.stdout.write('-')
    sys.stdout.write('|\n')
    sys.stdout.flush()


def createEmptyFiles(paths):
    for path in paths:
        f = open(path, "w")
        f.close()


def testPath(path):
    if not os.path.exists(path):
        os.mkdir(path)


def appendFile(path, content):
    f = open(path, "a")
    f.write(content)
    f.close()


def readFile(path):
    f = open(path, "r")
    content = f.read()
    f.close()
    return content


# processFusedGenes functions

# takes the protein, taxa, and sequence information and produce a FASTA format string for that sequence
def toFASTA(prot, taxa, seq):
    return ">" + prot + " [" + taxa + "]\n" + seq + "\n\n"


# A custom parsing function that uses different delims for each entry
# input: a list of delimiters, taken in order
# output: a list of string version of each values
def custParse(inStr, delims):
    curStr = inStr
    outList = []
    for delim in delims:
        delLoc = curStr.find(delim)
        value = curStr[0:delLoc]
        if len(value) > 0:
            outList.append(value)
        # remove the value and delim from the curStr
        curStr = curStr[delLoc + len(delim):]
    if len(curStr) > 0:
        outList.append(curStr)
    return outList


# cpickle write data to file
def dataToFile(data, filename):
    with open(filename, 'wb') as fout:
        dump(data, fout, HIGHEST_PROTOCOL)


# cpickle read file to data
def fileToData(filename):
    data = load(open(filename, "r"))
    return data


def percentSimple(i, length, numberNotification):
    scale = length / numberNotification
    if (i % scale == 0):
        print int(float(i) / float(length) * 100), "%"


def percent(i, length, numberNotification, header="", footer="", percentRange=(0, 100)):
    scale = length / numberNotification
    if (i % scale == 0):
        progress = str(
            percentRange[0] + int(float(i) / float(length) * (percentRange[1] - percentRange[0]) * 100) / float(
                100)) + "%"
        write = header + progress + footer
        print write
