import configurations as conf
from src import GenFasta, util
from src import blast_suite as blast
from src import build_HSPInt_graph as buildGraph
from src import defineBordersFromGraph as findBorders
from src import visualization as vis
from src import PFAMComparison as pfamComp
from cPickle import dump, load
import os
import datetime
import shutil
import time


def read_families():
    families = []
    with open(conf.familiesFile, "r") as f:
        for i, line in enumerate(f):
            # format:
            # 2840	C1-set
            arr = line.split("\t")
            numprot = int(arr[0].strip())
            fam_name = arr[1].strip()
            if numprot > 50:
                families.append((numprot, fam_name))
    return families


def runAlg(FamNames, filename):
    # fam_name=family[1]
    # numprot=family[0]

    # generate fasta sequence
    outfolder = conf.fastaFolder
    util.generateDirectories(outfolder)
    # filename = str(numprot) + "_" + fam_name

    print " Generating sequence file..."
    outdir = os.path.join(outfolder, filename)
    GenFasta.GenerateFastaInputForMultiFamilies(FamNames, outdir)

    # generate protein lengths
    plenFolder = conf.proteinLenFolder
    util.generateDirectories(plenFolder)
    plenDict = blast.generateProtLenDict(conf.fastaFolder, filename)
    dump(plenDict, open(os.path.join(plenFolder, filename), "wb"))

    # create blast databases
    print " Conducting BLASTp all-to-all..."
    blast.makeblastdb(conf.fastaFolder, filename)

    # conduct all to all BLASTp
    alltoallFolder = conf.alltoallFolder
    util.generateDirectories(alltoallFolder)
    blast.alltoallBlastP(conf.fastaFolder, filename,os.path.join(alltoallFolder, filename))

    # This is where my algorithm starts and also where I'll start timing
    print " Conducting my algorithm..."
    startTime = time.time()
    # build HSPIntGraph
    seqSimGraph, numBlastLines, numIntEdge = buildGraph.build_graph(filename, conf.alltoallFolder)

    # identify protein module borders
    # putative domains
    numModules, moduleFamilyInfo = findBorders.generatePutativeModules(seqSimGraph)
    putativeResult = vis.visualizeModuleFamilyInfo(moduleFamilyInfo)

    # remove submodules
    findBorders.removeSubModules(moduleFamilyInfo)
    moduleResult = vis.visualizeModuleFamilyInfo(moduleFamilyInfo)
    # print moduleResult

    # rename modules to have lower numbers
    numModulesAfterFilter = findBorders.renameModules(moduleFamilyInfo)
    moduleResultRenamed = vis.visualizeModuleFamilyInfo(moduleFamilyInfo)

    endTime = time.time()
    # calculate elapsed time
    timediff = endTime - startTime

    # output the results
    util.generateDirectories(conf.textResultsFolder)

    consizePath = os.path.join(conf.textResultsFolder, filename + "_Modules.txt")
    with open(consizePath, "w") as f:
        f.write(moduleResultRenamed)

    detailedPath = os.path.join(conf.textResultsFolder, filename + "_detailedResults.txt")
    with open(detailedPath, "w") as f:
        f.write("number of Blast Edges: "+str(numBlastLines)+"\n")
        f.write("time elapsed: "+str(timediff))
        util.printL("Completed in "+str(timediff)+" seconds\n")
        f.write("number of IntervalEdges added: "+str(numIntEdge)+"\n")
        f.write("Putative Domains: " + str(numModules) + "\n" + putativeResult + "\n")
        f.write("RemoveSubModules: \n" + moduleResult + "\n")
        f.write("Final Module Definition: " + str(numModulesAfterFilter) + "\n" + moduleResultRenamed + "\n")

    # Write down the timing results
    with open(conf.runTimeFile, "a") as f:
        f.write(str(numBlastLines)+"\t"+str(timediff)+"\t"+filename+"\n")

    # compare the borders with pfam definitions side by side
    pFamDict = pfamComp.correspondingPFamDict(moduleFamilyInfo)
    pfamCompPath = os.path.join(conf.textResultsFolder, filename + "_pFamSideBySide.txt")
    with open(pfamCompPath, "w") as f:
        f.write(vis.visualizePFamComparison(moduleFamilyInfo,pFamDict))

    # dump the border files for future comparison
    util.generateDirectories(conf.pickledResultsFolder)
    myBordersPath = os.path.join(conf.pickledResultsFolder, filename + "_myBorders.cpickle")
    pFamBordersPath = os.path.join(conf.pickledResultsFolder, filename + "_pfamBorders.cpickle")
    dump(moduleFamilyInfo, open(myBordersPath, "wb"))
    dump(pFamDict, open(pFamBordersPath, "wb"))

    # remove extra folders to safe disk space
    if conf.deleteFolders:
        shutil.rmtree(conf.proteinLenFolder)
        shutil.rmtree(conf.blastdbFolder)
        shutil.rmtree(conf.alltoallFolder)


def main():
    # load all families in PFam
    families = read_families()
    families.sort(key=lambda x: x[0])

    famNames = []
    famNumProt = []
    util.generateDirectories(conf.textResultsFolder)

    numProtCount = 0

    # create the time info file
    util.generateDirectories(conf.runTimefolder)

    progressNumProt = 0
    loadIndex = 0

    if os.path.isfile(conf.savedProgress):
        loadIndex, progressNumProt = load(open(conf.savedProgress))
        userPrompt = raw_input("progressed file detected, do you wish to start from index "+str(loadIndex)+"? (y/n):")
        if "y" not in userPrompt and "Y" not in userPrompt:
            loadIndex = 0
            progressNumProt = 0

    # open the timing file if this is the first time we run this program
    if loadIndex == 0:
        with open(conf.runTimeFile, "w") as f:
            f.write("#BLAST Edges\ttime\tsampleName\n")

    # calculate total number of proteins in all of analysis (for progress searching)
    totNumProt = 0
    for family in families:
        if family[0] >= conf.minNumProteinInFamily:
            totNumProt += family[0]

    print "len(families):", len(families)
    for j in range(loadIndex, len(families)):
        family = families[j]
        if family[0] >= conf.minNumProteinInFamily:
            numProtCount += family[0]
            progressNumProt += numProtCount
            famNames.append(family[1])
            famNumProt.append(family[0])

            if numProtCount > conf.numProteinInDataset:
                filename = datetime.datetime.now().strftime("%Y%m%d_%I%p") + "_" + "%05d" % j + "_" + famNames[0]
                filedir = os.path.join(conf.textResultsFolder, filename + "_FamilyInfo.txt")
                with open(filedir, "w") as f:
                    f.write("NumProt\tFamily Name")
                    for i, famName in enumerate(famNames):
                        f.write(str(famNumProt[i])+"\t"+famName+"\n")

                # run my algorithm on the set of families
                util.printL("\nRunning Custom Dataset: "+filename+"\n")
                util.printL("Progress: "+str(int(float(progressNumProt)*100/totNumProt)/float(100))+" %")
                util.printL(" ( Family info located in: "+filedir+")\n")
                runAlg(famNames, filename)
                numProtCount = 0
                famNames = []
                famNumProt = []

                # save the starting point of the next experiment
                if j < len(families):
                    util.generateDirectories(conf.savedProgressFolder)
                    dump((j+1, progressNumProt), open(conf.savedProgress, "wb"))

    #runAlg([families[1], families[2]])

    # famNames = ["Neur_chan_memb"]
    # filename = datetime.datetime.now().strftime("%Y%m%d_%I%p") + "_" + famNames[0]
    # runAlg(famNames, filename)


def test():
    pfamComp.generatePFamInfoByProtein()


if __name__ == '__main__':
    # test()
    main()

