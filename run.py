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


def runAlg(filename):
    # fam_name=family[1]
    # numprot=family[0]

    # # generate fasta sequence
    # outfolder = conf.fastaFolder
    # util.generateDirectories(outfolder)
    # # filename = str(numprot) + "_" + fam_name
    #
    # print " Generating sequence file..."
    # outdir = os.path.join(outfolder, filename)
    # GenFasta.GenerateFastaInputForMultiFamilies(FamNames, outdir)

    if not conf.skipBLAST:
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
        blast.alltoallBlastP(conf.fastaFolder, filename, os.path.join(alltoallFolder, filename))

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
    findBorders.removeSuperModules(moduleFamilyInfo)
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
        f.write("Putative Modules: " + str(numModules) + "\n" + putativeResult + "\n")
        f.write("RemoveSuperModules: \n" + moduleResult + "\n")
        f.write("Final Module Definition: " + str(numModulesAfterFilter) + "\n" + moduleResultRenamed + "\n")

    # Write down the timing results
    with open(conf.runTimeFile, "a") as f:
        f.write(str(numBlastLines)+"\t"+str(timediff)+"\t"+filename+"\n")

    # # compare the borders with pfam definitions side by side
    # pFamDict = pfamComp.correspondingPFamDict(moduleFamilyInfo)
    # pfamCompPath = os.path.join(conf.textResultsFolder, filename + "_pFamSideBySide.txt")
    # with open(pfamCompPath, "w") as f:
    #     f.write(vis.visualizePFamComparison(moduleFamilyInfo,pFamDict))
    #

    # dump the border files for future comparison
    util.generateDirectories(conf.pickledResultsFolder)
    myBordersPath = os.path.join(conf.pickledResultsFolder, filename + "_myBorders.cpickle")
    # pFamBordersPath = os.path.join(conf.pickledResultsFolder, filename + "_pfamBorders.cpickle")
    dump(moduleFamilyInfo, open(myBordersPath, "wb"))
    # dump(pFamDict, open(pFamBordersPath, "wb"))

    # remove extra folders to safe disk space
    if conf.deleteFolders:
        shutil.rmtree(conf.proteinLenFolder)
        shutil.rmtree(conf.blastdbFolder)
        shutil.rmtree(conf.alltoallFolder)


def main():
    seqFiles = os.listdir(conf.inputFolder)
    # generate runtime output file
    util.generateDirectories(conf.runTimefolder)
    open(conf.runTimeFile, "w")
    for i, seqFile in enumerate(seqFiles):
        util.printL("\nAnalyzing "+seqFile+"("+str(i)+"/"+str(len(seqFiles))+")\n")
        runAlg(seqFile)


def test():
    pfamComp.generatePFamInfoByProtein()


if __name__ == '__main__':
    # test()
    main()

