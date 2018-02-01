import configurations as conf
from cPickle import load, dump
import util


def generatePFamInfoByProtein():
    with open(conf.FamToArrDictLoc) as f:
        famToArrDict = load(f)
    PFamInfoByProtein = {}

    open(conf.pFamToLabelFile, "w")
    length = len(famToArrDict.keys())
    util.progressbarGuide(20)
    for i, family in enumerate(famToArrDict.keys()):
        util.progressbar(i, length, 20)
        arrs = famToArrDict[family]
        for arr in arrs:
            # format:
            # 0      1           2               3               4           5           6           7
            # PDB_ID CHAIN_ID    PdbResNumStart  PdbResNumEnd    PFAM_ACC    PFAM_Name   PFAM_desc   eValue
            pid = arr[0]
            start = arr[2]
            end = arr[3]
            if pid in PFamInfoByProtein.keys():
                PFamInfoByProtein[pid].append([i, start, end])
            else:
                PFamInfoByProtein[pid] = [[i, start, end]]
        with open(conf.pFamToLabelFile, "a") as f:
            f.write(str(i)+"\t"+family+"\n")
    with open(conf.PFamInfoByProteinFile, "wb") as f:
        dump(PFamInfoByProtein, f)


# assumes that PFamInfoByProtein is already generated
def correspondingPFamDict(moduleFamilyInfo):
    PFamInfoByProtein = load(open(conf.PFamInfoByProteinFile,"rb"))
    pFamDict = {}
    for protein in moduleFamilyInfo.keys():
        modules = PFamInfoByProtein[protein]
        pFamDict[protein] = modules

    return pFamDict
