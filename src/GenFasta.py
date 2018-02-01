
import configurations as conf
from cPickle import load
import util


def toFastaString(desc, seq):
    return ">" + str(desc) + "\n" + str(seq) + "\n\n"


# assumes FamilyToArr dict is generated
# assumes pidToBioSeqDict is also generated
def GenerateFastaInputForGivenFamily(famName, outfiledir):
    # import the FamToArrDict
    with open(conf.FamToArrDictLoc, "rb") as f:
        FamToArrDict = load(f)
    with open(conf.FastaSeqDict, "rb") as f:
        pidToBioSeqDict = load(f)

    pfamRows = FamToArrDict[famName]

    proteinIDs = []
    for row in pfamRows:
        # format:
        # 0      1           2               3               4           5           6           7
        # PDB_ID	CHAIN_ID	PdbResNumStart	PdbResNumEnd	PFAM_ACC	PFAM_Name	PFAM_desc	eValue

        proteinIDs.append(row[0])

    # remove duplicates
    proteinIDs = list(set(proteinIDs))

    # write the sequences to file
    open(outfiledir, "w")
    for pid in proteinIDs:
        with open(outfiledir, "a") as f:
            f.write(toFastaString(pid, pidToBioSeqDict[pid]))


# assumes FamilyToArr dict is generated
# assumes pidToBioSeqDict is also generated
def GenerateFastaInputForMultiFamilies(famNames, outfiledir):
    # import the FamToArrDict
    with open(conf.FamToArrDictLoc, "rb") as f:
        FamToArrDict = load(f)
    with open(conf.FastaSeqDict, "rb") as f:
        pidToBioSeqDict = load(f)

    pfamRows = []
    for famName in famNames:
        pfamRows += FamToArrDict[famName]

    proteinIDs = []
    for row in pfamRows:
        # format:
        # 0      1           2               3               4           5           6           7
        # PDB_ID	CHAIN_ID	PdbResNumStart	PdbResNumEnd	PFAM_ACC	PFAM_Name	PFAM_desc	eValue

        proteinIDs.append(row[0])

    # remove duplicates
    proteinIDs = list(set(proteinIDs))

    # write the sequences to file
    open(outfiledir, "w")
    for pid in proteinIDs:
        with open(outfiledir, "a") as f:
            if pid in pidToBioSeqDict.keys():
                f.write(toFastaString(pid, pidToBioSeqDict[pid]))
            else:
                util.printL("protein, "+pid+", is not found in the fasta sequence provided.\n")


def main():
    GenerateFastaInputForGivenFamily("Stap_Strp_toxin", "Stap_Strp_toxin.fasta")


if __name__ == '__main__':
    main()
