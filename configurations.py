# parameters
# the overlap ratio cutoff to determine if two BLASTp hits on the same protein is considered overlapped
cutoffRatio = .60
evalueCutoff = 1e-4
# The amount of overlap two proteins need to have to be consider the "same" protein
simularProteinRatio = .90

numProteinInDataset = 1000

minNumProteinInFamily = 50

# delete generated folder from 1-3 to safe disk space
deleteFolders = False

# Location of runtime information
runTimefolder = "Results"
runTimeFile = "Results/$RuntimeInfo.txt"
logFolder = "Results"
logFile = "Results/log.txt"

# Location of saved progress (in the event of crashes)
savedProgressFolder = "Generated"
savedProgress = "Generated/savedProgress.cpickle"

# main core folders
fastaFolder = "Generated/0_FastaSequence"
proteinLenFolder = "Generated/1_ProteinLength"
blastdbFolder = "Generated/2_BlastDB"
alltoallFolder = "Generated/3_BlastAllToAll"
textResultsFolder = "Results/Text"
pickledResultsFolder = "Results/cpickle"

# Resources
FastaSeqDict = "Resources/FastaSeqDict.cpickle"
FamToArrDictLoc = "Resources/FamToLineDict.cpickle"
familiesFile = "Resources/familyCount.txt"
pFamToLabelFile = "Resources/pFamToLabels.txt"
PFamInfoByProteinFile = "Resources/PFamInfoByProtein.cpickle"

# Blast all to all
outputformat = 6
alltoallLogFolder = "Generated/log/BlastAllToAll"
blastdbLogFolder = "Generated/log/2_BlastDB"

