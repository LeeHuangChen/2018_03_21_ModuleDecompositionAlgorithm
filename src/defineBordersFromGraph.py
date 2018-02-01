#import os
#from cPickle import dump, load, HIGHEST_PROTOCOL
#import configurations as conf
import networkx as nx
import util


def overlap(s1start, s1end, s2start, s2end):
    maxstart = max(s1start, s2start)
    minend = min(s1end, s2end)
    if maxstart < minend:
        return minend - maxstart
    else:
        return 0


# borders is a list of list consists of elements: [ModuleID,start,end]
# this returns true if there is a border that intersects with another border in borders
def containsOverlapBorders(borders):
    for i in range(len(borders)):
        for j in range(i + 1, len(borders)):
            if borders[i][0] == borders[j][0]:
                overlapped = overlap(borders[i][1], borders[i][2], borders[j][1], borders[j][2]) > 0
                if overlapped:
                    return (i, j)
    return False


def collapsOverlappingBorders(borders):
    overlapInfo = containsOverlapBorders(borders)
    while overlapInfo != False:
        i = overlapInfo[0]
        j = overlapInfo[1]
        border1 = borders[i]
        border2 = borders[j]
        newstart = min(border1[1], border2[1])
        newend = max(border1[2], border2[2])
        # newborder=[border1[0],newstart,newend]

        # delete the overlapped border and add the merged border
        deleteindex = None
        otherindex = None
        if i > j:
            deleteindex = i
            otherindex = j
        else:
            deleteindex = j
            otherindex = i

        del borders[deleteindex]

        borders[otherindex][1] = newstart
        borders[otherindex][2] = newend
        # border.append(newborder)

        overlapInfo = containsOverlapBorders(borders)


def generatePutativeModules(g):
    # # generate the directories
    # util.generateDirectories(borderInfodir)
    # util.generateDirectories(borderResultdir)

    # with open(os.path.join(hspIntGraphdir, graphFile)) as fin:
    #     g = load(fin)

    # find the connected components of the graph:
    CCgraphs = list(nx.connected_component_subgraphs(g))
    CCgraphs.sort(key=lambda tup: len(tup.nodes()))

    # create a dictionary where key=protein name, val:border information
    modulefamilyinfo = {}
    for moduleID, CCgraph in enumerate(CCgraphs):
        for node in CCgraph.nodes():

            proteinName = node[0]
            start = int(node[1])
            end = int(node[2])

            if proteinName not in modulefamilyinfo:
                modulefamilyinfo[proteinName] = [[moduleID, start, end]]
            else:
                # if we already have information on this protein
                modulefamilyinfo[proteinName].append([moduleID, start, end])
                collapsOverlappingBorders(modulefamilyinfo[proteinName])

    # cleanup borders that ended up developing to overlap borders
    for protein in modulefamilyinfo.keys():
        collapsOverlappingBorders(modulefamilyinfo[protein])

    return len(CCgraphs), modulefamilyinfo


def checkForSubModules(modules):

    for i, M1 in enumerate(modules):
        for j, M2 in enumerate(modules):
            if i != j:
                s1 = M1[1]
                e1 = M1[2]

                s2 = M2[1]
                e2 = M2[2]

                M1_in_M2 = (s2 <= s1) and (e2 >= e1)
                M2_in_M1 = (s1 <= s2) and (e1 >= e2)

                if M1_in_M2:
                    return i  # returns the index of the module that is to be removed
                if M2_in_M1:
                    return j  # returns the index of the module that is to be removed
    return -1


def removeSubModules(moduleFamilyInfo):
    for protein in moduleFamilyInfo.keys():
        modules = moduleFamilyInfo[protein]
        removeIndex = checkForSubModules(modules)
        while removeIndex != -1:
            modules.pop(removeIndex)
            removeIndex = checkForSubModules(modules)


# renames all of the module labels from 1 so the numbers are smaller
# returns the number of modules identified
def renameModules(moduleFamilyInfo):
    moduleCounter = 1
    moduleDict = {}
    for protein in moduleFamilyInfo.keys():
        modules = moduleFamilyInfo[protein]
        for module in modules:
            moduleName = module[0]
            if moduleName in moduleDict.keys():
                module[0] = moduleDict[moduleName]
            else:
                module[0] = moduleCounter
                moduleDict[moduleName] = moduleCounter
                moduleCounter += 1

    return moduleCounter - 1
