# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 20:35:55 2018

@author: HuangMing
"""
import pandas as pd
import numpy as np
import os
import itertools
import random

import AtomClass as Atom
import MoleculeClass as Mol
import BondClass as Bond
import readMol2

def InfoInput(n, a, AtomInfo):
    atom = AtomInfo
    if n == 0:
        a.setIndex(atom[0])
    elif n == 1:
        a.setAtomName(atom[1])
    elif n == 2:
        pos = [atom[2], atom[3], atom[4]]
        a.setPos(pos)
    elif n == 3:
        pass
    elif n == 4:
        pass
    elif n == 5:
        a.setAtomType(atom[5])
    elif n == 6:
        a.setSubID(atom[6])
    elif n == 7:
        a.setSubName(atom[7])
    elif n == 8:
        a.setCharge(atom[8])

def BondInfoInput(n, b, BondInfo):
    bond = BondInfo
    if n == 0:
#        print(bond[0])
        index = bond[0]
        atom1 = bond[1]
        atom2 = bond[2]
        b.setIndex(index)
        b.setAtom(atom1, atom2)
        b.setBondType(bond[3])
    else:
        pass
        
def InfoOutput(a):
    index = a.getIndex()
    name = a.getAtomName()
    pos = a.getPos()
    atomType = a.getAtomType()
    subID = a.getSubID()
    subName = a.getSubName()
    charge = a.getCharge()
    content = [index, name, pos, atomType, subID, 
               subName, charge]
    for i in range(len(content)):
        if i == 2:
            print("X= %s, Y= %s, Z=%s" % (pos[0], pos[1], pos[2]))
        else:
            print(content[i])   

def AtomToMolecule(index, start, end, name, atomList): #update atoms content to the mol list

    mol = Mol.MoleculeInfo()
    mol.setIndex(index)
    mol.setName(name)
    for i in range (start, end):
        mol.setAtoms(atomList[i])
    return mol
        
def SortSubID(atomList): #split atoms to molecule
    ini_ID = atomList[0].getSubID()
    MoleculeList = []
    index = 0
    ini = 0
    for i in range (len(atomList)):
        subID = atomList[i].getSubID()
        if subID == ini_ID:
            if i == len(atomList) - 1:
#                print("tst-1: should generate a new mol here", i)
#                print(index)
                index +=1
                name = "mol" +str(index)
                mol = AtomToMolecule(index, ini, len(atomList), name, atomList)
                MoleculeList.append(mol)
#                print("MOl length: ", len(MoleculeList))
                return MoleculeList
            else:
                continue
        else:
            index +=1
            name = "mol" + str(index)
            mol = AtomToMolecule(index, ini, i, name, atomList)
            MoleculeList.append(mol)
            ini = i
            ini_ID = subID
    return MoleculeList

def MonCros(monLen, crosLen, moleculeList): #separate molecules into monomer and crosslinker
    monList = []
    crosList = []
    
    for i in range (len(moleculeList)):
        if len(moleculeList[i].getAtoms()) == monLen:
            monList.append(moleculeList[i])
        elif len(moleculeList[i].getAtoms()) == crosLen:
            crosList.append(moleculeList[i])
        else:
            print("Wired molecules doesn't belong neither monomer or crosslinker!")
    return monList, crosList

def SetReact(list, index): #list no matter monomer list or cros list, index includes all the reactive atom
    #label all the react atoms in each molecule
    for i in range (len(list)):
        start = list[i].getAtoms()[0].getIndex()
        for ii in range (len(index)):
            rxIndex = str(int(start)+int(index[ii]))
#            print('{},{}'.format(start, rxIndex))
            list[i].setReact(rxIndex)

def ExportMOL2(name, outputName, content, atomList, BondList):
    f = open(outputName, 'w')
    f.write('@<TRIPOS>MOLECULE\n')
    f.write(name)
    tmp_str = '\n' + content[0]+ '\n' + content[1] + '\n' + content[2] + '\n'
    f.write(tmp_str)
    f.write('\n@<TRIPOS>ATOM\n')
    
    for i in range (len(atomList)):
        atom = atomList[i]
        str1 = atom.outputData()
        f.write(str1)
    f.write('@<TRIPOS>BOND\n')
    for ii in range (len(BondList)):
        bond = BondList[ii]
#        print('bond: ', bond)
        str1 = bond.outputData()

        f.write(str1)
    f.write("\n")   
def GetLabelAtom(list):
    reactIndex = []
    for i in range(len(list)):
        for ii in range(len(list[i].getReact())):
            reactIndex.append(list[i].getReact()[ii]) #Structure for reactIndex: [molIndex, [react atoms index]]
    return reactIndex
def CalDist(aPos, bPos):
    x2 = np.power(float(aPos[0])-float(bPos[0]),2)
    y2 = np.power(float(aPos[1])-float(bPos[1]),2)
    z2 = np.power(float(aPos[2])-float(bPos[2]),2)
    dist = np.sqrt(x2+y2+z2)
    return dist
def GenBond(atom1, atom2, atomList, bondList, molList, bondLimit):
    index = len(bondList) + 1
    atom1Index = atom1.getIndex()
    atom2Index = atom2.getIndex()
    atom1Bond = CountBond(atom1Index, bondList) 
    atom2Bond = CountBond(atom2Index, bondList)
#    print("atoms we want to bond with: {:}, {:}".format(atom1Index, atom2Index))
    #print("mol list", molList)
    cond1 = Criteria(atom1Bond, atom1.getAtomName(), atomList)#atom1, atom2 which index starts from 0
    cond2 = Criteria(atom2Bond, atom2.getAtomName(), atomList)
    cond3 = Criteria2(atom1, atom2, molList, bondList)
    cond4 = Criteria3(bondLimit, bondList)
    print(cond1, cond2, cond3, cond4, len(bondList))
    if cond1 and cond2 and cond3 and cond4:
        b = Bond.BondsInfo()
        b.setIndex(index)
        b.setAtom(atom2.getIndex(), atom1.getIndex())
        b.setBondType(1) #Bondtype 1 means C3 bond
        bondList.append(b)
    
def CountBond(atomIndex, bondList):
    bonds = []
    tmp = []
    for i in range (len(bondList)):
        tmp.append(bondList[i].getAtom())
    for it in range (len(tmp)):
        for iit in range(2):
            bonds.append(tmp[it][iit])
    num = bonds.count(str(atomIndex))
    return num

def ChkAtomBonds(atomIndex, bondList): #check atoms generate which bonds. it will return the bonds index with this atom
    indexes = []
    bonds = []
    for i in range(len(bondList)): #bondList here just a list of bonds not bonds class
        bonds.append(bondList[i].getAtom())

    for i in range(len(bonds)):
        if atomIndex in bonds[i]:
            indexes.append(i)
    return indexes #returns a bond index

def ChkAtomMol(atomIndex, molList):
    molIndex = 0
    for i in range(len(molList)):
        atoms = []
        for ii in range (len(molList[i].getAtoms())):
            atoms.append(molList[i].getAtoms()[ii].getIndex())
        if atomIndex in atoms:
            molIndex = i
        else:
            continue
    return molIndex

def Criteria(bondsNum, atomName, atomList): #Bonds limitation
    if atomName == 'C':
        if bondsNum >=3:
            return False
    if atomName == 'N':
        if bondsNum >=2:
            return False
    if atomName == 'O':
        if bondsNum >=1:
            return False
    else:
        return True
    
def Criteria2(atom1, atom2, molList, bondList): 
#I wonder to make a change in this criteria
#Parameters left: atom1 is the atom which will be the head of the generating bond, and atom2 is the end
#
#This should become:
#1. From atom1, we get its molecule index, and should also get the reactive atoms on that molecules.
#2. Then we can check, if atom2 is belongs to any of them, if it is then this bond will not be generated


    molIndex = []
    indexMol = []
    #atomIndex includes all atoms bonded with atom1
    atom1BondAtoms = bondAtoms(atom1, bondList)
    for i in range (len(atom1BondAtoms)):
        indexMol.append(ChkAtomMol(atom1BondAtoms[i], molList))
        for ii in range (len(indexMol)):
            molIndex.append(indexMol[ii])
    atomMolIndex = []
    #atomsCrit = [atom1.getIndex()] + [i for i in atomsIndex]
#    print("Tst-AtomsIndex which we want to check: ", atom1BondAtoms)
    atom2Mol = ChkAtomMol(atom2.getIndex(), molList)
    ############Just for test to clearly see which mol does atoms belongs to#########
    tmp = atom1BondAtoms
    for i in range(len(tmp)):
        atomMolIndex.append(ChkAtomMol(tmp[i], molList))
#    print("Tst-atomMolList: ", atomMolIndex)
#    print("Atom2 mol#: ", atom2Mol)
    if atom2Mol in atomMolIndex:
#        print("Repeat mol, no good")
        return False
    else:
#        print("Good")
        return True
    #############Finish test################
    
#    if atom2.getIndex() in atomsCrit:
#        return False
#    else:
#        return True
    
def Criteria3(bondLimit, bondList):
    print("Bonds Limit: ", bondLimit, "Bonds number now: ", len(bondList))
    if len(bondList) == bondLimit:
        return False
    else:
        return True

def bondAtoms(atom, bondList):##return atoms index which connect to the bond
    bonds = []
    atoms = []
    for i in range(len(bondList)):
        bonds.append(bondList[i].getAtom())
    
    bondIndex = ChkAtomBonds(atom.getIndex(), bondList) #Index starts from 0
#    print("Chk BondIndex: ", bondIndex)
    for i in range (len(bondIndex)):
        atoms.append(bonds[bondIndex[i]]) #merge bonds' atoms together
    
    #print("atomsCK: ", len(atoms)) #atoms list [[atoms1, atoms2],[atoms1, atoms2s]]
    atomsIndex = list(itertools.chain.from_iterable(atoms))
#    print("Chk AtomIndex: ", atomsIndex)
#    print("Delete Index: ", atom.getIndex())
    atomsIndex = [i for i in atomsIndex if i != atom.getIndex()] #remove selected atoms
#    print("After delete atoms: ", atomsIndex) 
    return atomsIndex
#    for i in range(len(molList)):
#        atomsList = [] #Include all atoms in mol[i]
#        atoms = molList[i].getAtoms()
#        for ii in range(len(atoms)):
#            atomsList.append(atoms[ii].getIndex()) #all indexes which mol1's atoms
#        if any([atomTmp for atomTmp in atomsList if atomTmp in atomsIndex])
#            if atom2.getIndex in atomsList:
#                return False
#        else:
#            return True
    
def Crosslink(cutoff, atomList, bondList, molList, monList, crosList, bondsNum):
    mon = GetLabelAtom(monList)
    cros = GetLabelAtom(crosList) #I suppose that crosslinker govern the crosslinking process
    random.shuffle(mon)
    random.shuffle(cros)
    print("Bonds length: ", len(bondList))
    bondLimit = len(bondList) + bondsNum
    for i in range(len(cros)):
        for ii in range(len(mon)):
            crosIndex = int(cros[i])-1
            if crosIndex >= len(atomList):
                break
            else:                
                monIndex = int(mon[ii])-1
                atomCross = atomList[crosIndex]
                atomMon = atomList[monIndex]
                crosCoord = atomCross.getPos()
                monCoord = atomMon.getPos()
                dist = CalDist(crosCoord, monCoord)
                if dist <= cutoff:
                    GenBond(atomCross, atomMon, atomList, bondList, molList, bondLimit)
                else:
                    continue
                
def main(fileName, outputName, monLen, crosLen, monR, crosR, cutoff, bondsNum, cycle, monoNum, crosNum):    
    info = readMol2.InfoInput(fileName, monLen, crosLen, cycle, monoNum, crosNum)
    name = info[5][0]
    basicPara = info[5][1:]
    SetReact(info[3], monR)
    SetReact(info[4], crosR)
    Crosslink(cutoff, info[0], info[1], info[2], info[3], info[4], bondsNum)
    
    atoms = len(info[0])
    bonds = len(info[1])
    content_1 = '{:} {:} {:} {:} {:}'.format(str(atoms), str(bonds), str(0), str(0), str(0))
    basicPara[0] = content_1
    print("tst-2: ",basicPara)
    ExportMOL2(name, outputName, basicPara, info[0], info[1])
    
    return info[0], info[1]

#monR = [3, 21]
#crosR = [0, 14]
#a = main('system.mol2', 'tst.mol2', 25, 15, monR, crosR, 15., 2)
#newInfo = readMol2.InfoInput('sys-bond.mol2', 25, 15)
#newAtoms = newInfo[0]
#
#readMol2.AtomInfoUpdate(a[0], newAtoms, pos=True)
#ExportMOL2('123123', 'tst-2.mol2', 'tst', a[0],a[1])
