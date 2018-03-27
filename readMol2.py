# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 08:44:52 2018

@author: HuangMing
"""

import pandas as pd
import numpy as np
import sys
import os
import itertools
import random
import copy

import AtomClass as Atom
import BondClass as Bond
import MoleculeClass as Mol



def GetIndex(string, df):
    idx = df[df[0].str.contains(string)].index[0]
    return idx

def AtomInfoInput(a, atomInfo):
#    print(atomInfo)
    pos = [atomInfo[2], atomInfo[3], atomInfo[4]]
    a.setIndex(atomInfo[0])
    a.setAtomName(atomInfo[1])
    a.setPos(pos)
    a.setAtomType(atomInfo[5])
    a.setSubID(atomInfo[6])
    a.setSubName(atomInfo[7])
    a.setCharge(atomInfo[8])

def AtomInfoUpdate(atomsList_old, atomsList_new, idx=False, name=False, pos=False, atomType=False, subID=False, subName=False, charge=False):
    print('old_list length: {}'.format(len(atomsList_old)))
    print('new_list length: {}'.format(len(atomsList_new)))
    
    for i in range(len(atomsList_new)):
        if idx:
            pass
        if name:
            pass
        if pos:
            pos = atomsList_new[i].getPos()
            atomsList_old[i].setPos(pos)
        if atomType:
            pass
        if subID:
            pass
        if subName:
            pass
        if charge:
            atomsList_old[i].setCharge(atomsList_new[i].getCharge())
        
def BondInfoInput(b, bondInfo):
#    print(bondInfo)
    bond = bondInfo

    b.setIndex(bond[0])
    b.setAtom(bond[1], bond[2])
    b.setBondType(bond[3])

def MolInfoInput(index, start, end, name, atomsList):
    mol = Mol.MoleculeInfo()
    mol.setIndex(index)
    mol.setName(name)
    for i in range(start, end):
        mol.setAtoms(atomsList[i])
    return mol
    
def Atom2Mol(atomsList):
    moleculeList = []
    
    atoms = atomsList #List content is atoms
    iniID = atoms[0].getSubID()
    index = 0
    ini = 0
    for i in range(len(atoms)):
        subID = atoms[i].getSubID()
        if subID == iniID:
            if i == len(atomsList) - 1:
                index +=1
                name = "mol" +str(index)
                mol = MolInfoInput(index, ini, len(atomsList), name, atomsList)
                moleculeList.append(mol)
                return moleculeList
            else:
                continue
        else:
            index +=1
            name = "mol" + str(index)
            mol = MolInfoInput(index, ini, i, name, atomsList)
            moleculeList.append(mol)
            ini = i
            iniID = subID
    return moleculeList    
    
def ReadMol2(fileName):
    df = pd.read_csv(fileName, sep='\n', header=None)
    return df

def InfoInput(fileName, monLen, crosLen, monoNum, crosNum, dih=True):
    atomsList = []
    bondsList = []
    moleculesList = []
    monList = []
    crosList = []
    
    df = ReadMol2(fileName)
    sysName = df.iloc[GetIndex('MOLECULE', df)+1][0]
    sysInfo = df.iloc[2][0]
    content_1 = df.iloc[3][0]
    content_2 = df.iloc[4][0]
    basicPara = [sysName, sysInfo, content_1, content_2]
    atomStartIdx = GetIndex('ATOM', df) + 1
    bondStartIdx = GetIndex('BOND', df)
    atomDF = df.iloc[atomStartIdx:bondStartIdx][0].str.split().reset_index()
    bondDF = df.iloc[bondStartIdx+1:][0].str.split().reset_index()
    
    for i in range(len(atomDF)):
        a = Atom.AtomsInfo() #init a Atom structure
        atom = atomDF.iloc[i]#Atom information
        AtomInfoInput(a, atom[0])
        atomsList.append(a)
    
    for i in range(len(bondDF)):
        b = Bond.BondsInfo()
        bond = bondDF.iloc[i]
        BondInfoInput(b, bond[0])
        bondsList.append(b)
    
    if dih:
        moleculesList = Atom2Mol(atomsList)
        if len(moleculesList) < int(monoNum) + int(crosNum):
            atomsList = []
            for i in range(len(moleculesList)):    
                atomsNum = len(moleculesList[i].getAtoms())
                if atomsNum != monLen & atomsNum != crosLen:
                    mol = SplitMolecule(moleculesList[i], monLen, crosLen)
                    atomsList.append(mol.getAtoms())
                    print("after split")
                    moleculesList[i] = mol
                    moleculesList[i].outputInfo()
                else:
                    atomsList.append(moleculesList[i].getAtoms())
            atomsList = list(itertools.chain.from_iterable(atomsList))
            moleculesList = Atom2Mol(atomsList)    


        for i in range(len(moleculesList)):
            length = len(moleculesList[i].getAtoms())
            if length == monLen:
                monList.append(moleculesList[i])
            elif length == crosLen:
                crosList.append(moleculesList[i])
            else:
                print("Calculation of the molecule length met some error, please check!")
                print('monLen: ', monLen)
                print('crosLen: ', crosLen)
                print('length: ', length)
                sys.exit()         
        return atomsList, bondsList, moleculesList, monList, crosList, basicPara

def SplitMolecule(molecule, monLen, crosLen):
#    molecule.outputInfo()
    mol = Mol.MoleculeInfo()
    mol = copy.copy(molecule)
    mol.resetAtom()
    atoms = molecule.getAtoms()
    length = len(atoms)
    if length == 40:
        if atoms[0].getAtomName() == 'C': #Assume, monomer end is C, since it is epoxy, need to modify
            for i in range(0, monLen):
                id_old = atoms[i].getSubID()
                id_new = int(id_old) + 1
                atoms[i].setSubID(id_new)
                mol.setAtoms(atoms[i])
            for i in range(monLen, length):
                mol.setAtoms(atoms[i])
        elif atoms[0].getAtomName() == 'N':
            for i in range(0, crosLen):
                id_old = atoms[i].getSubID()
                id_new = id_old + 1
                atoms[i].setSubID(str(id_new))
                mol.setAtoms(atoms[i])
            for i in range(monLen, length):
                mol.setAtoms(atoms[i])
#    mol.outputInfo()
    return mol
#def main(filename, outputName, monLen, crosLen, monR, crosR, cutoff, bondsNum):
#monLen = 25
#crosLen = 15
#monNum = 4
#crosNum = 2
#a = InfoInput('tmp.mol2', monLen, crosLen, monNum, crosNum)
