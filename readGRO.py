# -*- coding: utf-8 -*-
"""
Created on Sat Mar 17 14:34:14 2018

@author: HuangMing
"""

import pandas as pd
import AtomClass as Atom
import BondClass as Bond
import crosslinking as MOL2
def CheckIndex(df, keyword):
    index = df.index[df[0].str.contains(keyword) == True].tolist()
    return index

def ChkAtomType(atomName, atomType):
    if atomName == 'C':
        if atomType == 'opls_145' or atomType == "opls_166":
            return 'C.ar'
        else:
            return 'C.3'
    elif atomName == 'N':
        return 'N.3'
    elif atomName == 'O':
        return 'O.3'
    elif atomName == 'H':
        return 'H'

def RemoveRows(df, keyword):
    tmp = df[df[0].str.contains(keyword) == False].reset_index(drop=True)
    tmp = tmp.iloc[:len(df)-1].reset_index(drop=True)
    return tmp

def InitData(GRO, TOP):
   
    atomStartIdx = CheckIndex(TOP, 'atoms')
    atomEndIdx = CheckIndex(TOP, '; total molecule charge')
    
    bondStartIdx = CheckIndex(TOP, 'bonds')
    bondEndIdx = CheckIndex(TOP, 'pairs') #Since I have delete the constraints part, need to modify here
    
    atomsTop = TOP.iloc[atomStartIdx[0]+2:atomEndIdx[0]]
    atomsGro = GRO.iloc[2:].iloc[:-1]
    
    bondsTop = TOP.iloc[bondStartIdx[0]+2:bondEndIdx[0]].reset_index(drop=True)
#Remove the hydrogen needs to update the atom index in the bond session, may be complexed    
#    atomsTopNoH = RemoveRows(atomsTop, " H ")
#    atomsGroNoH = RemoveRows(GRO.iloc[2:], "H ").iloc[:-1]
 
    return atomsGro, atomsTop, bondsTop

def AtomInfoInput(groData, topData):
    atoms = []
    arIndex = []
    for i in range(len(groData)):
        oriGRO = groData.iloc[i].str.split().tolist()[0] #from pandas datafrom to list
        oriTOP = topData.iloc[i].str.split().tolist()[0]
        if len(oriGRO) == 8:
#            print('tst-1: ', oriGRO[1])
            tmp1 = [oriGRO[1][0], oriGRO[1][1:]]
            oriGRO[1]=tmp1[0]
            oriGRO.insert(2, tmp1[1])
        #From gro file input
        atom = Atom.AtomsInfo()
        atom.setIndex(oriGRO[2])
        atom.setAtomName(oriGRO[1])
        atom.setPos([str(round(float(oriGRO[3])*10,2)),str(round(float(oriGRO[4])*10,2)),str(round(float(oriGRO[5])*10,2))])
        atom.setSubName(oriGRO[0])
        
        #From top file input
        atomType = ChkAtomType(oriGRO[1], oriTOP[1])
        atom.setAtomType(atomType)
        atom.setSubID(oriTOP[2])
        atom.setCharge(oriTOP[6])        
        if atomType == 'C.ar':
            arIndex.append(atom.getIndex())
        atoms.append(atom)
    
#    tmp = atoms[0].outputData()
#    print('tst-1: ', tmp)
#    for i in range(len(atoms)):
#        print('index: ', i, '\ttst-1: ', atoms[i].outputData())

    return atoms, arIndex

def BondInfoInput(topData, arIndex):
    bonds=[]
    for i in range(len(topData)):
        oriTOP = topData.iloc[i].str.split().tolist()[0]
        bond = Bond.BondsInfo()
        bond.setIndex(i+1)
        bond.setAtom(oriTOP[0], oriTOP[1])
        atomTmp = [oriTOP[0], oriTOP[1]] 
        if any(i in atomTmp for i in arIndex):
            #print("right!")
            bond.setBondType("ar")
        else:
            bond.setBondType(oriTOP[2])
        bonds.append(bond)
    return bonds

def ExportMol2(TOP, atomsNum, bondsNum, atoms, bondList, cycle, fileName):
    atomList = atoms[0]
    name = TOP.iloc[-1][0].split()[0]
    outputName = '{}-stp-{}.mol2'.format(fileName, cycle)
    content = [' {:} {:} {:} {:} {:}'.format(atomsNum, bondsNum, 0, 0, 0), 'SMALL', 'GASTEIGER']
    MOL2.ExportMOL2(name, outputName, content, atomList, bondList)
    
###############################################################################
def Main(fileName, topName, cycle):
    groName = '{}.gro'.format(fileName)
    topName = topName
    
    GRO = pd.read_csv(groName, sep="\n", header=None)
    TOP = pd.read_csv(topName, sep="\n", header=None)
    
    a = InitData(GRO, TOP)
    atoms = AtomInfoInput(a[0],a[1])
    bonds = BondInfoInput(a[2], atoms[1])
    atomsNum = len(atoms[0])
    bondsNum = len(bonds)
    ExportMol2(TOP, atomsNum, bondsNum, atoms, bonds, cycle, fileName)
    
