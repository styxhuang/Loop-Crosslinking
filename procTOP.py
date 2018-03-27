# -*- coding: utf-8 -*-
"""
v0.0: Created on Sat Mar 17 17:44:39 2018
After create the bond, since the bond length sometimes is longer than usual, and gromacs cannot recognize.
After we transfer mol2 file to gro file, some bond/angle/dihedral coefficient may miss. This file helps
us to fill in those blank.

@author: HuangMing
"""

import pandas as pd

#Following are Gaff data
dictBond = {
        'C-N': '1, 0.14700, 268278.'
        }

dictAngle = { #TODO: keep update angle coefficient
        'N-C-N': '1, 110.380, 553.9616',
        'C-O-C': '1, 117.600, 522.1632',
        'C-C-N': '1, 110.380, 553.9616',
        'C-N-C': '1, 110.900, 535.5520',
        'N-C-C': '1, 110.380, 553.9616',
        'H-C-N': '1, 109.920, 413.3792',
        'C-N-H': '1, 109.920, 394.1328'
        }

dictDihedral = { #Has little different on the 3 column, some combination is 2 and 3'
        'C-C-O-C': '1, 180.000, 3.766, 2',
        'C-O-C-C': '1, 180.000, 4.602, 2',
        'C-N-C-C': '1, 180.000, 2.008, 2',
        'C-C-N-C': '1, 180.000, 2.008, 2'
        }

def CheckIndex(df, keyword):
    index = df.index[df[0].str.contains(keyword) == True].tolist()
    return index

def ChkAnglePair(anglesIndex, TOP):
    string = TOP.iloc[anglesIndex].str.split()[0]
    pairs = [string[-3], string[-2], string[-1]]
    return pairs

def ChkErRow(df):
    index = []
    criNum = df.iloc[0][0].count(' ')
    
    for i in range(len(df)):
        spaceNum = df.iloc[i][0].count(' ')
        if spaceNum > criNum + 10: #If doesn't update the data right, modify this number
#            print(df.iloc[i]['index'])
            index.append(df.iloc[i]['index'])
#            index.append(df.iloc[i].index.values[0])
    return index
def ProcessString(index, df, category='bond'):
    a = '1'             
    b = '0.14700'
    c = '268278.'
    for i in range(len(index)):
        tmp = df.iloc[index[i]].str.split()[0]
        if category == 'bond':
            str1 = "{:>8}{:>6}{:>4}{:>12}{:>13}{:>8}{:>7}{:>6}".format(tmp[0],tmp[1],a , b, c, tmp[2], tmp[3],tmp[4])
            print(str1)
            df.iloc[index[i]] = str1
        elif category == 'angle':
            tmp1 = tmp[-3:]
            tmp1 = ''.join(tmp1)
            if tmp1 in dictAngle:
                coeff = dictAngle[tmp1].split(',')
            str1 = "{:>6}{:>6}{:>6}{:>4}{:>12}{:>12}{:>6}{:>7}{:>7}{:>6}".format(tmp[0],tmp[1],tmp[2], coeff[0], coeff[1], coeff[2], 
                    tmp[3], tmp[-3], tmp[-2],tmp[-1])
            print(str1)
            df.iloc[index[i]] = str1
        elif category == 'dihedral':
            print(index[i])
            tmp2 = tmp[-4:]
            tmp2 = ''.join(tmp2)
            if tmp2 in dictDihedral:
                coeff = dictDihedral[tmp2].split(',')
                str2 = "{:>6}{:>6}{:>6}{:>6}{:>4}{:>12}{:>12}{:>10}{:>3}{:>4}{:>8}{:>7}{:>7}{:>6}".format(tmp[0],tmp[1],tmp[2], tmp[3], coeff[0], coeff[1], coeff[2], coeff[3],
                    tmp[-6], tmp[-5], tmp[-4], tmp[-3], tmp[-2],tmp[-1])
                print(str2)
                df.iloc[index[i]] = str2
            else:
                df = df.drop(index[i])
        else:
            print("Unknow data type")

def UpdateRow(index, df, category='bond'): #For now only bond section needs to be updated, or the coefficient is for bond
    if category == 'bond':
        ProcessString(index, df, 'bond')
    elif category == 'angle':
        ProcessString(index, df, 'angle')
    elif category == 'dihedral':
        ProcessString(index, df, 'dihedral')
    return df

def DropRows(index, df): #For the empty dihedral, for now I just delete them
    tmp = df
    for i in range(len(index)):
        tmp = tmp.drop([index[i]])
    return tmp
def BondProc(TOP):
    bondStartIdx = CheckIndex(TOP, 'bonds')
    if CheckIndex(TOP, 'constraints') == []:
        bondEndIdx = CheckIndex(TOP, 'pairs')
    else:
        bondEndIdx = CheckIndex(TOP, 'constraints')
    bondsTOP = TOP.iloc[bondStartIdx[0]+2: bondEndIdx[0]].reset_index()
    bondIdx = ChkErRow(bondsTOP)
    df = UpdateRow(bondIdx, TOP)
    
    return df

def AngleProc(TOP):
    angleStartIdx = CheckIndex(TOP, 'angles')
    angleEndIdx = CheckIndex(TOP, 'dihedrals')
    
    anglesTOP = TOP.iloc[angleStartIdx[0]+2:angleEndIdx[0]].reset_index()
    anglesIndex = ChkErRow(anglesTOP)
    finalDF = UpdateRow(anglesIndex, TOP, 'angle')
    return finalDF

def DihedralProc(df):
    dihIdx = CheckIndex(df, 'dihedral')
    if len(dihIdx) >= 2:
        print('Tst start')
        index = list(range(dihIdx[0]+2, dihIdx[1]))
        print(index)      
        df = UpdateRow(index, df, 'dihedral')
        dihTOP = df.iloc[dihIdx[0]+2:dihIdx[1]].reset_index()
        print(dihTOP)
    else:
        dihStart = dihIdx[0]
        dihEnd = CheckIndex(df, '; Include Position restraint file')
        dihTOP = df.iloc[dihStart+2:dihEnd[0]].reset_index()
    
    idx = ChkErRow(dihTOP)
    finalDF = DropRows(idx, df)
    return finalDF

def DelEmptyRow(df): #Doesn't work right, need to check again
    idx = ChkErRow(df)
    finalDF = DropRows(idx, df)
    return finalDF

def ConstraintsProc(df):
    constrainStartIdx = CheckIndex(df, 'constraints')
    if constrainStartIdx == []:
        return df
    else:
        constrainEndIdx = CheckIndex(df, 'pairs')
        for i in range(constrainStartIdx[0], constrainEndIdx[0]):
            df = df.drop([i])
        return df    
    

def ExportTOP(filename, df):
    f = open(filename, 'w')
    for i in range(len(df)):
        line = df.iloc[i][0] + '\n'
        if '[ atoms ]' in line:
            line = '\n' + line
        f.write(line)
    f.close()

def Main(topName):
#topName = 'topol.top'
    TOP = pd.read_csv(topName, sep="\n", header=None)
    bondProc = BondProc(TOP)
    angleProc = AngleProc(bondProc)
    dihedralProc = DihedralProc(angleProc)
    #delEmptyRow = DelEmptyRow(dihedralProc)
    finalProc = ConstraintsProc(dihedralProc)
    ExportTOP('topol.top-end', finalProc)
