# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 10:41:31 2018

@author: huang
"""

class BondsInfo(object):
    def __init__(self):
        self.__Index__ = 0
        self.__Atom1__ = 0
        self.__Atom2__ = 0
        self.__BondType = ""
    
    def setIndex(self, index):
        self.__Index__ = index
        
    def setAtom(self, atomIndex1, atomIndex2):
        self.__Atom1__ = atomIndex1
        self.__Atom2__ = atomIndex2
        
    def setBondType(self, Type):
        self.__BondType__ = Type
        
    def getIndex(self):
        return self.__Index__
    
    def getAtom(self):
        return self.__Atom1__, self.__Atom2__
    
    def getBondType(self):
        return self.__BondType__
    
    def outputData(self):
        str1 = '{:>7}{:>6}{:>6}{:>5}\n'.format(str(self.__Index__), str(self.__Atom1__), str(self.__Atom2__), str(self.__BondType__))
        return str1
        