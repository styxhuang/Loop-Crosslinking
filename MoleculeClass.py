# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 21:28:26 2018

@author: HuangMing
"""

class MoleculeInfo(object):
    def __init__(self):
        self.__index__ = 0
        self.__name__ = ""
        self.__atoms__ = []
        self.__reactAtoms__ = []
    
    def resetAtom(self):
        self.__atoms__ = []
        
    def setIndex(self, index):
        self.__index__ = index
    
    def setName(self, name):
        self.__name__ = name
        
    def setAtoms(self, atom):
        self.__atoms__.append(atom)
    
    def setReact(self, index):
        self.__reactAtoms__.append(index)
        
    def getIndex(self):
        return self.__index__
    
    def getName(self):
        return self.__name__
    
    def getAtoms(self):
        return self.__atoms__
    
    def getReact(self):
        return self.__reactAtoms__
    
    def outputInfo(self, index=1):
#        if (len(self.__reactAtoms__) == 0):
#            print("Index: ", self.__index__, "\tName: ", self.__name__, "\t", len(self.__atoms__), "\t")
#        else:
#            print("Index: ", self.__index__, "\tName: ", self.__name__, "\t", len(self.__atoms__), "\t", "Reactive atoms: ", self.__reactAtoms__[0])
#            if index == 1:
        for i in range (len(self.__atoms__)):
            index = self.__atoms__[i].getIndex()
            name = self.__atoms__[i].getAtomName()
            subName = self.__atoms__[i].getSubName()
            subID = self.__atoms__[i].getSubID()
            print(index, "\t", name, "\t", subName, "\t", subID)
                
            
    
    