# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 11:15:53 2018

@author: huang
"""

import pandas as pd

class AtomsInfo(object):
    def __init__(self):
        self.__Index__ = 0.
        self.__AtomName__ = ""
        self.__X__ = 0.
        self.__Y__ = 0.
        self.__Z__ = 0.
        self.__AtomType__ = ""
        self.__sub_id__ = 0
        self.__sub_name__ = ""
        self.__charge__ = 0.
    def setIndex(self, index):
        self.__Index__ = index
        
    def setAtomName(self, name):
        self.__AtomName__ = name
    
    def setPos(self, pos):
        self.__X__ = pos[0]
        self.__Y__ = pos[1]
        self.__Z__ = pos[2]
    
    def setAtomType(self, Type):
        self.__AtomType__ = Type
    
    def setSubID(self, id):
        self.__sub_id__ = id      
  
    def setSubName(self, subName):
        self.__sub_name__ = subName
    
    def setCharge(self, charge):
        self.__charge__ = charge
    
    def getIndex(self):
        return self.__Index__
    
    def getAtomName(self):
        return self.__AtomName__
                
    def getPos(self):
        return [self.__X__, self.__Y__, self.__Z__]
    
    def getAtomType(self):
        return self.__AtomType__
    
    def getSubID(self):
        return self.__sub_id__
        
    def getSubName(self):
        return self.__sub_name__
    
    def getCharge(self):
        return self.__charge__
    
    def outputData(self):
        str1 = '{:>7}{:>3}\t\t{:>9}{:>9}{:>9}\t{:<6}\t{:>3}\t{:>6}{:>10}\n'.format(
                str(self.__Index__), self.__AtomName__, self.__X__, self.__Y__, self.__Z__, self.__AtomType__, 
                str(int(self.__sub_id__)), self.__sub_name__, self.__charge__)
        return str1
        
        
        
    