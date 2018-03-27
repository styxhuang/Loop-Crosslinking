# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 21:45:28 2018

@author: HuangMing
"""

class PolymerInfo(object):
    def __init__(self):
        self.__name__ = ""
        self.__molecule = []
        
    def setName(self, name):
        self.__name__ = name
        
    def setMol(self, molList):
        self.__molecule__.append(molList)
        
    def getName(self):
        return self.__name__
    
    def getMol(self):
        return self.__molecule__
    
    