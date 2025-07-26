import os, sys, json
import pandas as pd
import numpy as np
from dotenv import load_dotenv, set_key
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, AllChem, MACCSkeys
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.model_selection import train_test_split
from tabulate import tabulate
from loguru import logger
RDLogger.DisableLog('rdApp.*')

# logger.add(os.getenv('CONSOLE_LOG_PATH)', fileName=__name__, )

class Clean:
    def __init__(self): pass

    def removeNaN(self): pass
    def standardizeSmiles(self): pass
    def removeDuplicates(self): pass
    def removeOutliers(self): pass
    def addBinLabel(self): pass

class Features:
    def __init__(self): pass

    def getFeatures(self): pass
    def evalFeatures(self): pass
    def splitSets(self): pass

if __name__ == "__main__": pass