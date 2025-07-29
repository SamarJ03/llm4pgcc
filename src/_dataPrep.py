import os, sys, requests
import polars as pl
import numpy as np
import pubchempy as pcp
from tabulate import tabulate
from loguru import logger

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, AllChem, MACCSkeys
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.model_selection import train_test_split

featureColumns = {
    'rdkit': [desc[0] for desc in Descriptors.descList],
    'maccs': [f'MACCS_{i}' for i in range(167)],
    'ecfp4': [f'ECFP4_{i}' for i in range(2048)]
}
calc = MoleculeDescriptors.MolecularDescriptorCalculator(featureColumns['rdkit'])
uncharger = rdMolStandardize.Uncharger()
te = rdMolStandardize.TautomerEnumerator()

RDLogger.DisableLog('rdApp.*')
from utils.logger import setLog

class DataPrep:
    def __init__(self, rawData_path:str=None):
        if rawData_path: self.rawData = rawData_path

    def standardize(self, db:bool=False, tm:bool=False):
        log = setLog(name="dataPrep.standardize", debug=db, trackMem=tm)
        parsed = {'nan': [], 'non-standardized': [], 'duplicate': [], 'invalid': []}
        parsed = {
            'nan': {
                'compound_name': [],
                'SMILES': [],
                'PGCC_score': []
            },
            'non-standardized': {},
            'duplicate': {},
            'invalid_score': {}
        }
        df = pl.DataFrame(self.rawData).cast({'compound_name': str, 'SMILES': str, 'PGCC_score': float})
        lib = pl.DataFrame(schema=df.columns); rows = []
        log.info(f'standardizing raw compound library: {df.shape}')

        index = 0
        for row in df.iter_rows():
            index += 1
            name = str(row['compound_name']).strip()
            smi = str(row['SMILES']).strip()
            pScore = row['PGCC_score']
            npScore = row['nonPGCC_score'] if 'nonPGCC_score' in df.columns else None
            pathway = row['pathway'] if 'pathway' in df.columns else None
            target = row['target'] if 'target' in df.columns else None
            info = row['info'] if 'info' in df.columns else None

            # NaN and outlier removal:
            if name is None or name.lower() in ['', 'nan', 'none', 'na']: name = '(default)'
            if smi is None or smi.lower() in ['', 'nan', 'none', 'na']: # remove nan SMILES value
                parsed['nan']['SMILES'].append(f'[{index}]_{name}')
                log.debug(f'!-[{index}]_{name}: NaN removal (nan smiles)\n')
                continue
            if pScore is None or str(pScore).lower() in ['', 'nan', 'none', 'na']: # remove nan PGCC_score value
                parsed['nan']['PGCC_score'].append(f'[{index}]_{name}')
                log.debug(f'!-[{index}]_{name}: NaN removal (nan PGCC_score)\n')
                continue
            if pScore<=0: # remove invalid PGCC_score value (0 <= score <= 1)
                parsed['invalid_score'][f'[{index}]_{name}'] = pScore
                log.debug(f'!-[{index}]_{name}: invalid PGCC_score: {pScore}\n')
                continue
            if npScore is not None or str(npScore).lower() in ['', 'nan', 'none', 'na'] or npScore<=0:
                npScore = pl.Null # reset for invalid nonPGCC_score
            nameKey = f'[{index}]_{name}'

            # SMILES normalization:
            norm = self.normalize(nameKey, smi, db=db, tm=tm)
            if norm[0]: newsmi = norm[1].strip()
            else: parsed['non-standardized'][nameKey] = norm[1]

            # search for name on pubchempy if not listed
            if name=='(default)':
                url = f"https://cactus.nci.nih.gov/chemical/structure/{smi}/iupac_name"
                response = requests.get(url)
                response.raise_for_status()
                name = response.text

            # remove duplicates
            dupRow = lib.filter(pl.col("SMILES")==newsmi)
            if not dupRow.is_empty or dupRow is not None:
                oldName = str(dupRow['compound_name'].item()); oldSmi = str(dupRow['SMILES'].item()); oldScore = float(dupRow['PGCC_score'].item())
                if pScore < oldScore: # drop previous iteration
                    lib = lib.remove(compound_name=oldName, SMILES=oldSmi, PGCC_score=oldScore)
                    parsed['duplicate'][f'{oldName}_{oldSmi}'] = oldScore
                    log.debug(f'!-[{oldName}_{oldSmi}]: duplicate removal (previous)')
                elif pScore > oldScore:
                    parsed['duplicate'][f'{name}_{smi}'] = pScore
                    logger.debug(f'!-{nameKey}: duplicate removal (current)')
                    continue
                elif name==oldName and smi==oldSmi and pScore==oldScore:
                    logger.debug(f'!-{nameKey}: exact duplicates; skipping current')
                    continue

            # add binary labels (and optional int labels)
            binLabel = 1 if pScore < 1 else 0
            intLabel = None
            if npScore is not None:
                npBinary = 1 if npScore < 1 else 0
                intLabel = binLabel + npBinary

            

    @logger.catch
    def normalize(self, nameKey, smi, db:bool=False, tm:bool=False) -> tuple[bool, str]:
        log = setLog(name='dataPrep.normalize', debug=db, trackMem=tm)
        log.debug(f'{nameKey}: Normalizing ')
        newsmi = None

        try: mol = Chem.MolFromSmiles(smi)
        except Exception as e:
            log.debug(f'!-{nameKey}: failed normalization (smi -> rdkit.Chem.Mol)')
            return False, str(e)

        try: cleanedMol = rdMolStandardize.Cleanup(mol)
        except Exception as e:
            log.debug(f'!-{nameKey}: failed normalization (H+ removal, metal atom disconnection, mol normalization+reionization)')
            return False, str(e)

        try: parentMol = rdMolStandardize.FragmentParent(cleanedMol)
        except Exception as e:
            log.debug(f'!-{nameKey}: failed normalization (retrieving parent fragment)')
            return False, str(e)

        try: unchargedMol = uncharger.uncharge(parentMol)
        except Exception as e:
            log.debug(f'!-{nameKey}: failed normalization (mol neutralization)')
            return False, str(e)

        try: tautMol = te.Canonicalize(unchargedMol)
        except Exception as e:
            log.debug(f'!-{nameKey}: failed normalization (tautomer canonicalization)')
            return False, str(e)

        try: newsmi = Chem.MolToSmiles(tautMol)
        except Exception as e:
            log.debug(f'!-{nameKey}: failed normalization (rdkit.Chem.Mol -> smi)')
            return False, str(e)

        return True, str(newsmi)


    # def getFeatures(self): pass
    # def evalFeatures(self): pass
    # def splitSets(self): pass

def setupSourceData(dataPath=os.path.join(os.getenv("BASE_PATH"), 'data')):
    log = setLog(name='dataPrep.setupSourceData')
    df = None
    for file in os.listdir(dataPath):
        name = str(os.fsdecode(file))
        if name=='raw_lib.csv': return pl.DataFrame(name)
        if name.startswith("source"):
            if name.endswith(".csv"): df = pl.read_csv(name, has_header=True)
            elif name.endswith(".xlsx"): df = pl.read_excel(name, has_header=True)
            else: log.warning(f'Incompatable types for source data: {name.split('.')[-1]}')
    if not isinstance(df, pl.DataFrame): log.error('No source data found in data/ directory..\n(Ensure that input file name begins with "source" and is either CSV or XLSX.)')

    rawData_df = pl.DataFrame(df).select(['compound_name', 'SMILES', 'PGCC_score'])
    toIncl = pl.DataFrame(df).filter(optCol for optCol in ['nonPGCC_score', 'pathway', 'target', 'info'] if optCol in df.columns)
    rawData_df = pl.concat([rawData_df, toIncl], how="horizontal")
    rawData_df.write_csv(os.path.join(dataPath, 'raw_lib.csv'))
    return rawData_df

if __name__ == "__main__":

    dataPath = os.path.join(os.getenv("BASE_DIR"), 'data')
    featPath = os.path.join(dataPath, 'features'); os.makedirs(featPath, exist_ok=True)
    scafPath = os.path.join(dataPath, 'scaffolds'); os.makedirs(scafPath, exist_ok=True)

    os.makedirs(os.path.join(featPath, 'raw'), exist_ok=True)
    os.makedirs(os.path.join(scafPath, 'ecfp4'), exist_ok=True)
    os.makedirs(os.path.join(scafPath, 'maccs'), exist_ok=True)
    os.makedirs(os.path.join(scafPath, 'rdkit'), exist_ok=True)

    rawData = setupSourceData(dataPath=dataPath)

    cleanedDataDF = DataPrep()
