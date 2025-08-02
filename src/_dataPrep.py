import os, requests, json
import polars as pl
import numpy as np
import pubchempy as pcp
from loguru import logger
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, AllChem, MACCSkeys
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.model_selection import train_test_split
from utils import setLog, Data

featureColumns = {
    'rdkit': [desc[0] for desc in Descriptors.descList],
    'maccs': [f'MACCS_{i}' for i in range(167)],
    'ecfp4': [f'ECFP4_{i}' for i in range(2048)]
}
calc = MoleculeDescriptors.MolecularDescriptorCalculator(featureColumns['rdkit'])
uncharger = rdMolStandardize.Uncharger()
te = rdMolStandardize.TautomerEnumerator()
RDLogger.DisableLog('rdApp.*')

class DataPrep:
    def __init__(self, rawData: pl.DataFrame=None, rs:int=None, vt:float=None):
        self.base_path = os.getenv(
            "BASE_PATH",
            os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        )
        self.data_path = os.path.join(self.base_path, 'data')
        self.featPath = os.path.join(self.data_path, 'features'); os.makedirs(self.featPath, exist_ok=True)
        self.scafPath = os.path.join(self.data_path, 'scaffolds'); os.makedirs(self.scafPath, exist_ok=True)
        self.featMetaPath = os.path.join(self.featPath, 'meta')
        os.makedirs(self.featMetaPath, exist_ok=True)
        os.makedirs(os.path.join(self.scafPath, 'ecfp4'), exist_ok=True)
        os.makedirs(os.path.join(self.scafPath, 'maccs'), exist_ok=True)
        os.makedirs(os.path.join(self.scafPath, 'rdkit'), exist_ok=True)
        self.rs = rs if rs is not None else int(os.getenv("RANDOM_STATE", '15'))
        self.vt = vt if vt is not None else float(os.getenv("VARIANCE_THRESHOLD", '0.1'))
        self.logger = setLog(
            name="dataPrep",
            debug=os.getenv("DEBUG", 'false').lower()=='true',
            trackMem=os.getenv("TRACK_MEM", 'false').lower()=='true'
        )

        libDF, features = self.standardize(rawData=rawData)
        for key, val in features.items(): self.curateFeatures(df=pl.DataFrame(val), type=str(key))

    @logger.catch
    def normalize(self, nameKey, smi) -> tuple[bool, str]:
        log = self.logger
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

    @logger.catch
    def addFeatures(self, smi:str, nameKey:str) -> dict:
        log = self.logger
        try: mol = Chem.MolFromSmiles(smi)
        except Exception as e: log.debug(f'!-{nameKey}: feature extraction: could not configure Chem.Mol')
        if not mol:
            log.debug(f'!-{nameKey}: feature extraction: could not configure Chem.Mol')
            return None

        rdkitSer, maccsSer, ecfp4Ser = None, None, None
        try: rdkitSer = pl.Series(
            values=list(calc.CalcDescriptors(mol)),
            index=featureColumns['rdkit'],
            dtype=float
        )
        except Exception as e: log.debug(f'!-{nameKey}: feature extraction [rdkit]: {e}')
        try: maccsSer = pl.Series(
            values=list(MACCSkeys.GenMACCSKeys(mol)),
            index=featureColumns['maccs'],
            dtype=float
        )
        except Exception as e: log.debug(f'!-{nameKey}: feature extraction [maccs]: {e}')
        try: ecfp4Ser = pl.Series(
            values=list(AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nbits=2048)),
            index=featureColumns['ecfp4'],
            dtype=float
        )
        except Exception as e: log.debug(f'!-{nameKey}: feature extraction [ecfp4]: {e}')

        features = {'rdkit': rdkitSer, 'maccs': maccsSer, 'ecfp4': ecfp4Ser}
        return features

    def standardize(self, rawData:pl.DataFrame=None):
        log = self.logger
        parsed = {
            'nan': {
                'compound_name': [],
                'SMILES': []
            },
            'non-standardized': {},
            'duplicate': {},
            'invalid_score': {}
        }
        if rawData is not None: df = pl.DataFrame(rawData)
        else: df = pl.read_csv(os.path.join(self.data_path, 'raw_lib.csv'))
        if not {'compound_name', 'SMILES', 'score'}.issubset(df.columns): log.error("Missing required columns")
        else: df = df.cast({'compound_name': str, 'SMILES': str, 'score': float})
        lib = pl.DataFrame(schema=df.columns)
        log.info(f'standardizing raw compound library: {df.shape}')

        rdkit, ecfp4, maccs = [], [], []

        index = 0
        for row in df.iter_rows():
            index += 1
            name = str(row['compound_name']).strip()
            smi = str(row['SMILES']).strip()
            pScore = row['score']
            pathway = row['pathway'] if 'pathway' in df.columns else None
            target = row['target'] if 'target' in df.columns else None
            info = row['info'] if 'info' in df.columns else None

            # NaN and outlier removal:
            if name is None or name.lower() in ['', 'nan', 'none', 'na']: name = '(default)'
            if smi is None or smi.lower() in ['', 'nan', 'none', 'na']: # remove nan SMILES value
                parsed['nan']['SMILES'].append(f'[{index}]_{name}')
                log.debug(f'!-[{index}]_{name}: NaN removal (nan smiles)\n')
                continue
            if pScore is None or str(pScore).lower() in ['', 'nan', 'none', 'na']: # remove nan score value
                parsed['nan']['score'].append(f'[{index}]_{name}')
                log.debug(f'!-[{index}]_{name}: NaN removal (nan score)\n')
                continue
            if pScore<=0: # remove invalid score value (0 <= score <= 1)
                parsed['invalid_score'][f'[{index}]_{name}'] = pScore
                log.debug(f'!-[{index}]_{name}: invalid score: {pScore}\n')
                continue
            nameKey = f'[{index}]_{name}'

            # SMILES normalization:
            norm = self.normalize(smi, nameKey)
            if norm[0]: newsmi = norm[1].strip()
            else:
                parsed['non-standardized'][nameKey] = norm[1]
                continue

            # search for name on pubchempy if not listed
            if name=='(default)':
                url = f"https://cactus.nci.nih.gov/chemical/structure/{smi}/iupac_name"
                try:
                    response = requests.get(url)
                    response.raise_for_status()
                    name = response.text
                except Exception as e: log.warning(f'Pubchem lookup failed: {e}')

            # remove duplicates
            dupRow = lib.filter(pl.col("SMILES")==newsmi)
            if not dupRow.is_empty or dupRow is not None:
                oldName = str(dupRow['compound_name'].item()); oldSmi = str(dupRow['SMILES'].item()); oldScore = float(dupRow['score'].item())
                if pScore < oldScore: # drop previous iteration
                    lib = lib.remove(compound_name=oldName, SMILES=oldSmi, score=oldScore)
                    parsed['duplicate'][f'{oldName}_{oldSmi}'] = oldScore
                    log.debug(f'!-[{oldName}_{oldSmi}]: duplicate removal (previous)')
                elif pScore > oldScore:
                    parsed['duplicate'][f'{name}_{smi}'] = pScore
                    logger.debug(f'!-{nameKey}: duplicate removal (current)')
                    continue
                elif name==oldName and smi==oldSmi and pScore==oldScore:
                    logger.debug(f'!-{nameKey}: exact duplicates; skipping current')
                    continue
            smi = newsmi

            # add binary labels (and optional int labels)
            binLabel = 1 if pScore < 1 else 0

            row = {
                'compound_name': name,
                'SMILES': smi,
                'score': pScore,
                'label': binLabel
            }
            if 'pathway' in df.columns: row['pathway'] = pathway
            if 'target' in df.columns: row['target'] = target
            if 'info' in df.columns: row['info'] = info
            lib = pl.concat([lib, pl.DataFrame(row)])

            features = self.addFeatures(smi, nameKey)
            rdkit.append(pl.concat([row, pl.Series(features['rdkit'])]))
            maccs.append(pl.concat([row, pl.Series(features['maccs'])]))
            ecfp4.append(pl.concat([row, pl.Series(features['ecfp4'])]))

        log.info(f'completed data standardization and SMILES normalization.')
        lib.write_csv(os.path.join(self.data_path, 'lib.csv'))
        with open(os.path.join(self.base_path, 'logs', 'parsed_compounds.json'), 'w') as f: json.dump(parsed, f, indent=4)

        log.info(f'completed feature extraction.')
        rdkitDF = pl.DataFrame(rdkit)
        rdkitDF.write_csv(os.path.join(self.featMetaPath, 'rdkit_meta.csv'))
        maccsDF = pl.DataFrame(maccs)
        maccsDF.write_csv(os.path.join(self.featMetaPath, 'maccs_meta.csv'))
        ecfp4DF = pl.DataFrame(ecfp4)
        ecfp4DF.write_csv(os.path.join(self.featMetaPath, 'ecfp4_meta.csv'))

        return lib, {
            'rdkitFeatures': rdkitDF,
            'maccsFeatures': maccsDF,
            'ecfp4Features': ecfp4DF
        }

    # drop low-var features and split into test/train/valid sets
    def curateFeatures(self, df: pl.DataFrame, type:str):
        log = self.logger
        if type not in ['rdkit', 'maccs', 'ecfp4']: log.error("dataset type not supported: {type}")
        colToDrop = ['compound_name', 'SMILES', 'score', 'label', 'pathway', 'target', 'info']
        # infoDF = pl.DataFrame(df.filter([col for col in colToDrop in df.columns]))
        infoDF = pl.DataFrame(df.select([col for col in colToDrop if col in df.columns]))
        out = pl.DataFrame(df.drop([col for col in colToDrop if col in df.columns])).cast(int if type!='rdkit' else float)
        outnp = out.to_numpy()
        var = np.var(outnp, axis=0)
        lowVarCol = [out.columns[i] for i, variance in enumerate(var) if variance < self.vt]
        meta = pl.concat([infoDF, out], how='horizontal').drop(lowVarCol)
        meta.write_csv(os.path.join(self.featPath, f'{type}.csv'))

        currPath = os.path.join(self.scafPath, type); os.makedirs(currPath, exist_ok=True)

        def split(df: pl.DataFrame, name:str, path, rs:int=15):
            drop_cols = ['compound_name', 'score', 'pathway', 'target', 'info']
            train, temp = train_test_split(df, test_size=0.2, random_state=rs, stratify=df['label'], shuffle=True)
            valid, test = train_test_split(temp, test_size=0.5, random_state=rs, stratify=temp['label'], shuffle=True)
            pl.DataFrame(train).drop(columns=[col for col in drop_cols if col in train.columns]).to_csv(f'{path}/{name}_train.csv', index=False)
            pl.DataFrame(test).drop(columns=[col for col in drop_cols if col in test.columns]).to_csv(f'{path}/{name}_test.csv', index=False)
            pl.DataFrame(valid).drop(columns=[col for col in drop_cols if col in valid.columns]).to_csv(f'{path}/{name}_valid.csv', index=False)

        def splitSubcats(name, path, sc:dict=None, rs:int=15):
            for key, val in sc.items():
                if str(key)=='all': split(df=val, name=name, path=path)
                else:
                    subPath = os.path.join(path, key); os.makedirs(subPath, exist_ok=True)
                    split(df=val, name=key, path=subPath)

        if type=='rdkit':
            descs = pl.Series(list(out.columns))
            estateCols = [col for col in descs.to_list() if "EState" in str(col)]
            fgcCols = [col for col in estateCols if col.startswith('fr_')]
            mtCols = [col for col in fgcCols if str(col).lower() in ['balabanj', 'bertzct', 'hallkieralpha', 'ipc', 'avgipc']]
            fbCols = [col for col in mtCols if any(str(col).startswith(prefix) for prefix in ['Fp', 'BCUT2D'])]
            saCols = [col for col in fbCols if any(x in str(col).lower() for x in ['peoe', 'smr', 'slogp']) or str(col).lower()=='labuteasa']
            sdcCols = [
                col for col in saCols if (
                    str(col).lower().startswith('n') and str(col).lower() not in [
                        'numvalenceelectrons', 'numradicalelectrons'
                    ]
                ) or str(col).lower() in ['fractioncsp3', 'ringcount']
            ]
            pcCol = [col for col in descs.to_list() if col not in estateCols+fgcCols+mtCols+fbCols+saCols+sdcCols]
            subcategories = {
                'all': pl.DataFrame(meta),
                'EState': pl.concat([infoDF, pl.DataFrame(out).filter(estateCols)], how='horizontal'),
                'functionalGroupCounts': pl.concat([infoDF, pl.DataFrame(out).filter(fgcCols)], how='horizontal'),
                'molecularTopology': pl.concat([infoDF, pl.DataFrame(out).filter(mtCols)], how='horizontal'),
                'fingerprintBased': pl.concat([infoDF, pl.DataFrame(out).filter(fbCols)], how='horizontal'),
                'surfaceArea': pl.concat([infoDF, pl.DataFrame(out).filter(saCols)], how='horizontal'),
                'structural': pl.concat([infoDF, pl.DataFrame(out).filter(sdcCols)], how='horizontal'),
                'physiochemical': pl.concat([infoDF, pl.DataFrame(out).filter(pcCol)], how='horizontal'),
            }
            splitSubcats(name=type, path=currPath, sc=subcategories, rs=self.rs)
        else: split(df=meta, name=type, path=currPath, rs=self.rs)

if __name__ == "__main__":
    # initData: pl.DataFrame = setupSourceData()
    # cleanedDataDF = DataPrep(rawData=initData)
