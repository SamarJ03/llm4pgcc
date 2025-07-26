import os, argparse
import pandas as pd
from loguru import logger
from dotenv import load_dotenv
load_dotenv()

class VerifyData:
    def __init__(self, initialParser: argparse.ArgumentParser):
        dataParser = initialParser.add_argument_group(title='input')
        dataParser.add_argument('-sl', '--set_labels', type=dict[str,str], default=None, help="")
        args, _ = dataParser.parse_known_args(['set_labels'])
        self.set_labels = args.set_labels
        self.sourceData_path = None

        dataPath = os.getenv("DATA_PATH")
        if not os.exists(dataPath): logger.error(f'Incorrect path to source data:\n{dataPath}')

        df = None
        for file in os.listdir(dataPath):
            name = str(os.fsdecode(file))
            if not name.endswith(['.csv', '.xlsx']) and name.startswith('source'): continue
            path = os.path.join(os.getenv("DATA_PATH"), name)
            df = pd.read_csv(path, nrows=0) if name.endswith('.csv') else pd.read_excel(path, nrows=0)
        if not df: logger.error('No source data found in data/ directory..\n(Ensure that input file name begins with "source" and is either CSV or XLSX.)')

        if isinstance(self.set_labels, dict[str,str]): df = self.changeLabels(df)

        rawData3 = self.sortData(df)
        rawData3.to_csv(os.path.join(os.getenv("DATA_PATH"), 'raw_lib.csv'))


    def sortData(self, input: pd.DataFrame):
        toInclude = ['compound_name'] if 'compound_name' in input.columns else []
        toInclude.append('SMILES', 'PGCC_score')
        toInclude.append(optCol for optCol in ['nonPGCC_score', 'pathway', 'target', 'info'] if optCol in input.columns)
        return pd.DataFrame(input, columns=toInclude, copy=True)

    def changeLabels(self, input: pd.DataFrame):
        df = pd.DataFrame(input, copy=True)
        new_labels = {obs: exp for exp, obs in dict(self.set_labels).items()}
        if 'SMILES' not in new_labels.keys() and 'SMILES' not in input.columns: logger.error('"SMILES" field not found in source data.')
        if 'PGCC_score' not in new_labels.keys() and 'PGCC_score' not in input.columns: logger.error('"PGCC_score" field not found in source data.')

        for obs, exp in new_labels.items():
            if obs not in df.columns:
                logger.info(f'{obs} not found in source data fields. Field rename skipped..')
                continue
            df.rename({obs:exp}, axis=1)
        return df