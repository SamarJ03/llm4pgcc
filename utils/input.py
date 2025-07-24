import os, sys, argparse
from loguru import logger
from dotenv import load_dotenv
load_dotenv()

class VerifyData:
    def __init__(self, initialParser: argparse.ArgumentParser):
        dataParser = initialParser.add_argument_group(title='input')

        dataParser.add_argument('-sl', '--set_labels', type=dict[str:str], default=None, help="")
        dataParser.add_argument('-vd', '--verify_data', action='store_true', help='')

        args, _ = dataParser.parse_known_args(['verify_data', 'set_labels'])
        self.verify_path = args.verify_path
        self.verify_data = args.verify_data
        self.set_labels = args.set_labels if isinstance(args.set_labels, dict[str:str]) else False
        self.data_path = args.data_path

    # PATH LAYOUT PLAN FOR DATA/:

    # source_lib.xlsx (user input csv or xlsx file)
    # raw_lib.csv (user input in CSV format, only required params kept)
    # clean_lib.csv (raw_lib.csv after data preprocessing steps)
    # lib.json (JSON file of clean_lib.csv compounds with BOTH required and optional fields)
    # features.csv (CSV of all compounds from clean_lib.csv and maccs/ecfp4/rdkit descriptors)
    # features/
    #  - ecfp4_meta.csv (ecfp4 descriptors prior to low-var removal)
    #  - maccs_meta.csv (maccs descriptors prior to low-var removal)
    #  - rdkit_meta.csv (rdkit descriptors prior to low-var removal)
    #  - ridLowVar/ (low var descriptors removed)
    #    - ecfp4.csv
    #    - maccs.csv
    #    - rdkit.csv
    # scaffold/
    #  - ecfp4/ (contains train/test/valid sets)
    #  - maccs/ (contains train/test/valid sets)
    #  - rdkit/ (contains subgroupings train/test/valid sets including file of all rdkit descriptors train/test/valid set)


    # def verifyData(self):
    #     dataPath = str(tempArgs['data_path'])
    #     configPath = os.path.join(os.getcwd(), 'config')

    #     if not os.path.exists(configPath):
    #         logger.error('config/ directory path not found..')
    #         return False, None

    #     for file in os.listdir(configPath):
    #         fileName = str(os.fsdecode(file))
    #         if dataPath=="": pass
    #             # if fileName.endswith(['.csv', '.xlsx']):
    #             #     filePath = os.path.join(configPath, fileName)
    #             #     import pandas as pd
    #             #     dataCol = pd.read_csv(filePath, nrows=0).columns if fileName.endswith('.csv') else pd.read_excel(filePath, nrows=0).columns
    #             # return verifyFileType()
    #         elif dataPath.endswith(['.csv', '.xlsx']):
    #             if fileName==dataPath: return True, os.path.join(configPath, fileName)
    #             else:
    #                 logger.error(f'No files match {dataPath} within config/')
    #                 return False, None
    #         else:
    #             logger.error(f'Invalid --data_path entry:\n{dataPath}')
    #             return False, None

    # def verifyFileType(fileName: str, configPath: str):
    #     filePath = os.path.join(configPath, fileName)
    #     dataCol = pd.read_csv(filePath, nrows=0).columns if fileName.endswith('.csv') else pd.read_excel(filePath, nrows=0).columns