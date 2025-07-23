import os, sys, argparse
from loguru import logger
from dotenv import load_dotenv
load_dotenv()

class PrepareData:
    def __init__(self, initialParser: argparse.ArgumentParser):
        dataParser = initialParser.add_argument_group(title='dataset')

        dataParser.add_argument('-p', '--data_path', type=str, default="", help='')
        dataParser.add_argument('-l', '--set_labels', type=dict[str:str], required=False, help="")

        dataType = dataParser.add_mutually_exclusive_group(required=True)
        dataType.add_argument('-f', '--full', action='store_true')
        dataType.add_argument('-s', '--split', action='store_true')
        dataType.add_argument('-r', '--ready', action='store_true')

        args, _ = dataParser.parse_known_args(['full', 'split', 'ready', 'set_labels', 'data_path'])
        if args.full: self.data_mode = 'full'
        elif args.split: self.data_mode = 'split'
        elif args.ready: self.data_mode = 'ready'
        if args.set_labels is not None: self.set_labels = args.set_labels
        self.data_path = args.data_path

    def verifyData():
        dataPath = str(tempArgs['data_path'])
        configPath = os.path.join(os.getcwd(), 'config')

        if not os.path.exists(configPath):
            logger.error('config/ directory path not found..')
            return False, None

        for file in os.listdir(configPath):
            fileName = str(os.fsdecode(file))
            if dataPath=="": pass
                # if fileName.endswith(['.csv', '.xlsx']):
                #     filePath = os.path.join(configPath, fileName)
                #     import pandas as pd
                #     dataCol = pd.read_csv(filePath, nrows=0).columns if fileName.endswith('.csv') else pd.read_excel(filePath, nrows=0).columns
                # return verifyFileType()
            elif dataPath.endswith(['.csv', '.xlsx']):
                if fileName==dataPath: return True, os.path.join(configPath, fileName)
                else:
                    logger.error(f'No files match {dataPath} within config/')
                    return False, None
            else:
                logger.error(f'Invalid --data_path entry:\n{dataPath}')
                return False, None

    def verifyFileType(fileName: str, configPath: str):
        filePath = os.path.join(configPath, fileName)
        dataCol = pd.read_csv(filePath, nrows=0).columns if fileName.endswith('.csv') else pd.read_excel(filePath, nrows=0).columns