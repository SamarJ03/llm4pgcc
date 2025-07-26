#!/usr/bin/env python3
import os, argparse, subprocess, sys
from pathlib import Path
import pandas as pd
from loguru import logger
from dotenv import load_dotenv, set_key
from utils import log

class Setup:
    def __init__(self):
        parser = argparse.ArgumentParser(
            prog='l4p', usage="llm4pgcc [-h, --help] [-v, --verbose] COMMAND",
            description='LLM4PGCC CLI TOOL'
        )
        parser.add_argument('-v', '--verbose',
            default="INFO", choices=[
                'TRACE', 'DEBUG', 'INFO', 'SUCCESS',
                'WARNING', 'ERROR', 'CRITICAL'
            ], help=''
        )
        vrb = str(parser.parse_known_args(['verbose']).verbose)

        # configure .env
        self.base_path = os.getcwd(); set_key(self.env_path, "BASE_PATH", self.base_path, export=True)
        self.env_path = os.path.join(self.base_path, '.env'); set_key(self.env_path, "ENV_PATH", self.env_path, export=True)
        if not os.path.exists(self.env_path): subprocess('touch $(pwd)/.env')
        if not load_dotenv(): raise Exception('Unable to load .env file.')

        self.verbose = vrb; set_key(self.env_path, "VERBOSE", vrb, export=True)
        for dirName in ['data', 'logs', 'resources', 'results', 'src', 'utils']:
            path = os.path.join(self.base_path, dirName); os.makedirs(path, exist_ok=True)
            setattr(self, f'{dirName}_path', path); set_key(self.env_path, f'{dirName.upper()}_PATH', path, export=True)

        log.Log(path=os.getcwd())
        self.setupSourceData()

        return str(self.verbose)

    def setupSourceData(self):
        # configure source input data
        df = None
        for file in os.listdir(self.data_path):
            name = str(os.fsdecode(file))
            if not name.endswith('.csv', '.xlsx') and name.startswith('source'): continue
            path = os.path.join(os.getenv("DATA_PATH"), name)
            df = pd.read_csv(path, nrows=0) if name.endswith('.csv') else pd.read_excel(path, nrows=0)
        if not df: logger.error('No source data found in data/ directory..\n(Ensure that input file name begins with "source" and is either CSV or XLSX.)')

        toInclude = ['compound_name'] if 'compound_name' in df.columns else []
        toInclude.append('SMILES', 'PGCC_score')
        toInclude.append(optCol for optCol in ['nonPGCC_score', 'pathway', 'target', 'info'] if optCol in df.columns)
        rawData = pd.DataFrame(df, columns=toInclude, copy=True)
        rawData.to_csv(os.path.join(self.data_path, 'raw_lib.csv'))

    # def setupLogger(self, vrb:str):
    #     # setup loguru.logger
    #     logger.remove()
    #     if vrb.lower() not in ['debug', 'info', 'warning', 'error', 'critical']:
    #         logger.warning(f'"{vrb}" is not a valid verbose level. Default set to "INFO".'); vrb = "INFO"
    #     logger.add(
    #         sys.stderr,
    #         format="\n<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
    #         colorize=True, level=vrb
    #     )
    #     appLog_path = os.path.join(self.logs_path, 'application_logs'); os.makedirs(appLog_path, exist_ok=True)
    #     for subDir in ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']:
    #         path = os.path.join(appLog_path, f'{subDir.lower()}_logs'); os.makedirs(path, exist_ok=True)
    #         logger.add(
    #             os.path.join(path, f'application_{subDir}_{{time}}.log'),
    #             level=subDir, rotation='10 MB'
    #         )
    #     perfLog_path = os.path.join(self.logs_path, 'performace_logs'); os.makedirs(perfLog_path, exist_ok=True)
    #     for subDir in ['dataPrep', 'createPrompts', 'synthesize', 'inference', 'summarize', 'codeGen']:
    #         path = os.path.join(perfLog_path, f'{subDir}_logs'); os.makedirs(path, exist_ok=True)
    #         logger.add(
    #             os.path.join(path, f'performance_{subDir}_{{time}}.log'),
    #             level="", rotation='10 MB'
    #         )

class Run: pass

if __name__ == "__main__":
    vrb = Setup()

