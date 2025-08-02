#!/usr/bin/env python3
import os, argparse, subprocess, sys
from pathlib import Path
import pandas as pd
from loguru import logger
from dotenv import load_dotenv, set_key
# from utils.logger import setLog
from utils import Data, setLog
setupData = Data()

class Run:
    @logger.catch
    def __init__(self, parser: argparse.ArgumentParser):
        parser.add_argument(
            '-f', '--source_path', type=str, required=True, default=None, metavar="PATH"
        )
        parser.add_argument('-v', '--verbose', default="INFO", choices=[
            'TRACE', 'DEBUG', 'INFO', 'SUCCESS', 'WARNING', 'ERROR', 'CRITICAL'
        ], help='')
        parser.add_argument('-db', '--debug', action='store_true', default=True)
        parser.add_argument('-tm', '--track_mem', action='store_true', default=True)
        parser.add_argument('-rs', '--random_state', type=int, default=15)
        parser.add_argument('-vt', '--variance_threshold', type=float, default=0.1)
        args = parser.parse_known_args(['source_path', 'verbose', 'debug', 'track_mem', 'random_state', 'variance_threshold'])
        self.debug = str(args.debug).lower()=='true'
        self.track_mem = str(args.track_mem).lower()=='true'
        self.log = setLog(name='llm4pgcc.py', debug=self.debug, trackMem=self.track_mem
        self.source_path = str(args.source_path)
        self.verbose = str(args.verbose).upper()
        self.random_state = int(args.random_state); self.variance_threshold = float(args.variance_threshold)
        self.base_path=os.getcwd(); self.env_path=os.path.join(os.getcwd(), '.env')

        self.setupProject()
        if setupData.ConfigureRawData(sourcePath=self.source_path): self.rawData = setupData.getRawData()
        else: self.rawData = None

    # setup .env so far
    @logger.catch
    def setupProject(self):
        try:
            load_dotenv()
            set_key(self.env_path, "BASE_PATH", self.base_path); set_key(self.env_path, "ENV_PATH", self.env_path)
            set_key(self.env_path, "SOURCE_PATH", self.source_path)
            set_key(self.env_path, "VERBOSE", self.verbose); set_key(self.env_path, "DEBUG", self.debug); set_key(self.env_path, "TRACK_MEM", self.track_mem)
            set_key(self.env_path, "DEBUG", self.debug); set_key(self.env_path, "TRACK_MEM", self.track_mem)
            set_key(self.env_path, "RANDOM_STATE", str(self.random_state))
            set_key(self.env_path, "VARIANCE_THRESHOLD", str(self.variance_threshold))
            set_key(self.env_path, "PYTHONPYCACHEPREFIX", os.path.join(self.base_path, 'logs', '.cache'))
            return True
        except Exception as e: raise Exception(f'Error seting up .env: {e}')

    def run(self): pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='l4p', usage="llm4pgcc [-h, --help] [-v, --verbose] COMMAND",
        description='LLM4PGCC CLI TOOL'
    )
    run = Run(parser)