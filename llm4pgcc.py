#!/usr/bin/env python3
import os, argparse, subprocess, sys
from pathlib import Path
import pandas as pd
from loguru import logger
from dotenv import load_dotenv, set_key
from utils.logger import setLog
log = setLog(name=__name__, debug=False, trackMem=False)

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
        load_dotenv()
        self.env_path = os.path.join(os.getcwd(), '.env')
        self.base_path = os.getcwd(); set_key(self.env_path, "BASE_PATH", self.base_path)
        self.verbose = vrb; set_key(self.env_path, "VERBOSE", vrb)
        set_key(self.env_path, "PYTHONPYCACHEPREFIX", os.path.join(self.base_path, 'logs', '.cache'))

        # self.setupSourceData()

        return str(self.verbose)

class Run: pass

if __name__ == "__main__":
    vrb = Setup()

