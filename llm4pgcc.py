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
        parser.add_argument('-db', '--debug', action='store_true', default=True)
        parser.add_argument('-tm', '--track_mem', action='store_true', default=True)
        vrb = str(parser.parse_known_args(['verbose']).verbose)
        args = parser.parse_known_args(['verbose', 'debug', 'track_mem'])
        self.verbose = str(args.verbose).upper(); self.debug = str(args.debug).lower()=='true'; self.track_mem = str(args.track_mem).lower()=='true'

        load_dotenv()
        self.env_path = os.path.join(self.base_path, '.env'); set_key(self.env_path, "ENV_PATH", self.env_path)
        self.base_path = os.getcwd(); set_key(self.env_path, "BASE_PATH", self.base_path)
        set_key(self.env_path, "VERBOSE", vrb); set_key(self.env_path, "DEBUG", self.debug); set_key(self.env_path, "TRACK_MEM", self.track_mem)
        set_key(self.env_path, "PYTHONPYCACHEPREFIX", os.path.join(self.base_path, 'logs', '.cache'))
        return self.verbose

    def setEnv(self): pass


class Run: pass

if __name__ == "__main__":
    vrb = Setup()

