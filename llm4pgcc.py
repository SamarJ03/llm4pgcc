#!/usr/bin/env python3
import os, argparse, subprocess, sys
from pathlib import Path
import pandas as pd
from loguru import logger
from dotenv import load_dotenv, set_key, get_key
logger.remove()
logger.add(sys.stderr, level="INFO")

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

        # set class and .env path variables
        self.base_path = os.getcwd(); set_key(self.env_path, "BASE_PATH", self.base_path, export=True)
        self.env_path = os.path.join(self.base_path, '.env'); set_key(self.env_path, "ENV_PATH", self.env_path, export=True)
        if not os.path.exists(self.env_path): subprocess('touch $(pwd)/.env')
        if not load_dotenv(): logger.critical('Unable to load .env file.')

        self.verbose = parser.parse_known_args(['verbose']).verbose; set_key(self.env_path, "VERBOSE", self.verbose, export=True)
        for dirName in ['data', 'logs', 'resources', 'results', 'results', 'src', 'tools']:
            path = os.path.join(self.base_path, dirName); os.makedirs(path, exist_ok=True)
            setattr(self, f'{dirName}_path', path); set_key(self.env_path, f'{dirName.upper()}_PATH', path, export=True)

        # initParser = parser.add_argument_group(title='init')

class Run(): pass

