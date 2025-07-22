#!/usr/bin/env python3
import os, argparse, subprocess, sys
from pathlib import Path
from loguru import logger
from dotenv import load_dotenv, set_key, get_key

class Setup:
    def setupEnvFile(self):
        base_path = self.base_dir
        env_path = os.path.join(base_path, '.env')
        if not os.path.exists(env_path): subprocess('touch $(pwd)/.env')
        if not load_dotenv(): raise Exception('.env not found..')

        # set paths in .env
        set_key(env_path, "BASE_PATH", base_path, export=True)
        set_key(env_path, "ENV_PATH", env_path, export=True)
        subdirs = ["config", "data", "docs", "logs", "resources", "results", "src", "tools"]
        for dir in subdirs: set_key(env_path, f'{dir.upper()}_PATH', os.path.join(base_path, dir.lower()), export=True)

        # set empty API secret variables in .env
        accepted_keys = ['huggingface', 'openrouter', 'openai', 'gemini', 'anthropic']
        for key in accepted_keys: set_key(env_path, f'{key.upper()}_API_KEY', value_to_set='', export=True)

        # add verbosity settings
        set_key(env_path, "VERBOSITY_DEBUG",value_to_set= "", export=True)
        set_key(env_path, "VERBOSITY_MODE", value_to_set="", export=True)

        return load_dotenv()

    def setupLogger(self):
        logDir = self.log_dir
        logger.remove()

        if os.getenv("LEVEL").lower()=='true': LEVEL = "DEBUG"
        elif os.getenv("LEVEL").lower()=='false': LEVEL = "WARNING"
        else: LEVEL = os.getenv("LEVEL").upper()

        MODE = os.getenv("MODE", "development")

        os.makedirs(logDir, exist_ok=True)
        os.makedirs(os.path.join(logDir, ))

        logger.add(
            sys.stderr,
            format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
            colorize=True,
            level=LEVEL
        )
        logger.add(
            os.path.join(logDir, f'app_log.json'), serialize=True
        )
        logger.add(
            os.path.join(logDir, f"{MODE}_debug_{{time}}.log"),
            level="DEBUG",
            compression="zip"
        )
        logger.add(
            os.path.join(logDir, f'{MODE}_info_{{time}}.log'),
            level="INFO"
        )
        logger.add(
            os.path.join(logDir, f'{MODE}_warning_{{time}}.log'),
            level="WARNING"
        )
        logger.add(
            os.path.join(logDir, f"{MODE}_error_{{time}}.log"),
            level="ERROR",
            backtrace=True,
            diagnose=True
        )
        logger.add(
            os.path.join(logDir, f'{MODE}_critical_{{time}}.log'),


        )
        logger.configure(extra={"mode": MODE, "level": LEVEL})
        return logger

    def format(): pass

    def __init__(self, baseDir:str, logDir:str):
        if os.path.exists(baseDir): self.base_dir = baseDir
        if os.path.exists(logDir): self.log_dir = logDir

class Run:
    def __init__(self):
        self.base_path = os.getcwd()
        self.app_setup = Setup(self.base_path)

    # def init(): pass

class Verification:
    def verifyDataPath(dataPath):
        dataPath = str(dataPath)
        configPath = os.path.join(os.getcwd(), 'config')

        if not os.path.exists(configPath):
            logger.error('config/ directory path not found..')
            return False, ''

        for file in os.listdir(configPath):
            fileName = str(os.fsdecode(file))
            if dataPath=="":
                if fileName.endswith(['.csv', '.xlsx']):
                    filePath = os.path.join(configPath, fileName)
                    import pandas as pd
                    dataCol = pd.read_csv(filePath, nrows=0).columns if fileName.endswith('.csv') else pd.read_excel(filePath, nrows=0).columns

            elif dataPath.endswith(['.csv', '.xlsx']):
                if fileName==str(args.data_path): return True, os.path.join(configPath, fileName)
                else:
                    logger.error(f'No files match {dataPath} within config/')
                    return False, ''
            else:
                logger.error(f'Invalid --data_path entry:\n{dataPath}')
                return False, ''

    def setRenames(alter: list[dict]): pass
    # take in list of dicts (single key/value pair per dict) and rename values within the data_path input

def parseArgs():
    parser = argparse.ArgumentParser(
        prog='llm4pgcc', usage="llm4pgcc [-h, --help] [-v, --verbose] COMMAND",
        description='LLM4PGCC CLI TOOL'
    )
    parser.add_argument('-f', '--full', type="store_true", help="")
    parser.add_argument('-s', '--split', type="store_true", help="")
    parser.add_argument('-r', '--ready', type="store_true", help="")
    parser.add_argument('-R', '--rename_columns', nargs="+", action='append')
    parser.add_argument('-R', '--rename_columns'
        nargs="+", action='append', default=None, type=dict,
        help=""
    )
    parser.add_argument('-d', '--data_path',
        type=str, default="",
        help=''
    )
    parser.add_argument('-v', '--verbose',
        default=True, choices=[
            True, False, 'TRACE', 'DEBUG', 'INFO', 'SUCCESS', 'WARNING',
            'ERROR', 'CRITICAL'
        ], help=''
    )
    args, _ = parser.parse_known_args(['full', 'split', 'ready', 'rename_columns', 'data_path', 'verbose'])

    if not any(args.full, args.split, args.ready): args.full = True

    return parser, {'data_path': args.data_path, 'process': args.process, 'verbose': args.verbose}

if __name__ == "__main__":
    parser, tempArgs = parseArgs()
