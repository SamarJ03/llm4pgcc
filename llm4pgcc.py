#!/usr/bin/env python3
import os, argparse, subprocess, sys
from pathlib import Path
import pandas as pd
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

    def __init__(self, baseDir:str, logDir:str):
        if os.path.exists(baseDir): self.base_dir = baseDir
        if os.path.exists(logDir): self.log_dir = logDir

class Verification:
    def verifyDataPath(tempArgs):
        dataPath = str(tempArgs['data_path'])
        configPath = os.path.join(os.getcwd(), 'config')

        if not os.path.exists(configPath):
            logger.error('config/ directory path not found..')
            return False, ''

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
                    return False, ''
            else:
                logger.error(f'Invalid --data_path entry:\n{dataPath}')
                return False, ''

    def verifyFileType(fileName: str, configPath: str):
        filePath = os.path.join(configPath, fileName)
        dataCol = pd.read_csv(filePath, nrows=0).columns if fileName.endswith('.csv') else pd.read_excel(filePath, nrows=0).columns


    def setRenames(alter: list[dict]): pass
    # take in list of dicts (single key/value pair per dict) and rename values within the data_path input

    def verifyDataMode(): pass
    # takes in args[full, split and ready]. Ensures only one mode is selected (defaults to 'ready').
    # it checks formatting to make sure "ready" is a dir path with proper subdirs/data types,
    #  and full/split are file paths with proper data type (delegate to self.verifyDataPath()??)

def parseInitArgs():
    parser = argparse.ArgumentParser(
        prog='llm4pgcc', usage="llm4pgcc [-h, --help] [-v, --verbose] COMMAND",
        description='LLM4PGCC CLI TOOL'
    )

    dataType = parser.add_mutually_exclusive_group(required=True)
    dataType.add_argument('-f', '--full', action='store_true')
    dataType.add_argument('-s', '--split', action='store_true')
    dataType.add_argument('-r', '--ready', action='store_true')
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
    return parser, vars(args)

class Setup:
    def __init__(self):
        # set class variables
        self.base_path = os.getcwd()
        for dirName in ['data', 'logs', 'resources', 'results']:
            path = os.path.join(self.base_path, dirName); os.makedirs(path, exist_ok=True)
            setattr(self, f'{dirName}_path', path)

        # parse init arguments from argparse
        self.verbose, self.data_path, self.data_mode, self.set_labels = None, None, None, None
        parser = self.parseInitArgs()

        # configure .env file
        self.env_path = os.path.join(self.base_path, '.env')
        if not os.path.exists(self.env_path): subprocess('touch $(pwd)/.env')
        if not load_dotenv(): pass

        self.setupEnvFile()
        self.setupLogger()

    def setupEnvFile(self):
        basePath = self.base_path; envPath = self.env_path
        if not os.path.exists(envPath): subprocess('touch $(pwd)/.env')
        if not load_dotenv(): raise Exception('.env not found..')

        # set paths in .env
        set_key(envPath, "BASE_PATH", basePath, export=True)
        set_key(envPath, "ENV_PATH", envPath, export=True)
        subdirs = ["config", "data", "docs", "logs", "resources", "results", "src", "tools"]
        for dir in subdirs: set_key(envPath, f'{dir.upper()}_PATH', os.path.join(basePath, dir.lower()), export=True)

        # set empty API secret variables in .env
        accepted_keys = ['huggingface', 'openrouter', 'openai', 'gemini', 'anthropic']
        for key in accepted_keys: set_key(envPath, f'{key.upper()}_API_KEY', value_to_set='', export=True)

        # add verbosity settings
        set_key(envPath, "VERBOSITY_DEBUG",value_to_set= "", export=True)
        set_key(envPath, "VERBOSITY_MODE", value_to_set="", export=True)

        return load_dotenv()

    def setupLogger(self):
        logDir = self.log_path
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
        logger.add(os.path.join(logDir, f'app_log.json'), serialize=True)
        logger.add(os.path.join(logDir, f"{MODE}_debug_{{time}}.log"),level="DEBUG",compression="zip")
        logger.add(os.path.join(logDir, f'{MODE}_info_{{time}}.log'),level="INFO")
        logger.add(os.path.join(logDir, f'{MODE}_warning_{{time}}.log'),level="WARNING")
        logger.add(os.path.join(logDir, f"{MODE}_error_{{time}}.log"),level="ERROR",backtrace=True,diagnose=True)
        logger.add(os.path.join(logDir, f'{MODE}_critical_{{time}}.log'),)
        logger.configure(extra={"mode": MODE, "level": LEVEL})
        return logger

    def verifyDataPath(tempArgs):
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

    def parseInitArgs(self):
        parser = argparse.ArgumentParser(
            prog='l4p', usage="llm4pgcc [-h, --help] [-v, --verbose] COMMAND",
            description='LLM4PGCC CLI TOOL'
        )
        parser.add_argument('-v', '--verbose',
            default=True, choices=[
                True, False, 'TRACE', 'DEBUG', 'INFO', 'SUCCESS', 'WARNING', 'ERROR', 'CRITICAL'
            ], help=''
        )

        initParser = parser.add_argument_group(title='init')
        dataType = initParser.add_mutually_exclusive_group(required=True)
        dataType.add_argument('-f', '--full', action='store_true')
        dataType.add_argument('-s', '--split', action='store_true')
        dataType.add_argument('-r', '--ready', action='store_true')
        initParser.add_argument('-d', '--data_path', type=str, default="", help='')
        initParser.add_argument('-S', '--set_labels', type=dict[str:str], required=False, help="")

        args, _ = parser.parse_known_args(['full', 'split', 'ready', 'set_labels', 'data_path', 'verbose'])
        if args.full: self.data_mode = 'full'
        elif args.split: self.data_mode = 'split'
        elif args.ready: self.data_mode = 'ready'
        if args.set_labels is not None: self.set_labels = args.set_labels
        self.data_path = args.data_path
        self.verbose = args.verbose

        return parser

class Run(): pass

