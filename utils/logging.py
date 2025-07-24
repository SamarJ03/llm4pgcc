import os, sys, argparse
from loguru import logger
from dotenv import load_dotenv
load_dotenv()

class Log:
    def __init__(self, initialParser: argparse.ArgumentParser, verbose=os.getenv("DEFAULT_VERBOSE")):
        logParser = initialParser.add_argument_group(title='--set_log')
        logParser.add_argument('')
        # retention, default_verbose,

    def initLogger(rtn='10 MB'):
        logDir = os.getenv("LOGS_PATH")
        logger.remove()
        if not os.path.exists(logDir): os.makedirs(logDir)

        appLog_path = os.path.join(logDir, 'application_logs'); os.makedirs(appLog_path, exist_ok=True)
        for subDir in ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']:
            path = os.path.join(appLog_path, f'{subDir.lower()}_logs'); os.makedirs(path, exist_ok=True)
            logger.add(
                os.path.join(path, f'application_{subDir}_{{time}}.log'),
                level=subDir, rotation=rtn
            )
        perfLog_path = os.path.join(logDir, 'performace_logs'); os.makedirs(perfLog_path, exist_ok=True)
        for subDir in ['dataPrep', 'createPrompts', 'synthesize', 'inference', 'summarize', 'codeGen']:
            path = os.path.join(perfLog_path, f'{subDir}_logs'); os.makedirs(path, exist_ok=True)
            logger.add(
                os.path.join(path, f'performance_{subDir}_{{time}}.log'),
                level="", rotation=rtn
            )

        return logger

    def setLogger(lvl: str):
        logger.remove()
        if lvl.lower() not in ['debug', 'info', 'warning', 'error', 'critical']:
            logger.warning(f'"{lvl}" is not a valid verbose level. Default set to "INFO".'); lvl = "INFO"
        logger.add(
            sys.stderr,
            format="\n<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
            colorize=True, level=lvl
        )
