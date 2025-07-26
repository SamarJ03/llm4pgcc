import os, sys, time, statistics, resource
from loguru import logger
from dotenv import set_key, load_dotenv

class Log:
    def __init__(self, path:str, verbose:str=os.getenv("VERBOSE")):
        load_dotenv()
        if not os.path.exists(path): logger.warning(f'Unknown directory: {path}')
        if verbose.lower() not in ['debug', 'info', 'warning', 'error', 'critical']:
            logger.warning(f'{verbose} is not a valid verbose level. Verbosity set to "INFO"..')
            verbose = "INFO"
        self.cwd = path
        self.vrb = verbose
        self.baseLog_dir = os.getenv("LOGS_PATH")

    def setupLog(self, lvl:str=None, rtn="10 MB"):
        if not lvl or lvl.lower() not in ['debug', 'info', 'warning', 'error', 'critical']: lvl = self.vrb.upper()
        logger.remove()

        logger.add(
            sys.stderr, level=lvl.upper(), colorize=True,
            format="<lvl>{level} | <g>{extra[fileName]}[<g>{function}] | <b>{time:ddd:D:M:YY h:mm:ss A} | {message}</>\n"
        )
        appLog_path = os.path.join(os.getenv("LOGS_PATH"), 'application.log')
        logger.add(
            appLog_path, level=lvl.upper(), rotation=rtn, catch=True,
            format="{level} | {extra[fileName]}[{function}] | {time:ddd:D:M:YY h:mm:ss A} | {message}"
        )
        errLog_path = os.path.join(os.getenv("LOGS_PATH"), 'error.log')
        logger.add(
            errLog_path, level="ERROR", rotation=rtn, catch=True,
            format="{extra[fileName]}[{function} @ ln:{line}] | {time:ddd:D:M:YY h:mm:ss A} | {message}"
        )


        logger.add(os.path.join(os.getenv("LOGS_PATH"), "filtered_data.log"), serialize=True)

def memTrack():
    resource.getrusage(resource.RUSAGE_SELF)

    # def setBase(self, cwd:str, appLog=False, perfLog=False):
    #     logger.remove()
    #     logger.add(
    #         sys.stderr,
    #         format="\n<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level}</level> | "
    #         "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
    #         colorize=True, level=self.vrb
    #     )

    # def setAppLog(self, rtn='10 MB'):
    #     appLog_path = os.path.join(self.baseLog_dir, 'application_logs'); os.makedirs(appLog_path, exist_ok=True)
    #     for subDir in ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']:
    #         path = os.path.join(appLog_path, f'{subDir.upper()}_application_{{time}}.log')
    #         logger.add(
    #             path, level=subDir, rotation=rtn,
    #             format=''
    #         )

    # def setAppLog(self, rtn='10 MB'):
    #     appLog_path = os.path.join(self.baseLog_dir, 'application_logs'); os.makedirs(appLog_path, exist_ok=True)
    #     for subDir in ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']:
    #         path = os.path.join(appLog_path, f'{subDir.upper()}_application_{{time}}.log')
    #         logger.add(
    #             path, level=subDir, rotation=rtn,
    #             format=''
    #         )

        # if appLog:
        #     appLog_path = os.path.join(self.baseLog_dir, 'application_logs')
        #     os.makedirs(appLog_path, exist_ok=True)
        #     for subDir in ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']:
        #         path = os.path.join(appLog_path, f'{subDir.lower()}_logs')
        #         os.makedirs(path, exist_ok=True)
        #         logger.add(
        #             os.path.join(path, f'application_{subDir}_{{time}}.log'),
        #             level=subDir, rotation='10 MB',

        #         )

        # if perfLog:
        #     perfLog_path = os.path.join(self.baseLog_dir, 'performace_logs')
        #     os.makedirs(perfLog_path, exist_ok=True)
        #     for subDir in ['dataPrep', 'createPrompts', 'synthesize', 'inference', 'summarize', 'codeGen']:
        #         path = os.path.join(perfLog_path, f'{subDir}_logs')
        #         os.makedirs(path, exist_ok=True)
        #         logger.add(
        #             os.path.join(path, f'performance_{subDir}_{{time}}.log'),
        #             level="", rotation='10 MB',
        #             format=''
        #         )

