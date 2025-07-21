import os
from pathlib import Path
import argparse
from loguru import logger
from dotenv import load_dotenv, set_key, get_key

def setupEnv():
    base_path = os.getcwd()
    env_path = os.path.join(base_path, '.env') # ~/.env
    set_key(env_path, "BASE_PATH", base_path) # <base directory>

    config_var = set_key(env_path, "CONFIG_PATH", os.path.join(base_path, 'config')) # ~/config
    if config_var[0]: os.makedirs(config_var[2], exist_ok=True)

    build_var = set_key(env_path, "BUILD_PATH", os.path.join(base_path, 'config', 'build')) # ~/config/build
    if build_var[0]: os.makedirs(build_var[2], exist_ok=True)

    data_var = set_key(env_path, "DATA_PATH", os.path.join(base_path, 'data')) # ~/data
    if data_var[0]: os.makedirs(data_var[2], exist_ok=True)

    docs_var = set_key(env_path, "DOCS_PATH", os.path.join(base_path, 'docs')) # ~/docs
    if docs_var[0]: os.makedirs(docs_var[2], exist_ok=True)

    log_var = set_key(env_path, "LOGS_PATH", os.path.join(base_path, 'logs')) # ~/logs
    if log_var[0]: os.makedirs(log_var[2], exist_ok=True)

    resource_var = set_key(env_path, "RESOURCES_PATH", os.path.join(base_path, 'resources')) # ~/resources
    if resource_var[0]: os.makedirs(resource_var[2], exist_ok=True)

    result_var = set_key(env_path, "RESULTS_PATH", os.path.join(base_path, 'results')) # ~/results
    if result_var[0]: os.makedirs(result_var[2], exist_ok=True)

    src_var = set_key(env_path, "SRC_PATH", os.path.join(base_path, 'src')) # ~/src
    if src_var[0]: os.makedirs(src_var[2], exist_ok=True)

    tools_var = set_key(env_path, "TOOLS_PATH", os.path.join(base_path, 'tools')) # ~/tools
    if tools_var[0]: os.makedirs(tools_var[2], exist_ok=True)
    
    return load_dotenv(env_path)

def setupLogger(verbosity):
    if isinstance(verbosity, bool): verbosity = "TRACE" if verbosity else "ERROR"
    elif isinstance(verbosity, int): verbosity = min(max(0, verbosity), 5)
    elif isinstance(verbosity, str): verbosity = str(verbosity).upper()

    logger.remove(0)
    logger.remove()
    logger.add(
        get_key("LOGS_PATH"), level=verbosity, colorize=True,
        format="<red>[{level}]</red> Message : <green>{message}</green> @ {time}"
    )
    # logger.add(get_key("LOGS_PATH"), format = "<red>[{level}]</red> Message : <green>{message}</green> @ {time}", colorize=True)

def setupConda():
    # meow
    return condaCheck

def process(): pass

def getParser(): 
    parser = argparse.ArgumentParser(
        prog='llm4pgcc', usage="llm4pgcc [-h, --help] [-v, --verbose] COMMAND",
        description='LLM4PGCC CLI TOOL'
    )
    parser.add_argument('--data_path', type=str, default=os.path.join(get_key("CONFIG_PATH"), 'config_library.xlsx'), help='')
    parser.add_argument('--verbose',
        default=True, choices=[
            True, False, 'TRACE', 'DEBUG', 'INFO', 'SUCCESS', 'WARNING',
            'ERROR', 'CRITICAL', list(range(6))
        ], help=''
    )
    parser.add_argument('--process',
        default=True, choices=[True, False, 'full', 'split', 'ready', list(range(2))], help=(
            ''
        )
    )
    args, _ = getParser().parse_known_args(['data_path', 'process', 'verbose'])
    return parser, {'data_path': args.data_path, 'process': args.process, 'verbose': args.verbose}